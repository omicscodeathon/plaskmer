import streamlit as st
import pandas as pd
import os
import time
from pathlib import Path
import regex
from huggingface_hub import HfApi, hf_hub_download
from dotenv import load_dotenv
from harvester import ParallelOmniSystem # YOUR REAL HARVESTER
from Bio import SeqIO
from kmer_tool import count_kmers, get_all_kmers, fasta_to_kmer_df # Untouched K-mer script
from sklearn.metrics.pairwise import cosine_similarity
import numpy as np
import plotly.express as px
from constants import AFRICA
from Bio import Restriction
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import matplotlib.pyplot as plt
from dna_features_viewer import GraphicFeature, CircularGraphicRecord
from gene_annotator import fetch_official_genes
from orf_blaster import identify_orf_via_blast
import pyarrow.parquet as pq
import config

# Load local .env file
load_dotenv()


# --- PATH CONFIGURATION (PARQUET ENGINE) ---
# --- PATH CONFIGURATION (PARQUET ENGINE) ---
SCRIPT_DIR = Path(__file__).parent.resolve()

def find_database_file(filename):
    # 1. Check local scripts folder
    path_same_folder = SCRIPT_DIR / filename
    if path_same_folder.exists():
        return path_same_folder
    
    # 2. Check local data folder
    path_github = SCRIPT_DIR.parent / "data" / filename
    if path_github.exists():
        return path_github
        
    # 3. CLOUD FALLBACK: Download directly from Hugging Face
    import streamlit as st
    try:
        hf_token = os.environ.get("HF_TOKEN")
        repo_id = getattr(config, 'HF_REPO_ID', "Jeffiq/Plaskmer") 
        
        # Try looking in the data/ folder first...
        try:
            cloud_path = hf_hub_download(
                repo_id=repo_id, 
                filename=f"data/{filename}", 
                repo_type="dataset",
                token=hf_token
            )
        except Exception:
            # If not in data/, try looking in the main root folder...
            cloud_path = hf_hub_download(
                repo_id=repo_id, 
                filename=filename, 
                repo_type="dataset",
                token=hf_token
            )
            
        return Path(cloud_path)
    except Exception as e:
        # Stop hiding the error! Show exactly what Hugging Face is complaining about.
        st.error(f"⚠️ Cloud Download Failed for '{filename}': {e}")
        return None

# The single, high-performance database file
PARQUET_FILE = find_database_file(config.MASTER_PARQUET)

# 
# --- HUGGING FACE HELPER ---
def push_to_huggingface(local_filepath, repo_filepath):
    hf_token = os.environ.get("HF_TOKEN")
    if not hf_token:
        st.error("⚠️ HF_TOKEN not found in .env or secrets. Cannot push to cloud.")
        return False
    try:
        api = HfApi()
        repo_id = getattr(config, 'HF_REPO_ID', "Jeffiq/Plaskmer") # Uses config or defaults
        api.upload_file(
            path_or_fileobj=str(local_filepath),
            path_in_repo=f"data/{repo_filepath}",
            repo_id=repo_id, 
            repo_type="dataset",
            token=hf_token,
            commit_message="Auto-update database via Streamlit Harvester"
        )
        return True
    except Exception as e:
        st.error(f"❌ HF Upload Failed: {e}")
        return False

# --- VIRTUAL PCR HELPERS ---
def reverse_complement(seq):
    base_dict = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
    return "".join([base_dict.get(i.upper(), "N") for i in seq[::-1]])

def find_matches(chrom_seq, f_primer, r_primer, mismatch, lower_limit, upper_limit):
    proper_matches = []
    f_primer_query = f"(?:{f_primer}){{s<={mismatch}}}"
    r_primer_query = f"(?:{r_primer}){{s<={mismatch}}}"
    f_matches = regex.finditer(f_primer_query, chrom_seq)
    r_matches = regex.finditer(r_primer_query, chrom_seq)
    f_starts = [i.start() for i in f_matches]
    r_stops = [i.end() for i in r_matches]
    for i in f_starts:
        for j in r_stops:
            size = j - i
            if lower_limit <= size <= upper_limit:
                proper_matches.append({"start": i + 1, "end": j, "size": size, "amplicon": chrom_seq[i:j]})
    return proper_matches

# --- RAM-SAFE SEQUENCE LOADER ---
def generate_chromosome_dict(parquet_path, filter_plasmids=False):
    """Loads sequences directly from Parquet. Uses filtering to save RAM."""
    chrom_dict = {}
    if parquet_path and parquet_path.exists():
        try:
            # Load ONLY the specific columns we need
            df = pd.read_parquet(parquet_path, columns=["Accession", "Sequence", "Type"])
            
            # Filter out SRA (which has blank sequences) so they don't crash the tools
            df = df[df['Sequence'].astype(str).str.len() > 0]
            
            # Optional: Further reduce RAM by only loading Plasmids for Tab 7
            if filter_plasmids:
                df = df[df['Type'].str.contains("Plasmid", case=False, na=False)]
                
            for _, row in df.iterrows():
                chrom_dict[row['Accession']] = row['Sequence']
        except Exception as e:
            st.error(f"Error loading sequences: {e}")
            
    return chrom_dict

# --- INCREMENTAL K-MER HELPER (BATCHED PARQUET MAGIC) ---
def update_kmer_database(parquet_path, k=6):
    csv_filename = f"kmer_{k}mer_vectors.csv"
    csv_path = find_database_file(csv_filename) or (SCRIPT_DIR / csv_filename)
    
    existing_ids = set()
    if csv_path.exists():
        # Load only the ID column to save RAM
        existing_df = pd.read_csv(csv_path, usecols=['sequence_id'])
        existing_ids = set(existing_df['sequence_id'].astype(str))
    
    if not parquet_path or not parquet_path.exists():
        return csv_path, 0
        
    new_records_count = 0
    parquet_file = pq.ParquetFile(parquet_path)
    
    # MAGIC TRICK: Process in batches of 100 to absolutely protect desktop RAM
    for batch in parquet_file.iter_batches(batch_size=100):
        df_batch = batch.to_pandas()
        
        # 1. Filter: Only get rows with sequences that aren't in the CSV yet
        new_seqs = df_batch[
            (~df_batch['Accession'].astype(str).isin(existing_ids)) & 
            (df_batch['Sequence'].astype(str).str.len() > 0)
        ]
        
        if not new_seqs.empty:
            # 2. Write Temporary FASTA
            temp_fasta = SCRIPT_DIR / "temp_batch.fasta"
            with open(temp_fasta, "w") as f:
                for _, row in new_seqs.iterrows():
                    f.write(f">{row['Accession']}\n{row['Sequence']}\n")
            
            # 3. Pass to untouched kmer_tool
            new_kmer_df = fasta_to_kmer_df(str(temp_fasta), k=k)
            
            # 4. Append to the master CSV
            if csv_path.exists():
                new_kmer_df.to_csv(csv_path, mode='a', header=False)
            else:
                new_kmer_df.to_csv(csv_path)
                
            new_records_count += len(new_kmer_df)
            
            # 5. Delete Temporary FASTA (Clean up the evidence!)
            os.remove(temp_fasta)
            
    return csv_path, new_records_count

# --- MATH HELPER ---
def find_closest_matches(query_seq, k=6, top_n=5):
    csv_filename = f"kmer_{k}mer_vectors.csv"
    csv_path = find_database_file(csv_filename) or (SCRIPT_DIR / csv_filename)
    
    # 1. CHECK AND CREATE IF MISSING
    if not csv_path.exists():
        if PARQUET_FILE and PARQUET_FILE.exists():
            with st.spinner(f"🚨 {k}-mer index not found. Generating safely from Parquet Database in batches..."):
                # Run the batched update
                csv_path, num_new = update_kmer_database(PARQUET_FILE, k=k)
                
                if num_new > 0:
                    st.success(f"✅ Created {csv_filename} locally with {num_new} sequences.")
                    
                    with st.spinner("☁️ Pushing new index to Hugging Face repository..."):
                        pushed = push_to_huggingface(csv_path, csv_path.name)
                        if pushed:
                            st.success(f"☁️ Upload Successful: {csv_filename} is now live on Hugging Face!")
                        else:
                            st.error("❌ Cloud Sync Failed. File created locally but not pushed.")
                else:
                    return None, "K-mer calculation returned 0 sequences. Are there valid assembled sequences in the database?"
        else:
            return None, "Database missing. Please run the Harvester (Tab 2) first."

    # 3. PROCEED WITH THE SEARCH
    try:
        db_df = pd.read_csv(csv_path, index_col='sequence_id')
        
        all_kmers = get_all_kmers(k)
        counts = count_kmers(query_seq, k)
        total = sum(counts.values())
        if total == 0: return None, "Query sequence contains no valid k-mers."
        
        query_vector = np.array([counts.get(k_str, 0) / total for k_str in all_kmers]).reshape(1, -1)
        db_vectors = db_df.values

        similarities = cosine_similarity(query_vector, db_vectors)[0]
        
        results_df = pd.DataFrame({
            'Sequence ID': db_df.index,
            'Similarity Score': similarities
        }).sort_values(by='Similarity Score', ascending=False).head(top_n)
        
        return results_df, None
        
    except Exception as e:
        return None, f"System Error: {str(e)}"

    

# ==========================================
# PAGE CONFIGURATION & MAGNIFICENT UI
# ==========================================
st.set_page_config(page_title="The PlasKmer Studio", page_icon="🧬", layout="wide")

st.markdown("""
<style>
    #MainMenu {visibility: hidden;}
    footer {visibility: hidden;}
    header {visibility: hidden;}
    .stApp { background-color: #f4f7f6; }
    .stTabs [data-baseweb="tab-list"] { gap: 8px; background-color: transparent; }
    .stTabs [data-baseweb="tab"] {
        height: 50px; background-color: #ffffff;
        border-radius: 8px 8px 0px 0px; padding: 10px 24px;
        font-weight: 600; color: #4a5568;
        border: 1px solid #e2e8f0; border-bottom: none;
    }
    .stTabs [aria-selected="true"] { background-color: #2e66ff !important; color: white !important; border-color: #2e66ff; }
    .The-subheader { color: #10b981; font-weight: 700; border-bottom: 2px solid #e2e8f0; padding-bottom: 10px; margin-bottom: 20px; }
</style>
""", unsafe_allow_html=True)

html_banner = """
<div style="background: linear-gradient(135deg, #1e3a8a 0%, #3b82f6 100%); padding: 30px; border-radius: 12px; color: white; text-align: center; margin-bottom: 20px; box-shadow: 0 10px 15px -3px rgba(0, 0, 0, 0.1);">
    <h1 style="margin: 0; font-size: 2.5em; font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;">🧬The PlasKmer Studio</h1>
    <p style="margin: 10px 0 0 0; font-size: 1.2em; opacity: 0.9;">The one-stop high-performance portal for African Plasmid and mRNA informatics.</p>
</div>
"""
st.markdown(html_banner, unsafe_allow_html=True)

# --- TAB NAVIGATION ---
tab1, tab2, tab3, tab4, tab5, tab6, tab7, tab8 = st.tabs([
    "🔍 1. K-mer Search", 
    "🚜 2. Data Harvester", 
    "🎯 3. Virtual PCR", 
    "📊 4. Analytics", 
    "🌍 5. Global Stats",
    "🚧 6. Future Tools",
    "🧬 7. Plasmid Hub",
    "🛡️ 8. Gene Vault" 
])


# ==========================================
# TAB 1: K-MER SEARCH 
# ==========================================
with tab1:
    st.markdown("<h3 class='The-subheader'>Pattern-Based Retrieval & Closest Match</h3>", unsafe_allow_html=True)
    
    with st.container(border=True):
        col1, col2, col3 = st.columns([2, 1, 1])
        with col1:
            st.markdown("**Sequence Input**")
            query_input = st.text_area(
                "Paste Accession Number, FASTA, or Raw Sequence (Supports ambiguous codes: S, M, N, R...):", 
                placeholder="e.g. CP000921\nOR\n>MySequence\nATGCGTMSN...",
                height=150
            )
            uploaded_kmer_file = st.file_uploader("Or upload FASTA/CSV (Max 200MB)", type=["fasta", "fa", "csv"])
            
        with col2:
            st.markdown("**Search Parameters**")
            kmer_length = st.slider("K-mer Length (k)", min_value=3, max_value=12, value=6)
            match_type = st.radio("Match Type:", ["Exact", "Approximate"])
            
        with col3:
            st.markdown("**Database Scope**")
            db_scope = st.selectbox("Search In:", ["Both (Plasmid & mRNA)", "Plasmid Only", "mRNA Only"])
        
        search_triggered = st.button("🚀 Find Closest Match", type="primary", use_container_width=True)

    if search_triggered:
        if not query_input:
            st.warning("Please paste a sequence first.")
        else:
            with st.spinner(f"Calculating {kmer_length}-mer similarities..."):
                clean_seq = "".join([line for line in query_input.splitlines() if not line.startswith(">")])
                results, error = find_closest_matches(clean_seq, k=kmer_length)
                
                if error:
                    st.error(error)
                else:
                    st.success(f"Top {len(results)} Matches Found!")
                    best_match = results.iloc[0]
                    st.metric("Top Match Score", f"{best_match['Similarity Score']:.4f}")
                    st.table(results)


# ==========================================
# TAB 2: DATA HARVESTER
# ==========================================
with tab2:
    st.markdown("<h3 class='The-subheader'>Database Management & Processing</h3>", unsafe_allow_html=True)
    colA, colB = st.columns([1.5, 1])
    with colA:
        # --- HARVESTER BLOCK ---
        with st.container(border=True):
            st.markdown("#### 🚜 Run Parallel Omni Harvester")
            target_org = st.text_input("Target Organism:", placeholder="e.g., Vibrio cholerae")
            
            st.markdown("**1. Select Databases:**")
            db_col1, db_col2, db_col3, db_col4 = st.columns(4)
            with db_col1: use_bp = st.checkbox("BioProject")
            with db_col2: use_bs = st.checkbox("BioSample")
            with db_col3: use_sra = st.checkbox("SRA")
            with db_col4: use_nuc = st.checkbox("Nucleotide", value=True)
            
            st.markdown("**2. Select Types:**")
            type_col1, type_col2, type_col3 = st.columns(3)
            with type_col1: inc_plasmids = st.checkbox("Plasmids", value=True)
            with type_col2: inc_wgs = st.checkbox("WGS")
            with type_col3: inc_genes = st.checkbox("Specific Genes")
            
            target_goal = st.number_input("Records per country:", min_value=1, value=10)
            run_kmer_calc = st.checkbox("🧬 Auto-Calculate K-mers after harvest?", value=True)
            k_size = st.number_input("K-mer Size (k):", min_value=2, max_value=12, value=6)
            
            # This is automatically ON because of value=True
            upload_hf = st.toggle("☁️ Push to Hugging Face after harvest?", value=True)
            
            if st.button("🚜 Start Harvester Pipeline", type="primary", use_container_width=True):
                if not target_org:
                    st.error("Please enter a Target Organism.")
                else:
                    selected_dbs = []
                    if use_bp: selected_dbs.append("bioproject")
                    if use_bs: selected_dbs.append("biosample")
                    if use_sra: selected_dbs.append("sra")
                    if use_nuc: selected_dbs.append("nucleotide")
                    
                    types = []
                    if inc_plasmids: types.append("plasmid")
                    if inc_wgs: types.append("WGS")
                    if inc_genes: types.append("gene")
                    
                    with st.spinner("Harvesting via Parallel OmniSystem..."):
                        bot = ParallelOmniSystem(target_org=target_org, selected_dbs=selected_dbs, types=types, target_goal=target_goal, push_to_hf=False)
                        bot.run_all_parallel()
                    st.success("✅ Local Harvest Complete!")
                                                                   
                    # Trigger Incremental K-mers using the new Batched logic
                    if run_kmer_calc and PARQUET_FILE and PARQUET_FILE.exists():
                        with st.spinner(f"🧬 Incrementally calculating {k_size}-mers for new sequences..."):
                            kmer_csv_path, num_new = update_kmer_database(PARQUET_FILE, k=int(k_size))
                            if num_new > 0:
                                st.success(f"✅ Added {num_new} new sequences to {kmer_csv_path.name}!")
                            else:
                                st.info(f"⚡ No new sequences needed k-mer calculation. Database is up to date.")
                    
                    # Upload updated files
                    if upload_hf:
                        with st.spinner("☁️ Pushing updated files to Hugging Face..."):
                            if PARQUET_FILE and PARQUET_FILE.exists():
                                push_to_huggingface(PARQUET_FILE, PARQUET_FILE.name)
                            if run_kmer_calc and 'kmer_csv_path' in locals() and kmer_csv_path.exists():
                                push_to_huggingface(kmer_csv_path, kmer_csv_path.name)
                            st.success("✅ Cloud Database Synced!")

        # --- MANUAL CLOUD SYNC BLOCK ---
        with st.container(border=True):
            st.markdown("#### ☁️ Manual Cloud Sync")
            st.caption("Push your existing local database directly to Hugging Face without running a new harvest.")
            
            if st.button("⬆️ Push Existing Database to Hugging Face", use_container_width=True):
                import os
                from huggingface_hub import HfApi
                import config
                
                parquet_file = config.MASTER_PARQUET
                repo_id = getattr(config, 'HF_REPO_ID', None)
                
                if not repo_id or repo_id == "YourUsername/PlasKmer-Database":
                    st.error("⚠️ Please update HF_REPO_ID in your config.py first!")
                elif not os.path.exists(parquet_file):
                    st.error(f"❌ Cannot find '{parquet_file}' locally. Run the Harvester at least once!")
                else:
                    with st.spinner(f"Uploading {parquet_file} to {repo_id}..."):
                        try:
                            api = HfApi() # Uses the token from your .env automatically
                            api.upload_file(
                                path_or_fileobj=parquet_file,
                                path_in_repo=f"data/{parquet_file}",
                                repo_id=repo_id,
                                repo_type="dataset"
                            )
                            st.success(f"✅ Successfully synced local database to {repo_id}!")
                            st.balloons()
                        except Exception as e:
                            st.error(f"❌ Upload failed: {e}")

                    

# ==========================================
# TAB 3: VIRTUAL PCR 
# ==========================================
with tab3:
    st.markdown("<h3 class='The-subheader'>In Silico Amplification</h3>", unsafe_allow_html=True)
    with st.container(border=True):
        pcr_col1, pcr_col2 = st.columns(2)
        with pcr_col1: fwd_primer = st.text_input("Forward Primer:", key="fwd")
        with pcr_col2: rev_primer = st.text_input("Reverse Primer:", key="rev")
            
        set_col1, set_col2, set_col3 = st.columns(3)
        with set_col1: pcr_mismatches = st.number_input("Mismatches", value=0)
        with set_col2: min_len = st.number_input("Min Size", value=50)
        with set_col3: max_len = st.number_input("Max Size", value=5000)
            
        if st.button("🔬 Run Virtual PCR", type="primary", use_container_width=True):
            if not fwd_primer or not rev_primer:
                st.error("Enter primers.")
            elif not PARQUET_FILE or not PARQUET_FILE.exists():
                st.error("Parquet database not found!")
            else:
                with st.spinner("🧬 Running Python PCR..."):
                    chrom_dict = generate_chromosome_dict(PARQUET_FILE)
                    all_results = []
                    for name, seq in chrom_dict.items():
                        for match in find_matches(seq, fwd_primer.upper(), reverse_complement(rev_primer.upper()), pcr_mismatches, min_len, max_len):
                            all_results.append({"Sequence Name": name, "Pattern": "Fwd Strand", "Position": f"{match['start']} - {match['end']}", "BP Size": match['size'], "Sequence": match['amplicon']})
                        for match in find_matches(seq, rev_primer.upper(), reverse_complement(fwd_primer.upper()), pcr_mismatches, min_len, max_len):
                            all_results.append({"Sequence Name": name, "Pattern": "Rev Strand", "Position": f"{match['start']} - {match['end']}", "BP Size": match['size'], "Sequence": match['amplicon']})

                if all_results:
                    st.success(f"Found {len(all_results)} amplicons.")
                    st.dataframe(pd.DataFrame(all_results), use_container_width=True)
                else:
                    st.warning("No amplicons found.")

# ==========================================
# TAB 4: DATABASE INSIGHTS & ANALYTICS
# ==========================================
with tab4:
    st.markdown("<h3 class='The-subheader'>Database Insights & Analytics</h3>", unsafe_allow_html=True)
    
    import config
    import os
    import pandas as pd
    from huggingface_hub import hf_hub_download
    import analytics_engine 

    # We will force a check for the file right here
    if PARQUET_FILE and PARQUET_FILE.exists():
        df_raw = pd.read_parquet(PARQUET_FILE)
    else:
        st.warning("⚠️ Local file missing. Attempting direct cloud connection to Hugging Face...")
        try:
            hf_token = os.environ.get("HF_TOKEN")
            repo_id = getattr(config, 'HF_REPO_ID', "Jeffiq/Plaskmer") # Check your capitals here!
            
            # Force the download directly inside Tab 4
            cloud_path = hf_hub_download(
                repo_id=repo_id, 
                filename="data/master_database.parquet", 
                repo_type="dataset",
                token=hf_token
            )
            df_raw = pd.read_parquet(cloud_path)
            st.success("☁️ Successfully connected to Hugging Face Analytics!")
            
        except Exception as e:
            st.error(f"❌ Hugging Face Connection Blocked: {e}")
            df_raw = pd.DataFrame() # Create an empty frame so the app doesn't crash

    # Now proceed with the analytics IF we successfully got data
    if df_raw.empty:
        st.warning("⚠️ Analytics cannot proceed until the cloud connection is fixed.")
    else:
        with st.spinner("Calculating genomic metrics..."):
            df = analytics_engine.prep_analytics_data(df_raw)
            
        # --- GLOBAL METRICS ---
        st.markdown("### 📊 Global Database Metrics")
        
        # PROPERLY INDENTED SECTION STARTS HERE
        m1, m2, m3, m4 = st.columns(4)
        m1.metric("Total Sequences", len(df))
        m2.metric("Unique Organisms", df['Organism'].nunique())
        m3.metric("African Countries", df['Country'].nunique())
        m4.metric("Avg GC Content", f"{df['GC_Content'].mean():.2f}%")
        
        st.markdown("---")
        
        # --- INTERACTIVE ORGANISM ANALYSIS ---
        st.markdown("### 🦠 Deep Organism Analysis")
        
        # Get a list of unique organisms for the dropdown
        available_organisms = df['Organism'].unique().tolist()
        selected_org = st.selectbox("Select an Organism to Analyze:", available_organisms)
        
        if selected_org:
            # 1. Map & Types (Side by side)
            col1, col2 = st.columns([1.5, 1])
            with col1:
                fig_map = analytics_engine.plot_organism_heatmap(df, selected_org)
                st.plotly_chart(fig_map, use_container_width=True)
            with col2:
                fig_types = analytics_engine.plot_type_distribution(df, selected_org)
                st.plotly_chart(fig_types, use_container_width=True)
            
            # 2. Quality Metrics (Lengths & GC)
            fig_quality = analytics_engine.plot_sequence_quality_metrics(df, selected_org)
            st.plotly_chart(fig_quality, use_container_width=True)
        
        st.markdown("---")
        
        # --- CROSS-COUNTRY COMPARISON ---
        st.markdown("### 🌍 Cross-Country Comparison Matrix")
        available_countries = sorted(df['Country'].unique().tolist())
        
        # Default to the first two countries if they exist
        default_countries = available_countries[:2] if len(available_countries) >= 2 else available_countries
        selected_countries = st.multiselect("Select Countries to Compare:", available_countries, default=default_countries)
        
        if selected_countries:
            fig_comparison = analytics_engine.plot_cross_country_comparison(df, selected_countries)
            st.plotly_chart(fig_comparison, use_container_width=True)

        st.markdown("---")
        
        # --- PLASMID & mRNA DEEP DIVE ---
        st.markdown("### 🧬 Plasmid & mRNA Deep Dive")
        st.caption("Isolate and analyze specific genetic elements across the entire database.")
        
        # Let the user toggle between Plasmid and mRNA analysis
        analysis_type = st.radio("Select Sequence Type to Analyze:", ["Plasmid", "mRNA"], horizontal=True)
        
        # Check if there is actually data for this type
        type_exists = df['Type'].str.contains(analysis_type, case=False, na=False).any()
        
        if not type_exists:
            st.info(f"No {analysis_type} records found in the current database.")
        else:
            # 1. Show the specific Heatmap
            fig_type_map = analytics_engine.plot_specific_type_heatmap(df, analysis_type)
            if fig_type_map:
                st.plotly_chart(fig_type_map, use_container_width=True)
            
            # 2. Show the Host Organism Bar and Signature Scatter Plot side-by-side
            fig_bar, fig_scatter = analytics_engine.plot_type_scatter_and_bar(df, analysis_type)
            
            if fig_bar and fig_scatter:
                type_col1, type_col2 = st.columns(2)
                with type_col1:
                    st.plotly_chart(fig_bar, use_container_width=True)
                with type_col2:
                    st.plotly_chart(fig_scatter, use_container_width=True)
        
# =============================================
# TAB 5: GLOBAL STATS
# =============================================
with tab5:
    st.markdown("<h3 class='The-subheader'>Geographic Distribution (Africa)</h3>", unsafe_allow_html=True)
    
    if PARQUET_FILE and PARQUET_FILE.exists():
        df_log = pd.read_parquet(PARQUET_FILE, columns=["Country", "Organism"])
        
        if 'Country' in df_log.columns and 'Organism' in df_log.columns:
            df_log['Clean_Country'] = df_log['Country'].apply(
                lambda x: str(x).split(':')[0].strip() if pd.notnull(x) else ""
            )
            df_africa = df_log[df_log['Clean_Country'].isin(AFRICA)]
            
            if not df_africa.empty:
                country_counts = df_africa.groupby('Clean_Country').size().reset_index(name='Sequence Count')
                org_counts = df_africa.groupby(['Clean_Country', 'Organism']).size().reset_index(name='Org_Count')
                top_orgs = org_counts.sort_values('Org_Count', ascending=False).drop_duplicates('Clean_Country')
                top_orgs = top_orgs.rename(columns={'Organism': 'Top Organism'})
                
                map_data = pd.merge(country_counts, top_orgs[['Clean_Country', 'Top Organism']], on='Clean_Country')
                
                fig_map = px.choropleth(
                    map_data,
                    locations="Clean_Country",
                    locationmode="country names",
                    color="Sequence Count",
                    hover_name="Clean_Country",
                    hover_data={"Clean_Country": False, "Sequence Count": True, "Top Organism": True},
                    scope="africa",
                    color_continuous_scale="Tealgrn",
                    title="Harvested Genomic Data Density"
                )
                
                fig_map.update_layout(
                    margin={"r": 0, "t": 40, "l": 0, "b": 0},
                    geo=dict(showframe=False, showcoastlines=True, projection_type='equirectangular', bgcolor='rgba(0,0,0,0)'),
                    paper_bgcolor='rgba(0,0,0,0)'
                )
                
                st.plotly_chart(fig_map, use_container_width=True)
                
                with st.expander("View Geographic Data Table"):
                    st.dataframe(map_data.sort_values('Sequence Count', ascending=False), use_container_width=True)
            else:
                st.info("No sequence data from African countries currently exists in the database.")
        else:
            st.warning("The current database is missing required columns ('Country' or 'Organism').")
    else:
        st.warning("No database found. Run the Harvester in Tab 2 to generate geographical analytics!")


# =============================================
# TAB 6: PENDING DEVELOPMENT
# =============================================
with tab6:
    st.markdown("<h3 class='The-subheader'>🚧 Future Tools & Modules</h3>", unsafe_allow_html=True)
    with st.container(border=True):
        st.info("**Status:** Pending Development")
        st.markdown("""
        This workspace is reserved for upcoming The PlasKmer integration tools. Future pipelines will include:
        * 🧬 **Advanced Assembly:** De novo plasmid reconstruction modules.
        * 🌳 **Phylogenetics:** Evolutionary mapping of harvested African variants.
        * 📉 **mRNA Folding:** Secondary structure predictions and stability scoring.
        """)

# ==========================================
# TAB 7: PLASMID INTELLIGENCE HUB
# ==========================================
with tab7:
    st.markdown("<h3 class='The-subheader'>Whole Plasmid Analytics</h3>", unsafe_allow_html=True)
    
    if PARQUET_FILE and PARQUET_FILE.exists():
        # Load only plasmids to save memory
        plasmid_dict = generate_chromosome_dict(PARQUET_FILE, filter_plasmids=True)
        
        if not plasmid_dict:
            st.warning("No valid assembled plasmids found in the current Parquet database.")
        else:
            selected_plasmid_id = st.selectbox("Select a Plasmid to Analyze:", list(plasmid_dict.keys()))
            selected_seq = plasmid_dict[selected_plasmid_id]
            seq_length = len(selected_seq)
            
            st.success(f"Loaded valid plasmid: {selected_plasmid_id} ({seq_length} bp)")
            
            colA, colB = st.columns([1.5, 1])
            
            with colA:
                st.markdown("#### 🔄 Circular Visualizer & Cut Sites")
                common_enzymes = ["EcoRI", "BamHI", "HindIII", "SmaI", "XhoI", "NotI", "PstI"]
                selected_enzymes = st.multiselect("Select Nucleases to Map:", common_enzymes, default=["EcoRI", "BamHI"])
                
                if selected_enzymes:
                    with st.spinner("Calculating restriction sites and drawing map..."):
                        rb = Restriction.RestrictionBatch(selected_enzymes)
                        my_seq = Seq(selected_seq)
                        cut_results = rb.search(my_seq)
                        
                        features = []
                        for enz, cuts in cut_results.items():
                            for cut in cuts:
                                features.append(GraphicFeature(start=cut, end=cut+1, strand=1, color="#e74c3c", label=str(enz)))
                        
                        if features:
                            record = CircularGraphicRecord(sequence_length=seq_length, features=features)
                            fig, ax = plt.subplots(figsize=(6, 6))
                            record.plot(ax=ax) 
                            st.pyplot(fig)
                        else:
                            st.info(f"None of the selected enzymes cut this {seq_length} bp sequence.")

            with colB:
                st.markdown("#### 🧬 Nearest Neighbors & Host Range")
                with st.spinner("Running K-mer similarity algorithm..."):
                    sim_results, error = find_closest_matches(selected_seq, k=6, top_n=6)
                    
                    if error:
                        st.error(error)
                    elif sim_results is not None:
                        sim_results = sim_results.iloc[1:] # Drop self match
                        if not sim_results.empty:
                            best_match_id = sim_results.iloc[0]['Sequence ID']
                            best_score = sim_results.iloc[0]['Similarity Score']
                            st.metric("Closest Relative", best_match_id)
                            st.metric("K-mer Similarity", f"{best_score*100:.2f}%")
                            st.markdown("**Other Similar Sequences:**")
                            st.dataframe(sim_results, use_container_width=True)
                        else:
                            st.info("No close neighbors found.")

                st.markdown("#### 🌍 Known Locations")
                df_meta = pd.read_parquet(PARQUET_FILE, columns=["Accession", "Organism", "Country"])
                plasmid_meta = df_meta[df_meta['Accession'] == selected_plasmid_id]
                
                if not plasmid_meta.empty:
                    st.write(f"**Source Organism:** {plasmid_meta.iloc[0].get('Organism', 'Unknown')}")
                    st.write(f"**Country:** {plasmid_meta.iloc[0].get('Country', 'Unknown')}")
                else:
                    st.info("Metadata not found.")
    else:
        st.warning("⚠️ Database missing. Please run the Harvester first.")


# ==========================================
# TAB 8: GENE & PROTEIN VAULT 
# ==========================================
with tab8:
    st.markdown("<h3 class='The-subheader'>Gene Extraction (Open Reading Frames)</h3>", unsafe_allow_html=True)
    
    if PARQUET_FILE and PARQUET_FILE.exists():
        if 'full_db_dict' not in locals():
            # Load all assembled sequences
            full_db_dict = generate_chromosome_dict(PARQUET_FILE)
        
        if not full_db_dict:
             st.warning("No valid sequences found in database for ORF scanning.")
        else:
            target_seq_id = st.selectbox("Select Sequence to Scan for Genes:", list(full_db_dict.keys()), key="gene_vault_select")
            raw_dna = full_db_dict[target_seq_id]
            
            min_gene_len = st.number_input("Minimum Gene Length (bp):", min_value=100, value=300, step=50)
            
            if st.button("🧬 1. Scan for Reading Frames", type="primary"):
                with st.spinner("Transcribing and Translating DNA..."):
                    seq_obj = Seq(raw_dna)
                    orfs = []
                    
                    translated = seq_obj.translate(table=11) 
                    proteins = str(translated).split("*") 
                    
                    current_pos = 0
                    for p in proteins:
                        if "M" in p:
                            start_idx = p.find("M")
                            actual_protein = p[start_idx:]
                            dna_len = len(actual_protein) * 3
                            
                            if dna_len >= min_gene_len:
                                orfs.append({
                                    "Start (bp)": (current_pos + start_idx) * 3,
                                    "Length (bp)": dna_len,
                                    "Protein Length (AA)": len(actual_protein),
                                    "Protein Sequence": actual_protein,
                                    "NCBI Identity": "Pending..." 
                                })
                        current_pos += len(p) + 1 
                    
                    if orfs:
                        st.session_state['current_orfs'] = pd.DataFrame(orfs).sort_values("Length (bp)", ascending=False)
                        st.success(f"✅ Found {len(orfs)} potential genes!")
                    else:
                        st.warning("No genes found matching the criteria.")

            # If ORFs have been found and saved in session state, show the BLAST UI
            if 'current_orfs' in st.session_state:
                st.dataframe(st.session_state['current_orfs'][['Start (bp)', 'Length (bp)', 'Protein Length (AA)', 'NCBI Identity']], use_container_width=True)
                
                st.markdown("### 🌐 Step 2: Identify ORFs via NCBI")
                st.warning("Note: Searching NCBI takes ~30-60 seconds per sequence. Do not close the window.")
                
                col1, col2 = st.columns([1, 4])
                start_blast = col1.button("🚀 Start NCBI Search")
                stop_blast = col2.button("🛑 Stop Search")
                
                if stop_blast:
                    st.error("Search Cancelled by User.")
                    st.stop() 
                    
                if start_blast:
                    df = st.session_state['current_orfs']
                    progress_bar = st.progress(0)
                    status_text = st.empty()
                    
                    for index, row in df.iterrows():
                        status_text.text(f"Searching ORF at {row['Start (bp)']} bp...")
                        
                        identity = identify_orf_via_blast(row['Protein Sequence'])
                        df.at[index, 'NCBI Identity'] = identity
                        
                        current_prog = (index + 1) / len(df)
                        progress_bar.progress(current_prog)
                        
                        st.session_state['current_orfs'] = df
                        time.sleep(1) 
                        
                    status_text.success("All ORFs processed!")
                    st.rerun() # Replaced deprecated experimental_rerun
    else:
        st.warning("⚠️ Database missing. Please run the Harvester first.")
