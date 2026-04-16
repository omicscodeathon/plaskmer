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
    # 1. Check local scripts folder (The directory where app.py lives)
    path_same_folder = SCRIPT_DIR / filename
    if path_same_folder.exists():
        return path_same_folder
    
    # 2. Check local data folder (One level up, in /data)
    path_local_data = SCRIPT_DIR.parent / "data" / filename
    if path_local_data.exists():
        return path_local_data
        
    # 3. CLOUD FALLBACK: Download directly from Hugging Face
    import streamlit as st
    try:
        hf_token = os.environ.get("HF_TOKEN")
        repo_id = getattr(config, 'HF_REPO_ID', "Jeffiq/Plaskmer") 
        
        # Try looking in the data/ folder first (your new structure)
        try:
            cloud_path = hf_hub_download(
                repo_id=repo_id, 
                filename=f"data/{filename}", 
                repo_type="dataset",
                token=hf_token
            )
        except Exception:
            # Fallback: Look in the main root folder (your old structure)
            cloud_path = hf_hub_download(
                repo_id=repo_id, 
                filename=filename, 
                repo_type="dataset",
                token=hf_token
            )
            
        return Path(cloud_path)
    except Exception as e:
        # We print to console to avoid cluttering the UI with 404s on first runs
        print(f"🔍 Cloud check for '{filename}': {e}")
        return None

# The single, high-performance database file
PARQUET_FILE = find_database_file(config.MASTER_PARQUET)

# 
# --- HUGGING FACE HELPER (Optimized) ---
def push_to_huggingface(local_filepath, repo_filepath):
    hf_token = os.environ.get("HF_TOKEN")
    if not hf_token:
        st.error("⚠️ HF_TOKEN not found. Cannot push to cloud.")
        return False
        
    try:
        api = HfApi(token=hf_token)
        repo_id = getattr(config, 'HF_REPO_ID', "Jeffiq/Plaskmer")
        
        # Ensure we don't end up with data/data/filename
        final_path = repo_filepath
        if not final_path.startswith("data/"):
            final_path = f"data/{final_path}"
            
        api.upload_file(
            path_or_fileobj=str(local_filepath),
            path_in_repo=final_path,
            repo_id=repo_id, 
            repo_type="dataset",
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

# --- UPGRADED K-MER HELPER (PARQUET ENGINE) ---
def update_kmer_database(parquet_path, k=6):
    """
    Finds existing IDs in Parquet, generates K-mers for new sequences, 
    and appends them to the Parquet database.
    """
    pq_filename = f"kmer_{k}mer_vectors.parquet"
    target_path = find_database_file(pq_filename) or (SCRIPT_DIR / pq_filename)
    
    existing_ids = set()
    if target_path.exists():
        try:
            # Load only ID column to check for duplicates
            existing_df = pd.read_parquet(target_path, columns=['sequence_id'])
            existing_ids = set(existing_df['sequence_id'].astype(str))
        except Exception as e:
            print(f"Error reading existing Parquet IDs: {e}")

    if not parquet_path or not parquet_path.exists():
        return target_path, 0
        
    new_records_count = 0
    master_pq = pq.ParquetFile(parquet_path)
    all_new_kmers = []
    
    for batch in master_pq.iter_batches(batch_size=100):
        df_batch = batch.to_pandas()
        
        # Filter sequences not already in K-mer DB
        new_seqs = df_batch[
            (~df_batch['Accession'].astype(str).isin(existing_ids)) & 
            (df_batch['Sequence'].astype(str).str.len() > 0)
        ]
        
        if not new_seqs.empty:
            temp_fasta = SCRIPT_DIR / f"temp_batch_{k}.fasta"
            with open(temp_fasta, "w") as f:
                for _, row in new_seqs.iterrows():
                    f.write(f">{row['Accession']}\n{row['Sequence']}\n")
            
            # Generate K-mers (using your existing tool)
            new_kmer_df = fasta_to_kmer_df(str(temp_fasta), k=k)
            
            if not new_kmer_df.empty:
                # 🚨 THE FIX: Handle the ID whether it's in the index or misnamed
                if 'sequence_id' not in new_kmer_df.columns:
                    # If it's in the index, pull it out into a column
                    new_kmer_df = new_kmer_df.reset_index()
                    # Rename the first column (which was the index) to 'sequence_id'
                    new_kmer_df.rename(columns={new_kmer_df.columns[0]: 'sequence_id'}, inplace=True)
                
                # Now it is guaranteed to exist, so we convert it to string
                new_kmer_df['sequence_id'] = new_kmer_df['sequence_id'].astype(str)
                
                all_new_kmers.append(new_kmer_df)
                new_records_count += len(new_kmer_df)
            
            if temp_fasta.exists():
                os.remove(temp_fasta)
            
    if all_new_kmers:
        combined_new = pd.concat(all_new_kmers, ignore_index=True)
        if target_path.exists():
            existing_full = pd.read_parquet(target_path)
            # Ensure types match before concat
            existing_full['sequence_id'] = existing_full['sequence_id'].astype(str)
            updated_df = pd.concat([existing_full, combined_new], ignore_index=True)
            updated_df.to_parquet(target_path, index=False, engine='pyarrow')
        else:
            combined_new.to_parquet(target_path, index=False, engine='pyarrow')
            
    return target_path, new_records_count

# --- INCREMENTAL ORF EXTRACTOR (BATCHED) ---
def update_orf_database(parquet_path):
    import re
    orf_filename = "orf_database.parquet"
    orf_path = find_database_file(orf_filename) or (SCRIPT_DIR / orf_filename)
    
    existing_ids = set()
    if orf_path.exists():
        # Load only the Accession column to save RAM
        existing_df = pd.read_parquet(orf_path, columns=['Accession'])
        existing_ids = set(existing_df['Accession'].astype(str))
        
    if not parquet_path or not parquet_path.exists():
        return orf_path, 0
        
    new_records_count = 0
    parquet_file = pq.ParquetFile(parquet_path)
    
    # Regex for basic ORF (Start codon, multiple of 3 bases, Stop codon), minimum 300bp
    orf_pattern = re.compile(r'(?=(ATG(?:...){100,}?(?:TAA|TAG|TGA)))')
    all_new_orfs = []
    
    # Process in safe batches of 100
    for batch in parquet_file.iter_batches(batch_size=100):
        df_batch = batch.to_pandas()
        
        # Filter sequences that haven't been processed for ORFs yet
        new_seqs = df_batch[
            (~df_batch['Accession'].astype(str).isin(existing_ids)) & 
            (df_batch['Sequence'].astype(str).str.len() > 0)
        ]
        
        for _, row in new_seqs.iterrows():
            seq = str(row['Sequence'])
            acc = str(row['Accession'])
            seq_len = len(seq)
            
            # Forward Strand ORFs
            for match in orf_pattern.finditer(seq):
                start = match.start()
                orf_seq = match.group(1)
                end = start + len(orf_seq)
                all_new_orfs.append({
                    "Accession": acc,
                    "ORF_ID": f"{acc}_F_{start}",
                    "Start": start, "End": end, "Strand": "Forward",
                    "Length": len(orf_seq), "ORF_Sequence": orf_seq
                })
                
            # Reverse Strand ORFs
            rev_seq = reverse_complement(seq)
            for match in orf_pattern.finditer(rev_seq):
                rev_start = match.start()
                orf_seq = match.group(1)
                actual_end = seq_len - rev_start
                actual_start = actual_end - len(orf_seq)
                all_new_orfs.append({
                    "Accession": acc,
                    "ORF_ID": f"{acc}_R_{actual_start}",
                    "Start": actual_start, "End": actual_end, "Strand": "Reverse",
                    "Length": len(orf_seq), "ORF_Sequence": orf_seq
                })
                
    if all_new_orfs:
        new_orf_df = pd.DataFrame(all_new_orfs)
        if orf_path.exists():
            # Append to the existing Parquet database
            existing_df = pd.read_parquet(orf_path)
            combined = pd.concat([existing_df, new_orf_df], ignore_index=True)
            combined.to_parquet(orf_path, engine='pyarrow')
        else:
            # Create a brand new Parquet database
            new_orf_df.to_parquet(orf_path, engine='pyarrow')
            
        new_records_count = len(new_orf_df)
        
    return orf_path, new_records_count

# --- MATH HELPER ---
def find_closest_matches(query_seq, k=6, top_n=50):
    pq_filename = f"kmer_{k}mer_vectors.parquet"
    # Use your robust find_database_file function
    pq_path = find_database_file(pq_filename) or (SCRIPT_DIR / pq_filename)
    
    # 1. CHECK AND CREATE IF MISSING (Auto-Repair Logic)
    if not pq_path.exists():
        if PARQUET_FILE and PARQUET_FILE.exists():
            with st.spinner(f"🚨 {k}-mer Parquet index not found. Generating from Master Database..."):
                # Use the updated update_kmer_database we built earlier
                pq_path, num_new = update_kmer_database(PARQUET_FILE, k=k)
                
                if num_new > 0:
                    st.success(f"✅ Created {pq_filename} locally with {num_new} sequences.")
                    
                    with st.spinner("☁️ Syncing Parquet index to Hugging Face..."):
                        # Ensure we push to the data/ folder
                        pushed = push_to_huggingface(pq_path, f"data/{pq_path.name}")
                        if pushed:
                            st.success(f"☁️ Cloud Sync Successful: {pq_filename} is now live!")
                        else:
                            st.warning("⚠️ Local file created, but Cloud Sync failed.")
                else:
                    return None, "K-mer calculation returned 0 sequences."
        else:
            return None, "Master Database missing. Please run Harvester first."

    # 2. PROCEED WITH HIGH-SPEED SEARCH
    try:
        # Load Parquet (Much faster than CSV)
        db_df = pd.read_parquet(pq_path)
        
        # Set sequence_id as index if it isn't already
        if 'sequence_id' in db_df.columns:
            db_df = db_df.set_index('sequence_id')
        
        # --- Your original K-mer calculation logic ---
        all_kmers = get_all_kmers(k)
        counts = count_kmers(query_seq, k)
        total = sum(counts.values())
        
        if total == 0: 
            return None, "Query sequence contains no valid k-mers."
        
        # Generate the query vector
        query_vector = np.array([counts.get(k_str, 0) / total for k_str in all_kmers]).reshape(1, -1)
        
        # The database vectors (everything except the index)
        db_vectors = db_df.values

        # Calculate Similarity
        similarities = cosine_similarity(query_vector, db_vectors)[0]
        
        # Build Results
        results_df = pd.DataFrame({
            'Sequence ID': db_df.index.astype(str),
            'Similarity Score': similarities
        }).sort_values(by='Similarity Score', ascending=False).head(top_n)
        
        return results_df, None
        
    except Exception as e:
        return None, f"System Error during similarity search: {str(e)}"
    

# ==========================================
# PAGE CONFIGURATION & MAGNIFICENT UI
# ==========================================
st.set_page_config(page_title="The PlasKmer Studio", page_icon="🧬", layout="wide")
# --- USER SETTINGS & PRIVACY (SIDEBAR) ---




# Update your Entrez email dynamically anywhere in the app like this:
from Bio import Entrez



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
        
        search_triggered = st.button("🚀 Find Closest Match", type="primary",width="stretch")

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
    
    # Initialize session state if not present
    if "user_email" not in st.session_state:
        st.session_state.user_email = ""

    # 🚨 THE GATEKEEPER: Ask for email directly in Tab 2!
    if not st.session_state.user_email:
        st.warning("🛑 **Action Required:** You must provide your NCBI email below to access the Harvester tools.")
        st.info("NCBI requires an email to identify traffic. This is kept strictly in your browser and is never stored on our servers.")
        
        # Put the input box right here in the middle of the screen
        user_email_input = st.text_input("Enter your NCBI registered email:", placeholder="email@example.com")
        
        if st.button("💾 Save & Unlock Harvester", type="primary"):
            if "@" in user_email_input and "." in user_email_input:
                st.session_state.user_email = user_email_input
                st.rerun() # This instantly refreshes the page to show the tools!
            else:
                st.error("❌ Please enter a valid email address.")
    else:
        # If email exists, show the tools!
        colA, colB = st.columns([1.5, 1])
        with colA:
            
            # --- NEW MAINTENANCE SECTION ---
            with st.container(border=True):
                st.markdown("#### 🛠️ Maintenance & Catch-Up")
                st.caption("Analyze existing records in your database that were downloaded before the ORF engine was added.")
                
                if st.button("🔍 Scan Existing Database for Missing ORFs", help="Analyze all saved sequences for ORFs.", width="stretch"):
                    
                    # INJECT EMAIL JUST IN CASE
                    from Bio import Entrez
                    Entrez.email = st.session_state.user_email
                    
                    if PARQUET_FILE and PARQUET_FILE.exists():
                        with st.spinner("Analyzing existing records... this may take a few minutes."):
                            orf_path, num_added = update_orf_database(PARQUET_FILE)
                            
                            if num_added > 0:
                                st.success(f"✅ Success! Found and analyzed {num_added} new ORF entries.")
                                try:
                                    push_to_huggingface(orf_path, orf_path.name)
                                    st.info("☁️ Global ORF Database updated on Hugging Face.")
                                except:
                                    st.warning("ORF Database updated locally, but failed to push to Hugging Face. Check your token.")
                            else:
                                st.info("✨ Everything is up to date! All existing records already have ORF data.")
                    else:
                        st.error("No master database found to scan.")

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
                
                upload_hf = st.toggle("☁️ Push to Hugging Face after harvest?", value=True)
                run_orf_calc = st.checkbox("🧬 Auto-Extract ORFs after harvest?", value=True)
                
                if st.button("🚜 Start Harvester Pipeline", type="primary", width="stretch"):
                    if not target_org:
                        st.error("Please enter a Target Organism.")
                    else:
                        # 🚨 INJECT USER EMAIL INTO GLOBAL ENTREZ
                        from Bio import Entrez
                        Entrez.email = st.session_state.user_email
                        
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
                                                                                    
                        # 1. Trigger Incremental K-mers
                        if run_kmer_calc and PARQUET_FILE and PARQUET_FILE.exists():
                            with st.spinner(f"🧬 Incrementally calculating {k_size}-mers for new sequences..."):
                                kmer_csv_path, num_new_kmers = update_kmer_database(PARQUET_FILE, k=int(k_size))
                                if num_new_kmers > 0:
                                    st.success(f"✅ Added {num_new_kmers} new sequences to K-mer index!")
                        
                        # 2. Trigger Incremental ORF Extraction
                        if run_orf_calc and PARQUET_FILE and PARQUET_FILE.exists():
                            with st.spinner(f"🧬 Extracting global ORFs for new sequences..."):
                                orf_parquet_path, num_new_orfs = update_orf_database(PARQUET_FILE)
                                if num_new_orfs > 0:
                                    st.success(f"✅ Extracted and saved {num_new_orfs} new ORFs to the database!")
                                else:
                                    st.info("⚡ No new ORFs detected.")
                        
                        # 3. Upload updated files to Hugging Face
                        if upload_hf:
                            with st.spinner("☁️ Pushing updated files to Hugging Face..."):
                                if PARQUET_FILE and PARQUET_FILE.exists():
                                    push_to_huggingface(PARQUET_FILE, PARQUET_FILE.name)
                                if run_kmer_calc and 'kmer_csv_path' in locals() and kmer_csv_path.exists():
                                    push_to_huggingface(kmer_csv_path, kmer_csv_path.name)
                                if run_orf_calc and 'orf_parquet_path' in locals() and orf_parquet_path.exists():
                                    push_to_huggingface(orf_parquet_path, orf_parquet_path.name)
                                st.success("✅ Cloud Databases Synced Successfully!")

            # --- MANUAL CLOUD SYNC BLOCK ---
            with st.container(border=True):
                st.markdown("#### ☁️ Manual Cloud Sync")
                st.caption("Push your existing local database directly to Hugging Face without running a new harvest.")
                
                if st.button("⬆️ Push Existing Database to Hugging Face", width="stretch"):
                    import os
                    from huggingface_hub import HfApi
                    import config
                    
                    parquet_file = config.MASTER_PARQUET
                    # 💡 Defaulting to your specified Hugging Face repo
                    repo_id = getattr(config, 'HF_REPO_ID', "Jeffiq/Plaskmer")
                    
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
            
        if st.button("🔬 Run Virtual PCR", type="primary", width="stretch"):
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
                    st.dataframe(pd.DataFrame(all_results), width="stretch")
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
                st.plotly_chart(fig_map, width="stretch")
            with col2:
                fig_types = analytics_engine.plot_type_distribution(df, selected_org)
                st.plotly_chart(fig_types, width="stretch")
            
            # 2. Quality Metrics (Lengths & GC)
            fig_quality = analytics_engine.plot_sequence_quality_metrics(df, selected_org)
            st.plotly_chart(fig_quality, width="stretch")
        
        st.markdown("---")
        
        # --- CROSS-COUNTRY COMPARISON ---
        st.markdown("### 🌍 Cross-Country Comparison Matrix")
        available_countries = sorted(df['Country'].unique().tolist())
        
        # Default to the first two countries if they exist
        default_countries = available_countries[:2] if len(available_countries) >= 2 else available_countries
        selected_countries = st.multiselect("Select Countries to Compare:", available_countries, default=default_countries)
        
        if selected_countries:
            fig_comparison = analytics_engine.plot_cross_country_comparison(df, selected_countries)
            st.plotly_chart(fig_comparison, width="stretch")

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
                st.plotly_chart(fig_type_map, width="stretch")
            
            # 2. Show the Host Organism Bar and Signature Scatter Plot side-by-side
            fig_bar, fig_scatter = analytics_engine.plot_type_scatter_and_bar(df, analysis_type)
            
            if fig_bar and fig_scatter:
                type_col1, type_col2 = st.columns(2)
                with type_col1:
                    st.plotly_chart(fig_bar, width="stretch")
                with type_col2:
                    st.plotly_chart(fig_scatter, width="stretch")

        # --- DATABASE EXPORT / DOWNLOAD ---
        st.markdown("### 💾 Export Database")
        with st.expander("Download Data to Local Disk"):
            dl_col1, dl_col2 = st.columns(2)
            
            with dl_col1:
                # CSV Export
                csv_data = df.to_csv(index=False).encode('utf-8')
                st.download_button(
                    label="📄 Download as CSV",
                    data=csv_data,
                    file_name="plaskmer_database.csv",
                    mime="text/csv",
                    width="stretch"
                )
                
            with dl_col2:
                # FASTA Export (Compressing to avoid massive memory usage)
                import gzip
                def generate_fasta(dataframe):
                    fasta_str = ""
                    for _, row in dataframe.iterrows():
                        fasta_str += f">{row['Accession']} | {row['Organism']} | {row['Country']}\n{row['Sequence']}\n"
                    return fasta_str.encode('utf-8')
                
                fasta_data = generate_fasta(df)
                st.download_button(
                    label="🧬 Download as FASTA",
                    data=fasta_data,
                    file_name="plaskmer_sequences.fasta",
                    mime="text/plain",
                    width="stretch"
                )
        
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
                
                st.plotly_chart(fig_map, width="stretch")
                
                with st.expander("View Geographic Data Table"):
                    st.dataframe(map_data.sort_values('Sequence Count', ascending=False), width="stretch")
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
# TAB 7: PLASMID INTELLIGENCE HUB (Multi-Table View)
# ==========================================
with tab7:
    st.markdown("<h3 class='The-subheader'>Whole Plasmid Analytics & Cross-Species Intelligence</h3>", unsafe_allow_html=True)
    
    if PARQUET_FILE and PARQUET_FILE.exists():
        df_meta = pd.read_parquet(PARQUET_FILE, columns=["Accession", "Sequence", "Organism", "Country", "Type"])
        
        # Robust filtering
        df_plasmids_only = df_meta[
            (df_meta['Type'].str.contains("plasmid", case=False, na=False)) | 
            (df_meta['Accession'].str.contains("plasmid", case=False, na=False)) |
            (df_meta['Accession'].str.startswith(("p", "P")))
        ].copy()
        
        if df_plasmids_only.empty:
            st.warning("⚠️ No records explicitly labeled 'Plasmid' found. Showing all sequences.")
            df_plasmids_only = df_meta.copy()

        # --- SELECTION ---
        col_sel1, col_sel2 = st.columns(2)
        with col_sel1:
            available_hosts = sorted(df_plasmids_only['Organism'].dropna().unique().tolist())
            selected_host = st.selectbox("🦠 Select Host Organism:", available_hosts, key="tab7_org_v2")
        with col_sel2:
            filtered_list = df_plasmids_only[df_plasmids_only['Organism'] == selected_host]['Accession'].tolist()
            selected_plasmid_id = st.selectbox("🧬 Select Specific Sequence:", filtered_list, key="tab7_acc_v2")
        
        current_data = df_plasmids_only[df_plasmids_only['Accession'] == selected_plasmid_id].iloc[0]
        selected_seq = current_data['Sequence']
        
        st.info(f"📍 **Origin:** {current_data['Country']} | **Length:** {len(selected_seq):,} bp")

       # --- VISUALIZER & ORFs ---
        col_vis, col_orf = st.columns([1.2, 1])
        
        with col_vis:
            st.markdown("#### 🔄 Circular Visualizer")
            if len(selected_seq) > 0:
                try:
                    # 1. Base Feature (The whole plasmid)
                    features = [
                        GraphicFeature(start=0, end=len(selected_seq), strand=+1, color="#cddc39", label="Backbone")
                    ]
                    
                    # 2. Add Restriction Sites to make the map look professional
                    seq_obj = Seq(selected_seq)
                    for site in Restriction.EcoRI.search(seq_obj):
                        features.append(GraphicFeature(start=site-1, end=site, strand=+1, color="#ff4b4b", label="EcoRI"))
                    for site in Restriction.BamHI.search(seq_obj):
                        features.append(GraphicFeature(start=site-1, end=site, strand=+1, color="#2e66ff", label="BamHI"))

                    # 3. Draw the Map
                    record = CircularGraphicRecord(sequence_length=len(selected_seq), features=features)
                    fig, ax = plt.subplots(figsize=(5, 5))
                    record.plot(ax=ax)
                    st.pyplot(fig)
                except Exception as e:
                    st.error(f"Could not generate visual: {e}")
            else:
                st.warning("Sequence length is 0. Cannot draw map.")
            
        with col_orf:
            st.markdown("#### 🧬 Sequence ORFs")
            
            # Look for the global ORF database
            orf_filename = "orf_database.parquet"
            orf_path = find_database_file(orf_filename) or (SCRIPT_DIR / orf_filename)
            
            if orf_path.exists():
                try:
                    df_orfs = pd.read_parquet(orf_path)
                    
                    # Filter for only the selected plasmid
                    plasmid_orfs = df_orfs[df_orfs['Accession'].astype(str) == str(selected_plasmid_id)]
                    
                    if not plasmid_orfs.empty:
                        # Display a clean table of the ORFs
                        st.dataframe(
                            plasmid_orfs[['Start', 'End', 'Strand', 'Length']].sort_values(by="Length", ascending=False), 
                            use_container_width=True, 
                            hide_index=True,
                            height=350 # Locks the height so it aligns nicely next to the circle
                        )
                        st.caption(f"Found **{len(plasmid_orfs)}** ORFs. Showing longest first.")
                    else:
                        st.info("No ORFs extracted for this sequence yet. Go to Tab 2 and run the Maintenance scan!")
                except Exception as e:
                    st.error(f"Error loading ORFs: {e}")
            else:
                st.warning("⚠️ ORF Database not found. Please run the Harvester or Maintenance scan in Tab 2 first.")

        st.markdown("---")

        
        # --- THE TWO TABLES: COMPARATIVE INTELLIGENCE ---
        st.markdown("### 🌍 Comparative Intelligence Analysis")
        
        with st.spinner("Analyzing genomic relationships..."):
            sim_results, error = find_closest_matches(selected_seq, k=6, top_n=50)
            
            if not error and sim_results is not None:
                sim_results['Sequence ID'] = sim_results['Sequence ID'].astype(str)
                df_meta['Accession'] = df_meta['Accession'].astype(str)
                
                full_sim = pd.merge(sim_results, df_meta[['Accession', 'Organism', 'Country']], 
                                    left_on='Sequence ID', right_on='Accession')
                full_sim = full_sim[full_sim['Accession'] != selected_plasmid_id]

                # --- TABLE 1: INTERNAL SIMILARITY ---
                st.markdown(f"#### 🔍 1. Relatives within **{selected_host}**")
                same_host_df = full_sim[full_sim['Organism'] == selected_host].head(5)
                if not same_host_df.empty:
                    st.dataframe(same_host_df[['Accession', 'Country', 'Similarity Score']], 
                                 use_container_width=True, hide_index=True)
                else:
                    st.info("No close relatives found in the same species.")

                st.markdown("<br>", unsafe_allow_html=True)

                # --- TABLE 2: CROSS-SPECIES INTELLIGENCE ---
                st.markdown("#### 🚀 2. Cross-Species Relatives (Potential Host Jumps)")
                diff_host_df = full_sim[full_sim['Organism'] != selected_host].head(10)
                if not diff_host_df.empty:
                    st.dataframe(diff_host_df[['Accession', 'Organism', 'Country', 'Similarity Score']], 
                                 use_container_width=True, hide_index=True)
                    st.caption("🚨 High similarity in different species suggests potential Horizontal Gene Transfer (HGT).")
                else:
                    st.info("No cross-species relatives detected in current database.")
            else:
                st.error(f"Similarity Engine Error: {error}")
    else:
        st.warning("⚠️ Database missing.")


        
# ==========================================
# TAB 8: GLOBAL GENE & PROTEIN VAULT (Batch Mode)
# ==========================================
with tab8:
    st.markdown("<h3 class='The-subheader'>Global Gene Extraction (Open Reading Frames)</h3>", unsafe_allow_html=True)
    
    if PARQUET_FILE and PARQUET_FILE.exists():
        # Load full metadata for filtering
        df_all = pd.read_parquet(PARQUET_FILE, columns=["Accession", "Organism", "Type", "Sequence"])
        
        # --- BATCH FILTERS ---
        st.markdown("#### 🛠️ Batch Processing Filters")
        c1, c2, c3 = st.columns(3)
        
        with c1:
            org_list = ["ALL Organisms"] + sorted(df_all['Organism'].unique().tolist())
            sel_org = st.selectbox("Filter by Organism:", org_list)
        
        with c2:
            type_list = ["ALL Types", "Plasmids Only", "Chromosomes Only"]
            sel_type = st.selectbox("Filter by Genomic Type:", type_list)
            
        with c3:
            min_gene_len = st.number_input("Min Gene Length (bp):", min_value=100, value=300)

        # Apply Filters to the dataframe
        df_filtered = df_all.copy()
        if sel_org != "ALL Organisms":
            df_filtered = df_filtered[df_filtered['Organism'] == sel_org]
        
        if sel_type == "Plasmids Only":
            df_filtered = df_filtered[df_filtered['Type'].str.contains("Plasmid", case=False, na=False)]
        elif sel_type == "Chromosomes Only":
            df_filtered = df_filtered[~df_filtered['Type'].str.contains("Plasmid", case=False, na=False)]

        st.metric("Sequences to Scan", len(df_filtered))

        if st.button("🚀 Run Global Scan on Filtered Database", type="primary"):
            all_discovered_orfs = []
            progress_bar = st.progress(0)
            status_text = st.empty()
            
            for i, (idx, row) in enumerate(df_filtered.iterrows()):
                acc = row['Accession']
                status_text.text(f"Scanning {acc} ({i+1}/{len(df_filtered)})...")
                
                # ORF Extraction Logic
                seq_obj = Seq(row['Sequence'])
                # Scan 6 frames (3 forward, 3 reverse)
                for strand, nuc in [(1, seq_obj), (-1, seq_obj.reverse_complement())]:
                    translated = nuc.translate(table=11)
                    proteins = str(translated).split("*")
                    
                    curr_pos = 0
                    for p in proteins:
                        if "M" in p:
                            start_idx = p.find("M")
                            actual_protein = p[start_idx:]
                            dna_len = len(actual_protein) * 3
                            
                            if dna_len >= min_gene_len:
                                all_discovered_orfs.append({
                                    "Accession": acc,
                                    "Organism": row['Organism'],
                                    "Start (bp)": (curr_pos + start_idx) * 3,
                                    "Length": dna_len,
                                    "Protein": actual_protein,
                                    "Identity": "Not Identified"
                                })
                        curr_pos += len(p) + 1
                
                progress_bar.progress((i + 1) / len(df_filtered))
            
            if all_discovered_orfs:
                st.session_state['global_orfs'] = pd.DataFrame(all_discovered_orfs)
                st.success(f"✅ Found {len(all_discovered_orfs)} total genes across the selection!")
            else:
                st.warning("No genes found matching those criteria.")

        # --- DISPLAY RESULTS & BLAST ---
        if 'global_orfs' in st.session_state:
            res_df = st.session_state['global_orfs']
            
            st.markdown("### 📊 Discovered Genes")
            # Sub-filter results table
            search_query = st.text_input("🔍 Search genes by Accession or Organism:")
            if search_query:
                res_df = res_df[res_df['Accession'].str.contains(search_query, case=False) | 
                                res_df['Organism'].str.contains(search_query, case=False)]
            
            st.dataframe(res_df[['Accession', 'Organism', 'Length', 'Identity']], use_container_width=True)
            
            if st.button("🌐 Identify Top 5 Longest Genes via NCBI"):
                # BLASTing a whole database is too slow for Streamlit, 
                # so we do the top N most interesting ones.
                top_orfs = res_df.sort_values("Length", ascending=False).head(5)
                for i, row in top_orfs.iterrows():
                    with st.spinner(f"Identifying {row['Accession']} gene..."):
                        id_result = identify_orf_via_blast(row['Protein'])
                        st.session_state['global_orfs'].at[i, 'Identity'] = id_result
                st.rerun()
    else:
        st.warning("⚠️ Database missing.")
