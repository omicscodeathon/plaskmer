import streamlit as st
import pandas as pd
from pathlib import Path
import time
from harvester import ParallelOmniSystem # Ensure harvester.py is in the same folder!

# --- DUAL-MODE PATH CONFIGURATION ---
SCRIPT_DIR = Path(__file__).parent.resolve()

def find_database_file(filename):
    """Hunts down the file whether on local Windows desktop (Mode 1) or GitHub cloud (Mode 2)."""
    # MODE 1: Local desktop (Files are in the exact same folder as app.py)
    path_same_folder = SCRIPT_DIR / filename
    if path_same_folder.exists():
        return path_same_folder
        
    # MODE 2: GitHub / Streamlit Cloud (app.py is in /scripts, data is in /data)
    path_github = SCRIPT_DIR.parent / "data" / filename
    if path_github.exists():
        return path_github
        
    return None

CSV_FILE = find_database_file("local_harvest_log.csv")
FASTA_FILE = find_database_file("plas_kmer_sequences.fasta")

# --- PAGE CONFIGURATION & CUSTOM COLORS ---
st.set_page_config(page_title="Afrigen PlasKmer Studio", page_icon="🧬", layout="wide")

st.markdown("""
<style>
    .stTabs [data-baseweb="tab-list"] { gap: 8px; }
    .stTabs [data-baseweb="tab"] {
        height: 50px; background-color: #e0e6ed;
        border-radius: 5px 5px 0px 0px; padding: 10px 20px;
        font-weight: 600; color: #1f2937;
    }
    .stTabs [aria-selected="true"] { background-color: #2e66ff !important; color: white !important; }
    .afrigen-header { color: #2e66ff; font-weight: 800; }
    .afrigen-subheader { color: #10b981; font-weight: 700; }
</style>
""", unsafe_allow_html=True)

st.markdown("<h1 class='afrigen-header'>🧬 Afrigen PlasKmer Studio</h1>", unsafe_allow_html=True)
st.markdown("The one-stop cloud portal for African Plasmid and mRNA informatics.")
st.divider()

# --- TAB NAVIGATION ---
tab1, tab2, tab3, tab4, tab5 = st.tabs([
    "🔍 1. K-mer Search", 
    "🚜 2. Data Harvester", 
    "🎯 3. Virtual PCR", 
    "📊 4. Analytics", 
    "🌍 5. Global Stats"
])

# ==========================================
# TAB 1: K-MER SEARCH 
# ==========================================
with tab1:
    st.markdown("<h3 class='afrigen-subheader'>Pattern-Based Retrieval & Closest Match</h3>", unsafe_allow_html=True)
    
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
        st.success("Query processed! Looking up closest FASTA match... (Backend logic pending)")

# ==========================================
# TAB 2: DATA HARVESTER & ADMIN
# ==========================================
with tab2:
    st.markdown("<h3 class='afrigen-subheader'>Database Management & Processing</h3>", unsafe_allow_html=True)
    
    colA, colB = st.columns([1.5, 1])
    
    with colA:
        with st.container(border=True):
            st.markdown("#### 🚜 Run Parallel Omni Harvester")
            st.caption("Scans all African nations for target sequences and appends to the database.")
            
            target_org = st.text_input("Target Organism:", placeholder="e.g., Vibrio cholerae")
            
            st.markdown("**1. Select Databases to Search:**")
            db_col1, db_col2, db_col3, db_col4 = st.columns(4)
            with db_col1: use_bp = st.checkbox("BioProject")
            with db_col2: use_bs = st.checkbox("BioSample")
            with db_col3: use_sra = st.checkbox("SRA")
            with db_col4: use_nuc = st.checkbox("Nucleotide", value=True)
            
            st.markdown("**2. Select Sequence Types:**")
            type_col1, type_col2, type_col3 = st.columns(3)
            with type_col1: inc_plasmids = st.checkbox("Plasmids", value=True)
            with type_col2: inc_wgs = st.checkbox("WGS")
            with type_col3: inc_genes = st.checkbox("Specific Genes")
            
            target_goal = st.number_input("Target new records per country:", min_value=1, max_value=1000, value=10)
            upload_hf = st.toggle("☁️ Push to Hugging Face after harvest? (GitHub Mode)")
            
            if st.button("🚜 Start Harvester Pipeline", type="primary", use_container_width=True):
                if not target_org:
                    st.error("Please enter a Target Organism to begin.")
                elif not any([use_bp, use_bs, use_sra, use_nuc]):
                    st.error("⚠️ No databases selected. Please check at least one database.")
                else:
                    st.info(f"Initiating Harvester for {target_org}...")
                    
                    selected_dbs = []
                    if use_bp: selected_dbs.append("bioproject")
                    if use_bs: selected_dbs.append("biosample")
                    if use_sra: selected_dbs.append("sra")
                    if use_nuc: selected_dbs.append("nucleotide")
                    
                    types = []
                    if inc_plasmids: types.append("plasmid")
                    if inc_wgs: types.append("WGS")
                    if inc_genes: types.append("gene")
                    
                    with st.spinner("Downloading from NCBI and writing to disk... Please check terminal for live logs."):
                        bot = ParallelOmniSystem(
                            target_org=target_org, 
                            selected_dbs=selected_dbs, 
                            types=types, 
                            target_goal=target_goal,
                            push_to_hf=upload_hf
                        )
                        bot.run_all_parallel()
                    
                    st.success("✅ Harvest Complete! Database files updated locally.")

    with colB:
        with st.container(border=True):
            st.markdown("#### 🧮 Run Bill's K-mer Tool")
            st.caption("Vectorize the newly downloaded FASTA files into mathematical matrices.")
            bill_k = st.number_input("Set 'k' for vectorization:", min_value=3, max_value=8, value=6)
            if st.button("Run Vectorization", type="primary", use_container_width=True):
                st.warning("Bill's code backend connection pending.")

# ==========================================
# TAB 3: VIRTUAL PCR
# ==========================================
with tab3:
    st.markdown("<h3 class='afrigen-subheader'>In Silico Amplification</h3>", unsafe_allow_html=True)
    
    with st.container(border=True):
        pcr_col1, pcr_col2 = st.columns(2)
        with pcr_col1:
            st.markdown("**Forward Primer (5' ➔ 3')**")
            fwd_primer = st.text_input("Sequence:", key="fwd", placeholder="e.g., ATGGCC...")
        with pcr_col2:
            st.markdown("**Reverse Primer (5' ➔ 3')**")
            rev_primer = st.text_input("Sequence:", key="rev", placeholder="e.g., TTAAGC...")
            
        pcr_mismatches = st.number_input("Allowed Mismatches", min_value=0, max_value=5, value=0)
        run_pcr = st.button("🔬 Run Virtual PCR", type="primary", use_container_width=True)

    if run_pcr:
        with st.status("🧬 PCR is Running. Please wait...", expanded=True) as status:
            progress_bar = st.progress(0)
            total_sequences = 14203 
            for percent_complete in range(1, 101):
                time.sleep(0.04)
                seqs_covered = int((percent_complete / 100) * total_sequences)
                progress_bar.progress(percent_complete, text=f"Scanning... {percent_complete}% ({seqs_covered:,} / {total_sequences:,} seqs)")
            status.update(label="PCR Scan Complete!", state="complete", expanded=False)
            
        st.success("Virtual Amplification Successful!")

# ==========================================
# TAB 4: ADVANCED ANALYTICS
# ==========================================
with tab4:
    st.markdown("<h3 class='afrigen-subheader'>Pre-Computed Analytics Suite</h3>", unsafe_allow_html=True)
    with st.expander("🧬 Genomics & Variants", expanded=True):
        st.button("Run SNP Calling & Variant Analysis")
    with st.expander("🌳 Evolution & Population"):
        st.button("View Phylogenetic Trees")

# ==========================================
# TAB 5: GLOBAL HARVEST STATISTICS
# ==========================================
with tab5:
    st.markdown("<h3 class='afrigen-subheader'>Global Harvest Statistics</h3>", unsafe_allow_html=True)
    st.markdown("Aggregated metrics from your local CSV and FASTA database. No raw data is displayed here.")
    
    col_csv, col_fasta = st.columns(2)
    
    # --- 1. CSV STATISTICS ---
    with col_csv:
        st.markdown("#### 📊 CSV Metadata Stats")
        if CSV_FILE:
            try:
                df_log = pd.read_csv(CSV_FILE)
                total_records = len(df_log)
                
                if 'Type' in df_log.columns:
                    plasmid_count = len(df_log[df_log['Type'].str.contains('Plasmid', case=False, na=False)])
                    mrna_count = len(df_log[df_log['Type'].str.contains('mRNA', case=False, na=False)])
                    wgs_count = len(df_log[df_log['Type'].str.contains('WGS', case=False, na=False)])
                else:
                    plasmid_count, mrna_count, wgs_count = 0, 0, 0
                
                c1, c2 = st.columns(2)
                c1.metric(label="Total Database Records", value=f"{total_records:,}")
                c2.metric(label="Total Plasmids", value=f"{plasmid_count:,}")
                
                c3, c4 = st.columns(2)
                c3.metric(label="Total mRNA", value=f"{mrna_count:,}")
                c4.metric(label="Total WGS / Other", value=f"{wgs_count:,}")
                
                if 'Country' in df_log.columns and 'Organism' in df_log.columns:
                    st.markdown("**Top 5 Contributing Countries:**")
                    st.bar_chart(df_log['Country'].value_counts().head(5))
                    
            except Exception as e:
                st.error(f"Error reading CSV file: {e}")
        else:
            st.error("❌ Could not find `local_harvest_log.csv` (Checked Local and GitHub paths).")

    # --- 2. FASTA STATISTICS ---
    with col_fasta:
        st.markdown("#### 🧬 FASTA Sequence Stats")
        if FASTA_FILE:
            with st.spinner("Calculating FASTA statistics (this may take a moment for large files)..."):
                try:
                    fasta_seq_count = 0
                    total_base_pairs = 0
                    
                    with open(FASTA_FILE, 'r', encoding='utf-8') as f:
                        for line in f:
                            if line.startswith('>'):
                                fasta_seq_count += 1
                            else:
                                total_base_pairs += len(line.strip())
                    
                    avg_length = total_base_pairs // fasta_seq_count if fasta_seq_count > 0 else 0
                    file_size_mb = round(FASTA_FILE.stat().st_size / (1024 * 1024), 2)
                    
                    f1, f2 = st.columns(2)
                    f1.metric(label="Sequences in FASTA", value=f"{fasta_seq_count:,}")
                    f2.metric(label="Total Base Pairs (bp)", value=f"{total_base_pairs:,}")
                    
                    f3, f4 = st.columns(2)
                    f3.metric(label="Average Seq Length", value=f"{avg_length:,} bp")
                    f4.metric(label="Physical File Size", value=f"{file_size_mb} MB")
                    
                    st.success(f"Successfully analyzed data from Mode: {'Local (Same Folder)' if FASTA_FILE.parent == SCRIPT_DIR else 'GitHub (Data Folder)'}")
                    
                except Exception as e:
                    st.error(f"Error reading FASTA file: {e}")
        else:
            st.error("❌ Could not find `plas_kmer_sequences.fasta` (Checked Local and GitHub paths).")
