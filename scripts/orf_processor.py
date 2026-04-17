import pandas as pd
import os
import gc
from huggingface_hub import HfApi
from Bio.Seq import Seq

# --- HELPER: PUSH TO HF ---
def push_to_huggingface(file_path, filename):
    """Local helper to push the checkpointed file to Hugging Face."""
    import config
    repo_id = getattr(config, 'HF_REPO_ID', None)
    if not repo_id or repo_id == "YourUsername/PlasKmer-Database":
        print("⚠️ No valid HF_REPO_ID to push to.")
        return
    try:
        api = HfApi()
        api.upload_file(
            path_or_fileobj=file_path,
            path_in_repo=f"data/{filename}",
            repo_id=repo_id,
            repo_type="dataset"
        )
    except Exception as e:
        print(f"❌ Failed to push {filename} to HF: {e}")


# --- HELPER: BIOPYTHON ORF FINDER ---
def extract_orfs_from_sequence(sequence, accession, min_protein_len=100):
    """Scans all 6 reading frames to find open reading frames (genes)."""
    orfs = []
    seq_obj = Seq(sequence)
    
    # Check both forward (+1) and reverse (-1) strands
    for strand, nuc in [(+1, seq_obj), (-1, seq_obj.reverse_complement())]:
        for frame in range(3):
            # Translate DNA to Protein using Bacterial table (11)
            trans = str(nuc[frame:].translate(table=11))
            aa_start = 0
            
            while aa_start < len(trans):
                aa_start = trans.find("M", aa_start) # Look for Start Codon (Methionine)
                if aa_start == -1: break
                
                aa_end = trans.find("*", aa_start)   # Look for Stop Codon
                if aa_end == -1: break
                
                # If the gene is long enough, save it
                if (aa_end - aa_start) >= min_protein_len:
                    protein_seq = trans[aa_start:aa_end]
                    
                    # Calculate exact DNA positions for the visualizer
                    if strand == 1:
                        start_nuc = frame + (aa_start * 3)
                        end_nuc = start_nuc + len(protein_seq) * 3 + 3
                    else:
                        end_nuc = len(seq_obj) - (frame + (aa_start * 3))
                        start_nuc = end_nuc - len(protein_seq) * 3 - 3
                    
                    orfs.append({
                        'Accession': str(accession),
                        'Start': start_nuc,
                        'End': end_nuc,
                        'Strand': strand,
                        'Length': len(protein_seq),
                        'Protein': protein_seq
                    })
                aa_start = aa_end + 1
    return orfs


# --- MAIN ENGINE: BATCH AND FLUSH ---
def safe_batch_orf_extractor(master_df, orf_parquet_path, batch_size=25, push_to_hf=True):
    """Extracts ORFs in batches of 25, saves to disk, syncs to HF, and flushes RAM."""
    
    print(f"🧬 Starting memory-safe ORF extraction (Batch size: {batch_size})...")
    batch_records = []
    
    # Keep track of Accessions we already processed to avoid duplicates
    existing_accs = set()
    if os.path.exists(orf_parquet_path):
        try:
            df_existing = pd.read_parquet(orf_parquet_path, columns=['Accession'])
            existing_accs = set(df_existing['Accession'].astype(str).unique())
        except:
            pass

    for index, row in master_df.iterrows():
        acc = str(row['Accession'])
        
        # Skip if we already processed this plasmid in a previous run
        if acc in existing_accs:
            continue
            
        # 1. Extract the ORFs using Biopython
        new_orfs = extract_orfs_from_sequence(row['Sequence'], acc)
        batch_records.extend(new_orfs)
        
        # 2. THE MAGIC BATCH CHECK
        if (index + 1) % batch_size == 0 and batch_records:
            print(f"📦 Reached {batch_size} sequences! Saving ORFs to disk...")
            
            df_batch = pd.DataFrame(batch_records)
            if os.path.exists(orf_parquet_path):
                df_existing = pd.read_parquet(orf_parquet_path)
                df_combined = pd.concat([df_existing, df_batch], ignore_index=True)
                df_combined.to_parquet(orf_parquet_path, index=False)
            else:
                df_batch.to_parquet(orf_parquet_path, index=False)
            
            if push_to_hf:
                print("☁️ Syncing ORF checkpoint to Hugging Face...")
                push_to_huggingface(orf_parquet_path, os.path.basename(orf_parquet_path))
            
            # 🧹 TAKE OUT THE TRASH
            print("🧹 Flushing RAM...")
            batch_records.clear()
            gc.collect()
            
    # 3. CLEAN UP THE REMAINDER
    if batch_records:
        df_batch = pd.DataFrame(batch_records)
        if os.path.exists(orf_parquet_path):
            df_existing = pd.read_parquet(orf_parquet_path)
            df_combined = pd.concat([df_existing, df_batch], ignore_index=True)
            df_combined.to_parquet(orf_parquet_path, index=False)
        else:
            df_batch.to_parquet(orf_parquet_path, index=False)
            
        if push_to_hf:
            push_to_huggingface(orf_parquet_path, os.path.basename(orf_parquet_path))
            
        batch_records.clear()
        gc.collect()
        
    print("✅ All ORFs extracted and safely stored!")
