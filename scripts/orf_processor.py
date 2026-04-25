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
    """Scans all 6 reading frames to find DNA Nucleotide ORF sequences."""
    orfs = []
    from Bio.Seq import Seq
    seq_obj = Seq(sequence)
    
    # Check both forward (+1) and reverse (-1) strands
    for strand, nuc in [(+1, seq_obj), (-1, seq_obj.reverse_complement())]:
        for frame in range(3):
            # Internal translation just to FIND the Start (M) and Stop (*) codons
            trans = str(nuc[frame:].translate(table=11))
            aa_start = 0
            
            while aa_start < len(trans):
                aa_start = trans.find("M", aa_start) 
                if aa_start == -1: break
                
                aa_end = trans.find("*", aa_start)   
                if aa_end == -1: break
                
                # Check if it meets your minimum length requirement
                if (aa_end - aa_start) >= min_protein_len:
                    # DNA length = (AA length * 3) + 3 for the stop codon
                    dna_len = ((aa_end - aa_start) * 3) + 3
                    
                    # Slicing the actual Nucleotide DNA
                    if strand == 1:
                        start_nuc = frame + (aa_start * 3)
                        end_nuc = start_nuc + dna_len
                        orf_dna = str(seq_obj[start_nuc:end_nuc])
                    else:
                        # Slice from reverse complement and handle coordinates
                        rev_start = frame + (aa_start * 3)
                        rev_end = rev_start + dna_len
                        orf_dna = str(seq_obj.reverse_complement()[rev_start:rev_end])
                        
                        # Forward coordinates for the visualizer
                        end_nuc = len(seq_obj) - rev_start
                        start_nuc = end_nuc - dna_len
                    
                    orfs.append({
                        'Accession': str(accession),
                        'Start': start_nuc,
                        'End': end_nuc,
                        'Strand': strand,  # FIXED: Matches old database format (1 or -1)
                        'Length': len(orf_dna),
                        'Sequence': orf_dna # <-- SAVING THE DNA NUCLEOTIDES
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
