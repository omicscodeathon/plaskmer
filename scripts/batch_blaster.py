import pandas as pd
import time
import os
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

ORF_FILE = "orf_database.parquet"

def clean_identity(title):
    """Extracts a readable gene name from the messy NCBI title."""
    parts = title.split("|")
    return parts[-1].strip() if len(parts) > 0 else title

def run_batch_blast():
    print("🚀 Initializing Offline Massive Database Updater...")
    
    if not os.path.exists(ORF_FILE):
        print(f"❌ Error: {ORF_FILE} not found.")
        return
        
    df = pd.read_parquet(ORF_FILE)
    
    # Find rows that need updating
    if 'Identity' not in df.columns:
        df['Identity'] = "Not Identified"
        
    unidentified_mask = df['Identity'] == "Not Identified"
    unidentified_count = unidentified_mask.sum()
    
    print(f"🔍 Found {unidentified_count} ORFs needing identification.")
    if unidentified_count == 0:
        print("✅ Database is fully identified. Exiting.")
        return

    # Process in batches to save to Parquet frequently
    processed = 0
    for idx, row in df[unidentified_mask].iterrows():
        protein_seq = row.get('Protein', '')
        if len(protein_seq) < 30: # Too short for meaningful blast
            df.at[idx, 'Identity'] = "Too Short"
            continue
            
        print(f"📡 BLASTing ORF from {row['Accession']} (Length: {len(protein_seq)} aa)...")
        
        success = False
        retries = 3
        while not success and retries > 0:
            try:
                # Use blastp (Protein to Protein)
                result_handle = NCBIWWW.qblast("blastp", "nr", protein_seq, hitlist_size=1)
                blast_record = NCBIXML.read(result_handle)
                
                if blast_record.alignments:
                    top_hit = blast_record.alignments[0]
                    gene_name = clean_identity(top_hit.title)
                    df.at[idx, 'Identity'] = gene_name
                    print(f"   ✅ Identified as: {gene_name}")
                else:
                    df.at[idx, 'Identity'] = "No Match Found"
                    print("   ❌ No Match Found")
                
                success = True
                
            except Exception as e:
                retries -= 1
                print(f"   ⚠️ NCBI Error: {e}. Retries left: {retries}")
                time.sleep(10) # Wait longer on error
                
        # Extremely conservative rate limit (NCBI will ban you if you hammer qblast)
        time.sleep(15) 
        
        processed += 1
        # Save every 50 records so progress isn't lost if it crashes
        if processed % 50 == 0:
            print(f"💾 Saving progress to {ORF_FILE}...")
            df.to_parquet(ORF_FILE, index=False)

    # Final save
    print(f"💾 Final Save to {ORF_FILE}...")
    df.to_parquet(ORF_FILE, index=False)
    print("🎉 Massive Update Complete!")

if __name__ == "__main__":
    run_batch_blast()
