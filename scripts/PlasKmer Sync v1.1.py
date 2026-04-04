import gspread
import time
import os
import re
from Bio import Entrez
from datetime import datetime

# --- CONFIGURATION ---
DEST_CONFIG = {
    "JSON": r"C:\Users\metha\OneDrive\Documents\CODEATHON 2026 APril\plaskmer-93160b8d056d.json",
    "URL": "https://docs.google.com/spreadsheets/d/1sejGoSyGfvsT1L3UXiUKt-gmb-dLrCGPbo52-UAbd0w/edit",
    "TAB": "PlasKmer"
}

# MUST set your email for NCBI
Entrez.email = "methajefferson@gmail.com"
SEQUENCE_FILE = "plas_kmer_sequences.fasta"

class PlasKmerMasterHarvester:
    def __init__(self):
        print("⌛ Connecting to Google Sheets...")
        try:
            self.gc = gspread.service_account(filename=DEST_CONFIG["JSON"])
            self.sheet = self.gc.open_by_url(DEST_CONFIG["URL"]).worksheet(DEST_CONFIG["TAB"])
        except Exception as e:
            print(f"❌ Connection Error: {e}"); exit()

    def get_ncbi_info(self, accession):
        """
        Fetches both the Sequence Type (Plasmid/Genomic) 
        and the DNA Sequence from NCBI using DB fallback.
        """
        for db in ["nuccore", "nucleotide"]:
            try:
                time.sleep(0.4) # Respect NCBI rate limits
                # 1. Fetch GenBank for Type detection
                handle_gb = Entrez.efetch(db=db, id=accession, rettype="gb", retmode="text")
                raw_gb = handle_gb.read().lower()
                
                seq_type = "Genomic DNA"
                if "plasmid" in raw_gb: seq_type = "Plasmid"
                elif "chromosome" in raw_gb: seq_type = "Chromosome"
                elif "wgs" in raw_gb or "shotgun" in raw_gb: seq_type = "WGS/Genomic"

                # 2. Fetch FASTA for the DNA string
                time.sleep(0.4)
                handle_fa = Entrez.efetch(db=db, id=accession, rettype="fasta", retmode="text")
                raw_fa = handle_fa.read().strip()
                
                # Extract DNA only (skip NCBI's default header)
                dna_lines = raw_fa.split('\n')
                dna_string = "\n".join(dna_lines[1:]) if dna_lines[0].startswith('>') else raw_fa
                
                return seq_type, dna_string
            except:
                continue
        return None, None

    def run(self):
        # 1. Read the sheet
        rows = self.sheet.get_all_values()
        if len(rows) < 2:
            print("No data in sheet."); return
        
        data = rows[1:] # Skip headers
        
        # 2. DELETE/WIPE the existing FASTA file to start fresh
        if os.path.exists(SEQUENCE_FILE):
            os.remove(SEQUENCE_FILE)
            print(f"🗑️ Old {SEQUENCE_FILE} deleted. Starting fresh.")

        print(f"🧬 Harvesting {len(data)} records from NCBI...")

        # 3. Create the new file in 'w' (write) mode - NO COMMENTS ALLOWED
        with open(SEQUENCE_FILE, "w", encoding="utf-8") as f:
            for row in data:
                # Map based on your Sheet columns
                u_id = row[0]
                acc = row[1]
                org = row[2]
                country = row[7] if len(row) > 7 else "Unknown"

                print(f"📡 Processing {acc}...", end="\r")

                # 4. Get the LIVE data from NCBI
                seq_type, dna = self.get_ncbi_info(acc)

                if dna and seq_type:
                    # 5. Write STRICT FASTA FORMAT
                    # The very first line of the file will be a '>' header.
                    f.write(f">{u_id} | {acc} | {org} | {seq_type} | {country}\n")
                    f.write(dna + "\n\n")
                else:
                    print(f"\n❌ Failed to harvest {acc} from NCBI.")

        print(f"\n\n✨ HARVEST COMPLETE.")
        print(f"📂 File saved: {SEQUENCE_FILE}")
        print("✅ Format: Pure FASTA (No comments, correct headers).")

if __name__ == "__main__":
    PlasKmerMasterHarvester().run()
