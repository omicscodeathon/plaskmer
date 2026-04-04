import gspread
from Bio import Entrez, SeqIO
import io, time, logging, re, xml.etree.ElementTree as ET
from datetime import datetime
import constants 

# --- CONFIG ---
CONFIG = {
    "JSON_PATH": r"C:\Users\metha\OneDrive\Documents\CODEATHON 2026 APril\hypnotic-epoch-452109-e0-d0449e621117.json",
    "SHEET_URL": "https://docs.google.com/spreadsheets/d/1viduBEfMkjZOhPafJaUIdV66vQS2sm8rO0PnRe1ryGQ/edit",
    "EMAIL": "methajefferson@gmail.com"
}

Entrez.email = CONFIG["EMAIL"]
logging.basicConfig(level=logging.INFO, format="%(message)s")

class PersistentAfricanHarvester:
    def __init__(self):
        try:
            self.gc = gspread.service_account(filename=CONFIG["JSON_PATH"])
            self.doc = self.gc.open_by_url(CONFIG["SHEET_URL"])
            self.master = self.doc.worksheet("Data Curation Sheet")
        except Exception as e:
            print(f"CRITICAL: Connection Error: {e}"); exit()

    def is_truly_african(self, raw_content, country_name):
        text = raw_content.decode('utf-8', errors='ignore').lower() if isinstance(raw_content, bytes) else str(raw_content).lower()
        country = country_name.lower()
        
        # 1. Traveler/Import Check
        if re.search(rf"(travel|imported|returned from|hospitalized in).{{0,40}}{country}", text):
            return False, None

        # 2. Specific Geographic Location Tag Check
        geo_tags = ["geographic location", "country", "geo_loc_name", "submerged in"]
        for tag in geo_tags:
            if tag in text and country in text:
                return True, country_name.capitalize()
                
        return False, None

    def run(self):
        print("\n" + "="*60)
        print("    v7.1 PERSISTENT AFRICAN HARVESTER (DB SELECTOR)    ")
        print("="*60)

        # --- STEP 1: DATABASE SELECTION ---
        print("\n[CONFIG] Choose databases to include in this run:")
        use_bp = input("Search BioProject? (y/n): ").lower() == 'y'
        use_bs = input("Search BioSample? (y/n): ").lower() == 'y'
        use_sra = input("Search SRA? (y/n): ").lower() == 'y'
        use_nuc = input("Search Nucleotide? (y/n): ").lower() == 'y'

        dbs = []
        if use_bp: dbs.append("bioproject")
        if use_bs: dbs.append("biosample")
        if use_sra: dbs.append("sra")
        if use_nuc: dbs.append("nucleotide")

        if not dbs:
            print("❌ Error: No databases selected. Exiting."); return

        # --- STEP 2: SEARCH PARAMETERS ---
        target_org = input("\nEnter Organism (e.g., Vibrio cholerae): ").strip()
        target_goal = int(input("How many NEW verified records do you want PER COUNTRY? (e.g., 20): ") or 20)

        print("\n⏳ Fetching existing accessions to prevent duplicates...")
        existing_accs = set(self.master.col_values(2))

        for country in constants.AFRICA:
            print(f"\n🌍 SCANNING: {country.upper()}")
            
            for db in dbs:
                found_for_this_db = 0
                # Search for a large pool (500) so we have enough to filter through
                query = f'("{target_org}"[Organism] OR "{target_org}") AND "{country}"'
                
                try:
                    search_handle = Entrez.esearch(db=db, term=query, retmax=500)
                    ids = Entrez.read(search_handle)["IdList"]
                    
                    if not ids:
                        continue

                    print(f"  └─ {db.upper()}: Scanning pool of {len(ids)} to find {target_goal} African samples...")

                    for acc in ids:
                        # Goal check
                        if found_for_this_db >= target_goal:
                            print(f"  └─ Goal reached for {db.upper()}.")
                            break
                        
                        # Duplicate check
                        if acc in existing_accs:
                            continue

                        # Persistent Fetch & Scan
                        try:
                            time.sleep(0.4) # Respect NCBI API
                            handle = Entrez.efetch(db=db, id=acc, rettype="xml" if db != "nucleotide" else "gb", retmode="text")
                            raw_data = handle.read()
                            
                            is_valid, verified_country = self.is_truly_african(raw_data, country)
                            
                            if is_valid:
                                main_url = f"https://www.ncbi.nlm.nih.gov/{db}/{acc}"
                                if db == "sra": main_url = f"https://trace.ncbi.nlm.nih.gov/Traces/sra/?run={acc}"
                                
                                row = [
                                    f"{verified_country[:2].upper()}_{acc}", 
                                    acc, target_org, db.capitalize(), 
                                    "N/A", "Verified", main_url, "N/A", 
                                    "v7.1_Persistent", "N/A", "N/A", "N/A", 
                                    datetime.now().strftime("%Y-%m-%d"), 
                                    "OmniBot_v7.1", verified_country
                                ]
                                
                                # Immediate write to sheet
                                self.master.append_row(row)
                                existing_accs.add(acc)
                                found_for_this_db += 1
                                print(f"     [+] Verified & Saved: {acc}")

                        except Exception:
                            continue # Move to next ID if fetch fails

                except Exception as e:
                    print(f"  └─ {db.upper()} Error: {e}")
                
        print("\n✓ MISSION COMPLETE. All selected databases processed.")

if __name__ == "__main__":
    PersistentAfricanHarvester().run()
