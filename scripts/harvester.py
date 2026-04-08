# harvester.py
import time
import os
import re
import csv
import threading
from queue import Queue
from Bio import Entrez
import config
import constants
from huggingface_hub import HfApi
import os

Entrez.email = config.NCBI_EMAIL


class ParallelOmniSystem:
    def __init__(self):
        print("🚀 Initializing Open Source Parallel Pipeline...")

        self.existing_accs = set()

        # 1. Initialize Local Backup File & Load Existing Accessions
        file_exists = os.path.exists(config.LOCAL_BACKUP)

        if file_exists:
            print(f"📊 Loading existing accessions from {config.LOCAL_BACKUP} to prevent duplicates...")
            with open(config.LOCAL_BACKUP, "r", encoding="utf-8") as f:
                reader = csv.reader(f)
                next(reader, None)  # Skip header
                for row in reader:
                    if len(row) > 1:
                        self.existing_accs.add(row[1])  # Accession is column 2
        else:
            print(f"🆕 Creating new master database: {config.LOCAL_BACKUP}")
            with open(config.LOCAL_BACKUP, "w", newline='', encoding="utf-8") as f:
                writer = csv.writer(f)
                writer.writerow(
                    ["Unique_ID", "Accession", "Organism", "Type", "Length", "File_Reference", "Source_URL", "Country"])

        # Queues for Parallel Processing
        self.search_to_harvest_q = Queue(maxsize=100)
        self.harvest_to_fasta_q = Queue(maxsize=100)

    def _save_local(self, row):
        """Writes to local CSV immediately after every harvest."""
        with open(config.LOCAL_BACKUP, "a", newline='', encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow(row)

    def producer_search(self):
        """Thread 1: Searches NCBI and finds candidates based on user selections."""
        print("\n" + "=" * 30 + "\nSEARCH PARAMETERS\n" + "=" * 30)

        use_bp = input("Search BioProject? (y/n): ").lower() == 'y'
        use_bs = input("Search BioSample? (y/n): ").lower() == 'y'
        use_sra = input("Search SRA? (y/n): ").lower() == 'y'
        use_nuc = input("Search Nucleotide? (y/n): ").lower() == 'y'

        selected_dbs = []
        if use_bp: selected_dbs.append("bioproject")
        if use_bs: selected_dbs.append("biosample")
        if use_sra: selected_dbs.append("sra")
        if use_nuc: selected_dbs.append("nucleotide")

        if not selected_dbs:
            print("⚠️ No databases selected. Exiting.");
            self.search_to_harvest_q.put(None);
            return

        target_org = input("Enter Organism (e.g., Vibrio cholerae): ").strip()

        print("\nSelect Sequence Types to include in query:")
        types = []
        if input("Include Plasmids? (y/n): ").lower() == 'y': types.append("plasmid")
        if input("Include Whole Genome (WGS)? (y/n): ").lower() == 'y': types.append("WGS")
        if input("Include Specific Genes? (y/n): ").lower() == 'y': types.append("gene")

        type_query = " OR ".join([f'"{t}"' for t in types]) if types else ""
        target_goal = int(input("Target new records per country: ") or 10)

        for country in constants.AFRICA:
            print(f"🔎 [Searcher] Scanning {country}...")
            for db in selected_dbs:
                query = f'("{target_org}"[Organism]) AND "{country}"'
                if type_query: query += f' AND ({type_query})'

                try:
                    search_handle = Entrez.esearch(db=db, term=query, retmax=100)
                    ids = Entrez.read(search_handle)["IdList"]

                    found_for_country = 0
                    for acc in ids:
                        if acc in self.existing_accs: continue
                        if found_for_country >= target_goal: break

                        # Strict African Verification (Anti-Traveler)
                        time.sleep(0.3)
                        v_handle = Entrez.efetch(db=db, id=acc, rettype="gb" if db == "nucleotide" else "xml",
                                                 retmode="text")
                        raw_data = v_handle.read()

                        # --- NEW FIX: Convert bytes to text if necessary ---
                        if isinstance(raw_data, bytes):
                            raw_data = raw_data.decode('utf-8', errors='ignore')

                        raw_data = raw_data.lower()
                        # -------------------------------------------------

                        if country.lower() in raw_data and not re.search(r"travel|imported|returned from", raw_data):
                            data = {
                                "uid": f"{country[:2].upper()}_{acc}",
                                "acc": acc,
                                "org": target_org,
                                "country": country.capitalize(),
                                "db": db
                            }
                            self.search_to_harvest_q.put(data)
                            found_for_country += 1
                            self.existing_accs.add(acc)
                            print(f"✅ [Searcher] Queued {acc} ({db}) in {country}")

                except Exception as e:
                    print(f"❌ [Searcher] Error in {db} for {country}: {e}")

        self.search_to_harvest_q.put(None)  # End of search signal

    def worker_harvest(self):
        """Thread 2: Downloads DNA, updates local CSV, passes DNA to Writer."""
        while True:
            item = self.search_to_harvest_q.get()
            if item is None: break

            acc = item['acc']
            db_type = item['db']
            print(f"📡 [Harvester] Extracting DNA for {acc}...")

            if db_type == "sra":
                main_url = f"https://trace.ncbi.nlm.nih.gov/Traces/sra/?run={acc}"
            elif db_type == "nucleotide":
                main_url = f"https://www.ncbi.nlm.nih.gov/nuccore/{acc}"
            else:
                main_url = f"https://www.ncbi.nlm.nih.gov/{db_type}/{acc}"

            try:
                time.sleep(0.5)
                h_gb = Entrez.efetch(db="nucleotide", id=acc, rettype="gb", retmode="text")
                raw_gb = h_gb.read().lower()
                stype = "Plasmid" if "plasmid" in raw_gb else "WGS" if "wgs" in raw_gb else "Genomic DNA"

                time.sleep(0.5)
                h_fa = Entrez.efetch(db="nucleotide", id=acc, rettype="fasta", retmode="text")
                raw_fa = h_fa.read().strip()
                dna = "\n".join(raw_fa.split('\n')[1:]) if raw_fa.startswith('>') else raw_fa

                if dna:
                    dna_length = str(len(dna))

                    row = [
                        item['uid'],
                        acc,
                        item['org'],
                        stype,
                        dna_length,
                        config.SEQUENCE_FILE,
                        main_url,
                        item['country']
                    ]

                    self._save_local(row)
                    print(f"🟢 [Harvester] Logged {acc} to local CSV Database!")

                    item['dna'] = dna
                    item['stype'] = stype
                    self.harvest_to_fasta_q.put(item)

            except Exception as e:
                print(f"⚠️ [Harvester] Failed {acc}: {e}")

            self.search_to_harvest_q.task_done()

        self.harvest_to_fasta_q.put(None)

    def consumer_fasta_writer(self):
        """Thread 3: Buffers DNA and batch writes to FASTA file."""
        fasta_buffer = []

        while True:
            item = self.harvest_to_fasta_q.get()

            if item is None:
                if fasta_buffer:
                    self._write_batch_to_file(fasta_buffer)
                self.harvest_to_fasta_q.task_done()
                break

            fasta_buffer.append(item)

            if len(fasta_buffer) >= config.BATCH_SIZE:
                self._write_batch_to_file(fasta_buffer)
                fasta_buffer.clear()

            self.harvest_to_fasta_q.task_done()

    def _write_batch_to_file(self, buffer_array):
        print(f"✍️ [Writer] Memory full! Batch writing {len(buffer_array)} sequences to {config.SEQUENCE_FILE}...")
        try:
            with open(config.SEQUENCE_FILE, "a", encoding="utf-8") as f:
                for item in buffer_array:
                    f.write(f">{item['uid']} | {item['acc']} | {item['org']} | {item['stype']} | {item['country']}\n")
                    f.write(item['dna'] + "\n\n")
            print(f"✅ [Writer] Successfully dumped batch to disk.")
        except Exception as e:
            print(f"❌ [Writer] Error batch writing to file: {e}")

    def run_all_parallel(self):
        start_time = time.time()

        search_thread = threading.Thread(target=self.producer_search)
        harvest_thread = threading.Thread(target=self.worker_harvest)
        writer_thread = threading.Thread(target=self.consumer_fasta_writer)

        search_thread.start()
        harvest_thread.start()
        writer_thread.start()

        search_thread.join()
        harvest_thread.join()
        writer_thread.join()

        # --- NEW HUGGING FACE UPLOAD STEP ---
        print("\n🚀 Uploading massive FASTA file to Hugging Face...")
        try:
            api = HfApi()
            api.upload_file(
                path_or_fileobj=config.SEQUENCE_FILE, # The local 700MB file
                path_in_repo="plas_kmer_sequences.fasta", # What to call it on Hugging Face
                repo_id="Jeffiq/Plaskmer", # UPDATE THIS to your HF username and dataset
                repo_type="dataset",
                token=os.environ.get("HF_TOKEN") # Safely grabs the password from GitHub Secrets
            )
            print("✅ Successfully pushed to Hugging Face!")
        except Exception as e:
            print(f"❌ Hugging Face upload failed: {e}")
        # ------------------------------------

        total_time = round((time.time() - start_time) / 60, 2)
        print(f"\n{'=' * 40}\n✨ MISSION COMPLETE\nTotal Time: {total_time} minutes\n{'=' * 40}")


if __name__ == "__main__":
    bot = ParallelOmniSystem()
    bot.run_all_parallel()
