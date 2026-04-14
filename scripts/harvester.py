import time
import os
import re
import threading
from queue import Queue
import pandas as pd
from Bio import Entrez
from huggingface_hub import HfApi
import config
import constants

# Ensure your NCBI email is set in your config.py
Entrez.email = config.NCBI_EMAIL

class ParallelOmniSystem:
    def __init__(self, target_org, selected_dbs, types, target_goal, push_to_hf=False):
        print("🚀 Initializing Parquet-Powered Parallel Pipeline...")
        
        self.target_org = target_org
        self.selected_dbs = selected_dbs
        self.types = types
        self.target_goal = target_goal
        self.push_to_hf = push_to_hf
        self.existing_accs = set()

        # 1. Initialize Local Parquet Backup & Load Existing Accessions
        self.parquet_file = config.MASTER_PARQUET
        
        # 🚨 THE DEBUG ALARM: THIS TELLS US EXACTLY WHERE THE FILE IS GOING
        print(f"\n🎯 ALERTTTT! I AM SAVING THE FILE EXACTLY HERE: {os.path.abspath(self.parquet_file)}\n")
        
        if os.path.exists(self.parquet_file):
            print(f"📊 Loading existing database from {self.parquet_file}...")
            try:
                # We only load the Accession column to save RAM!
                df_existing = pd.read_parquet(self.parquet_file, columns=["Accession"])
                self.existing_accs.update(df_existing["Accession"].dropna().tolist())
            except Exception as e:
                print(f"⚠️ Could not read existing Parquet file. It might be corrupted from a previous crash: {e}")
                # 🚨 Auto-delete corrupted files so we can start fresh!
                try:
                    os.remove(self.parquet_file)
                    print("🗑️ Deleted corrupted file. Starting fresh!")
                except:
                    print("⚠️ Could not delete corrupted file.")
        else:
            print(f"🆕 Master database not found. A new one will be created upon first batch.")

        # Queues for Parallel Processing
        self.search_to_harvest_q = Queue(maxsize=100)
        self.harvest_to_writer_q = Queue(maxsize=100)

    def producer_search(self):
        """Thread 1: Searches NCBI and finds candidates based on UI selections."""
        print(f"\n==============================\nSEARCH PARAMETERS: {self.target_org}\n==============================")

        if not self.selected_dbs:
            print("⚠️ No databases selected. Exiting.")
            self.search_to_harvest_q.put(None)
            return

        type_query = " OR ".join([f'"{t}"' for t in self.types]) if self.types else ""

        for country in constants.AFRICA:
            print(f"🔎 [Searcher] Scanning {country}...")
            for db in self.selected_dbs:
                query = f'("{self.target_org}"[Organism]) AND "{country}"'
                if type_query: query += f' AND ({type_query})'

                try:
                    search_handle = Entrez.esearch(db=db, term=query, retmax=100)
                    ids = Entrez.read(search_handle)["IdList"]

                    found_for_country = 0
                    for acc in ids:
                        # Convert to strict string immediately!
                        acc = str(acc) 
                        
                        if acc in self.existing_accs: continue
                        if found_for_country >= self.target_goal: break

                        time.sleep(0.3)
                        v_handle = Entrez.efetch(db=db, id=acc, rettype="gb" if db == "nucleotide" else "xml", retmode="text")
                        raw_data = v_handle.read()

                        if isinstance(raw_data, bytes):
                            raw_data = raw_data.decode('utf-8', errors='ignore')

                        raw_data = raw_data.lower()

                        if country.lower() in raw_data and not re.search(r"travel|imported|returned from", raw_data):
                            data = {
                                "Unique_ID": f"{country[:2].upper()}_{acc}",
                                "Accession": acc,
                                "Organism": str(self.target_org),
                                "Country": str(country.capitalize()),
                                "Source_URL": f"https://www.ncbi.nlm.nih.gov/nuccore/{acc}",
                                "db": str(db)
                            }
                            self.search_to_harvest_q.put(data)
                            found_for_country += 1
                            self.existing_accs.add(acc)
                            print(f"✅ [Searcher] Queued {acc} ({db}) in {country}")

                except Exception as e:
                    print(f"❌ [Searcher] Error in {db} for {country}: {e}")

        # Send poison pill to close workers
        self.search_to_harvest_q.put(None)

    def worker_harvest(self):
        """Thread 2: Consumes queued accessions, fetches FASTA sequences, and queues them for the Parquet writer."""
        while True:
            item = self.search_to_harvest_q.get()
            
            if item is None:
                self.search_to_harvest_q.put(None)
                break

            acc = item['Accession']
            db = item['db']
            
            try:
                time.sleep(0.4)
                
                if db == "sra":
                    record = {
                        "Unique_ID": item["Unique_ID"],
                        "Accession": acc,
                        "Organism": item["Organism"],
                        "Type": "SRA Raw Reads",
                        "Length": 0,
                        "Sequence": "",  
                        "Source_URL": item["Source_URL"],
                        "Country": item["Country"]
                    }
                    self.harvest_to_writer_q.put(record)
                    print(f"🧬 [Harvester] Logged SRA run {acc} (No assembled sequence)")
                    
                else:
                    handle = Entrez.efetch(db="nucleotide", id=acc, rettype="fasta", retmode="text")
                    fasta_data = handle.read().strip()
                    handle.close()
                    
                    if fasta_data:
                        lines = fasta_data.split('\n')
                        header = lines[0] if len(lines) > 0 else ""
                        sequence = "".join(lines[1:])
                        seq_length = len(sequence)
                        
                        seq_type = "WGS"
                        if "plasmid" in header.lower():
                            seq_type = "Plasmid"
                        elif "gene" in header.lower() or "mrna" in header.lower():
                            seq_type = "mRNA"
                            
                        record = {
                            "Unique_ID": item["Unique_ID"],
                            "Accession": acc,
                            "Organism": item["Organism"],
                            "Type": seq_type,
                            "Length": seq_length,
                            "Sequence": sequence,
                            "Source_URL": item["Source_URL"],
                            "Country": item["Country"]
                        }
                        
                        self.harvest_to_writer_q.put(record)
                        print(f"🧬 [Harvester] Downloaded sequence for {acc} ({seq_length} bp)")
                        
            except Exception as e:
                print(f"❌ [Harvester] Failed to fetch data for {acc}: {e}")
                
            self.search_to_harvest_q.task_done()

    def consumer_parquet_writer(self):
        """Thread 3: Listens for completed records and writes them to Parquet in batches."""
        batch = []
        
        while True:
            record = self.harvest_to_writer_q.get()
            
            if record is None:
                if batch:
                    self._write_batch_to_parquet(batch)
                break
                
            batch.append(record)
            
            if len(batch) >= getattr(config, 'BATCH_SIZE', 50):
                self._write_batch_to_parquet(batch)
                batch = [] 
                
            self.harvest_to_writer_q.task_done()

    def _write_batch_to_parquet(self, batch):
        """Helper function to append a batch of records to the Parquet file safely."""
        try:
            df_batch = pd.DataFrame(batch)
            
            # 🚨 THE FINAL FIX: Force Pandas native String format to absolutely guarantee compliance
            df_batch['Unique_ID'] = df_batch['Unique_ID'].astype("string")
            df_batch['Accession'] = df_batch['Accession'].astype("string")
            df_batch['Organism'] = df_batch['Organism'].astype("string")
            df_batch['Type'] = df_batch['Type'].astype("string")
            df_batch['Length'] = df_batch['Length'].astype("int32")
            df_batch['Sequence'] = df_batch['Sequence'].astype("string")
            df_batch['Source_URL'] = df_batch['Source_URL'].astype("string")
            df_batch['Country'] = df_batch['Country'].astype("string")
            
            if os.path.exists(self.parquet_file):
                df_batch.to_parquet(self.parquet_file, engine='fastparquet', append=True)
            else:
                df_batch.to_parquet(self.parquet_file, engine='fastparquet')
                
            print(f"💾 [Writer] Successfully saved batch of {len(batch)} records to Parquet.")
        except Exception as e:
            print(f"❌ [Writer] Failed to write batch to Parquet: {e}")

    def _upload_to_huggingface(self):
        """Uploads the unified Parquet file to Hugging Face."""
        print("☁️ [Cloud] Initiating Hugging Face Upload...")
        try:
            api = HfApi()
            repo_id = getattr(config, 'HF_REPO_ID', None)
            
            if not repo_id or repo_id == "YourUsername/PlasKmer-Database":
                print("⚠️ [Cloud] Valid HF_REPO_ID not configured. Skipping upload.")
                return

            if os.path.exists(self.parquet_file):
                api.upload_file(
                    path_or_fileobj=self.parquet_file,
                    path_in_repo=f"data/{self.parquet_file}",
                    repo_id=repo_id,
                    repo_type="dataset"
                )
                print("✅ [Cloud] Master Parquet Database successfully synced to Hugging Face!")
            else:
                print("⚠️ [Cloud] No Parquet file found to upload.")
        except Exception as e:
            print(f"❌ [Cloud] Hugging Face sync failed: {e}")

    def run_all_parallel(self):
        """Orchestrates all threads."""
        print("🚦 Starting Parallel Harvest (Parquet Engine)...")
        
        t_writer = threading.Thread(target=self.consumer_parquet_writer)
        t_writer.start()
        
        num_workers = 3
        workers = []
        for _ in range(num_workers):
            t = threading.Thread(target=self.worker_harvest)
            t.start()
            workers.append(t)
            
        self.producer_search()
        
        for t in workers:
            t.join()
            
        self.harvest_to_writer_q.put(None)
        t_writer.join()
        
        if self.push_to_hf:
            self._upload_to_huggingface()

        print("\n🎉 All parallel tasks completed successfully.")
