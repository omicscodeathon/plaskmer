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
from processor import parallel_process_records # <--- Moved to the top!

# Ensure your NCBI email is set in your config.py

import streamlit as st # Add this at the top

# NEW LOGIC: Pull from the user's browser session, not your config
user_email = st.session_state.get('user_email')
Entrez.api_key = os.getenv("NCBI_API_KEY")
if user_email:
    Entrez.email = user_email

class ParallelOmniSystem:
    def __init__(self, target_org, selected_dbs, types, target_goal, push_to_hf=False, db_path=None):
        print("🚀 Initializing Parquet-Powered Parallel Pipeline...")
        
        self.target_org = target_org
        self.selected_dbs = selected_dbs
        self.types = types
        self.target_goal = target_goal
        self.push_to_hf = push_to_hf
        self.existing_accs = set()

        # 1. Initialize Local Parquet Backup & Load Existing Accessions
        # 🚨 USE THE CACHED FILE APP.PY DOWNLOADED FROM HF, OTHERWISE FALLBACK
        self.parquet_file = db_path if db_path else config.MASTER_PARQUET
        
        print(f"\n🎯 ALERTTTT! I AM SAVING THE FILE EXACTLY HERE: {os.path.abspath(self.parquet_file)}\n")
        
        if os.path.exists(self.parquet_file):
            print(f"📊 Loading existing database from {self.parquet_file}...")
            try:
                # We only load the Accession column to save RAM!
                df_existing = pd.read_parquet(self.parquet_file, columns=["Accession"])
                self.existing_accs.update(df_existing["Accession"].dropna().tolist())
            except Exception as e:
                print(f"⚠️ Could not read existing Parquet file. It might be corrupted from a previous crash: {e}")
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

        self.search_to_harvest_q.put(None)

    def worker_harvest(self):
        """Thread 2: Consumes queued accessions, fetches FASTA sequences."""
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
                            "Country": item["Country"],
                            "db": db # Passing this so the processor keeps it
                        }
                        
                        self.harvest_to_writer_q.put(record)
                        print(f"🧬 [Harvester] Downloaded sequence for {acc} ({seq_length} bp)")
                        
            except Exception as e:
                print(f"❌ [Harvester] Failed to fetch data for {acc}: {e}")
                
            self.search_to_harvest_q.task_done()

    def consumer_parquet_writer(self):
        """Listens to the queue, writes to Parquet in batches, and syncs to HF to save RAM."""
        print("💾 Parquet Writer Thread Active.")
        import gc # Garbage collector to manually flush RAM
        
        batch_records = []
        BATCH_SIZE = 25  # 🚨 Save, Sync, and Flush RAM every 25 records
        
        while True:
            item = self.harvest_to_writer_q.get()
            
            # If we receive 'None', the harvest is completely finished
            if item is None:
                if batch_records:
                    self._save_batch_to_disk(batch_records)
                    if self.push_to_hf:
                        self.push_database_to_hf()
                self.harvest_to_writer_q.task_done()
                break
                
            # Add the new record to our temporary batch
            batch_records.append(item)
            
            # 🚨 THE MAGIC: IF BATCH IS FULL, SAVE, SYNC, AND DUMP MEMORY
            if len(batch_records) >= BATCH_SIZE:
                print(f"📦 Batch of {BATCH_SIZE} reached! Saving to disk and clearing memory...")
                
                # 1. Save to local Parquet
                self._save_batch_to_disk(batch_records)
                
                # 2. Push the checkpoint to Hugging Face
                if self.push_to_hf:
                    print("⬆️ Pushing mid-harvest checkpoint to Hugging Face...")
                    self.push_database_to_hf()
                
                # 3. 🧹 COMPLETELY FLUSH THE RAM
                batch_records.clear()
                gc.collect() 
                
            self.harvest_to_writer_q.task_done()

    def _save_batch_to_disk(self, batch):
        """Helper to safely append a batch to the Parquet file."""
        df_batch = pd.DataFrame(batch)
        if os.path.exists(self.parquet_file):
            # Read existing, append, and overwrite
            df_existing = pd.read_parquet(self.parquet_file)
            df_combined = pd.concat([df_existing, df_batch], ignore_index=True)
            df_combined.to_parquet(self.parquet_file, index=False)
        else:
            # Create new if it doesn't exist
            df_batch.to_parquet(self.parquet_file, index=False)

    def _write_batch_to_parquet(self, batch):
        """Helper function to append a batch of records to the Parquet file safely."""
        try:
            # 🔥 RUN THE HUNTER FIRST
            print(f"🔎 [Hunter] Deep-scanning and contacting NCBI for {len(batch)} sequences...")
            processed_batch = parallel_process_records(batch)
            
            df_batch = pd.DataFrame(processed_batch)
            
            # 🚨 ABSOLUTE SCHEMA COMPLIANCE
            df_batch['Unique_ID'] = df_batch.get('Unique_ID', '').astype("string")
            df_batch['Accession'] = df_batch.get('Accession', '').astype("string")
            df_batch['Organism'] = df_batch.get('Organism', '').astype("string")
            df_batch['Type'] = df_batch.get('Type', 'WGS').astype("string")
            df_batch['Length'] = df_batch.get('Length', 0).astype("int32")
            df_batch['Sequence'] = df_batch.get('Sequence', '').astype("string")
            df_batch['Source_URL'] = df_batch.get('Source_URL', '').astype("string")
            df_batch['Country'] = df_batch.get('Country', '').astype("string")
            df_batch['db'] = df_batch.get('db', 'unknown').astype("string")
            
            # New Enriched Fields
            df_batch['Selection_Marker'] = df_batch.get('Selection_Marker', 'None').astype("string")
            df_batch['Taxonomy_ID'] = df_batch.get('Taxonomy_ID', 'Pending').astype("string")
            df_batch['PubMed_IDs'] = df_batch.get('PubMed_IDs', 'None').astype("string")
            df_batch['Date_Harvested'] = df_batch.get('Date_Harvested', '').astype("string")
            
            # ENGINE SWITCH: 'pyarrow' handles schema evolution (new columns) better than 'fastparquet'
            if os.path.exists(self.parquet_file):
                # We read the old file, concat the new batch, and overwrite to ensure schema matches
                df_existing = pd.read_parquet(self.parquet_file, engine='pyarrow')
                df_combined = pd.concat([df_existing, df_batch], ignore_index=True)
                df_combined.to_parquet(self.parquet_file, engine='pyarrow', index=False)
            else:
                df_batch.to_parquet(self.parquet_file, engine='pyarrow', index=False)
                
            print(f"💾 [Writer] Successfully saved enriched batch of {len(batch)} records to Parquet.")
        except Exception as e:
            print(f"❌ [Writer] Failed to process/write batch to Parquet: {e}")

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
                    path_in_repo=f"data/{os.path.basename(self.parquet_file)}",
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
        
        num_workers = 1
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
