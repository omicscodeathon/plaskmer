# config.py
import os
from pathlib import Path

# --- DYNAMIC PATHING ---
# This finds the 'data' folder regardless of local or cloud environment
BASE_DIR = Path(__file__).parent.parent.resolve()
DATA_DIR = BASE_DIR / "data"

# Ensure the data directory exists
if not DATA_DIR.exists():
    os.makedirs(DATA_DIR, exist_ok=True)

NCBI_EMAIL = "methajefferson@gmail.com"
SEQUENCE_FILE = str(DATA_DIR / "plas_kmer_sequences.fasta")
LOCAL_BACKUP = str(DATA_DIR / "local_harvest_log.csv")
BATCH_SIZE = 50

# Optional: If you use Hugging Face, add your repo ID here
HF_REPO_ID = "Jeffiq/Plaskmer"
