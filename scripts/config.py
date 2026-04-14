# config.py
# ---------------------------------------------------------
# Configuration settings for the Afrigen Harvester
# ---------------------------------------------------------

NCBI_EMAIL = "methajefferson@gmail.com"  # Your NCBI email
MASTER_PARQUET = "master_database.parquet" # THE NEW BIG DATA FILE
BATCH_SIZE = 50  # Number of sequences to hold in RAM before writing to Parquet
HF_REPO_ID = "Jeffiq/Plaskmer" # Change to your actual HF Repo ID later
