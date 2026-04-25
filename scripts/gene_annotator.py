# gene_annotator.py
from Bio import Entrez
from Bio import SeqIO
import pandas as pd
import streamlit as st
import os
import time

# We DO NOT set Entrez.email here at the top. It crashes the app.


def fetch_official_genes(accession_id):
    """
    Connects to NCBI, downloads the GenBank record for the given ID,
    and extracts the official gene names and coordinates.
    """
    try:
        # 0. Safely grab the email and API key right before we make the call
        if 'user_email' in st.session_state and st.session_state.user_email:
            Entrez.email = st.session_state.user_email
        else:
            Entrez.email = "default_plaskmer@example.com" # Fallback just in case
            
        if os.getenv("NCBI_API_KEY"):
            Entrez.api_key = os.getenv("NCBI_API_KEY")

        # 1. Fetch the GenBank (.gb) file from NCBI
        time.sleep(0.35) # Respect NCBI 429 limits
            
        handle = Entrez.efetch(db="nucleotide", id=accession_id, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        handle.close()
        
        gene_list = []
        
        # 2. Loop through the features of the sequence
        for feature in record.features:
            # We are looking for "CDS" (Coding Sequences - which are the actual protein genes)
            if feature.type == "CDS":
                # Extract data safely (some genes might miss a label, so we use .get)
                gene_name = feature.qualifiers.get("gene", ["Unknown Gene"])[0]
                product = feature.qualifiers.get("product", ["Unknown Product"])[0]
                start = int(feature.location.start)
                end = int(feature.location.end)
                strand = feature.location.strand
                
                gene_list.append({
                    "Gene Name": gene_name,
                    "Product / Protein": product,
                    "Start (bp)": start,
                    "End (bp)": end,
                    "Strand": "+" if strand == 1 else "-"
                })
                
        return pd.DataFrame(gene_list)
        
    except Exception as e:
        return str(e) # Return the error message if it fails
