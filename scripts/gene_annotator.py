# gene_annotator.py
from Bio import Entrez
from Bio import SeqIO
import pandas as pd
import streamlit as st

# VERY IMPORTANT: NCBI requires an email address to use their API
Entrez.email = "methajefferson@gmail.com"  # <--- CHANGE THIS TO YOUR EMAIL

def fetch_official_genes(accession_id):
    """
    Connects to NCBI, downloads the GenBank record for the given ID,
    and extracts the official gene names and coordinates.
    """
    try:
        # 1. Fetch the GenBank (.gb) file from NCBI
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
