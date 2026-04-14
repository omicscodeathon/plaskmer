# orf_blaster.py
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import pandas as pd
import time

def identify_orf_via_blast(protein_sequence):
    """
    Takes a single translated protein sequence and BLASTs it against NCBI's database.
    Returns the top hit name.
    """
    try:
        # Run BLASTP (Protein BLAST) against the non-redundant database
        result_handle = NCBIWWW.qblast("blastp", "nr", protein_sequence, hitlist_size=1)
        blast_record = NCBIXML.read(result_handle)
        
        if blast_record.alignments:
            top_alignment = blast_record.alignments[0]
            # Clean up the long NCBI title to get just the protein name
            hit_def = top_alignment.hit_def.split(">")[0] 
            return hit_def
        else:
            return "No match found"
            
    except Exception as e:
        return f"Error: {str(e)}"
