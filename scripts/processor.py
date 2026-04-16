import pandas as pd
import time
from concurrent.futures import ThreadPoolExecutor
from Bio.Seq import Seq
from Bio import Entrez
import config
import streamlit as st

if 'user_email' in st.session_state and st.session_state.user_email:
    from Bio import Entrez
    Entrez.email = st.session_state.user_email


    
class GeneticHunter:
    """Parallel Engine to hunt for DNA features and enrich with NCBI PubMed/Taxonomy."""
    
    @staticmethod
    def fetch_ncbi_metadata(accession):
        """Reaches back into NCBI to find linked PubMed Articles and Taxonomy IDs."""
        tax_id = "Pending"
        pubmed_ids = []
        
        try:
            # Respect NCBI rate limits (max 3-10 requests per second)
            time.sleep(0.35) 
            
            # 1. Fetch Taxonomy ID via eSummary
            sum_handle = Entrez.esummary(db="nucleotide", id=accession)
            summary = Entrez.read(sum_handle)
            sum_handle.close()
            if summary and len(summary) > 0:
                tax_id = summary[0].get('TaxId', 'Unknown')
                
            # 2. Fetch linked PubMed Articles via eLink
            link_handle = Entrez.elink(dbfrom="nucleotide", db="pubmed", id=accession)
            links = Entrez.read(link_handle)
            link_handle.close()
            
            if links and links[0].get('LinkSetDb'):
                for link in links[0]['LinkSetDb'][0]['Link']:
                    pubmed_ids.append(link['Id'])
                    
        except Exception as e:
            # If NCBI times out, we just return what we have so we don't crash the pipeline
            pass 
            
        return str(tax_id), ", ".join(pubmed_ids) if pubmed_ids else "None"

    @staticmethod
    def hunt_features(record):
        """Deep analysis of a single sequence record, guaranteeing no data loss."""
        
        # 0. COPY ORIGINAL RECORD TO PREVENT DATA LOSS!
        enriched_data = record.copy()
        
        seq = str(record.get('Sequence', '')).upper() # .upper() ensures motifs match regardless of sequence case
        header = str(record.get('Type', '')).lower()
        accession = record.get('Accession', '')
        
        # 1. Offline Hunt: Refine Type using Biology Heuristics
        refined_type = "WGS"
        if "plasmid" in header or (0 < len(seq) < 200000): # Plasmids usually < 200kbp
            refined_type = "Plasmid"
        elif any(x in header for x in ["mrna", "transcript", "cdna"]):
            refined_type = "mRNA"
        elif len(seq) > 500000:
            refined_type = "Chromosome"
            
        # 2. Offline Hunt: Environmental & Clinical Markers
        markers = []
        marker_map = {
            # Clinical Antibiotic Resistance
            "AmpC (Ampicillin)": "ACTG", # Placeholder sequence motif
            "NDM-1 (Carbapenem)": "ATGCGT", 
            "mcr-1 (Colistin)": "CGTAGC", 
            "tet(X) (Tetracycline)": "CCAA",
            
            # Environmental Heavy Metal Resistance
            "merA (Mercury)": "GGCATT",
            "arsC (Arsenic)": "TTAACG",
            "copA (Copper)": "GGACCC"
        } 
        
        for name, motif in marker_map.items():
            if motif in seq: markers.append(name)

        # 3. Online Hunt: Go to NCBI for Metadata
        tax_id, pubmed_refs = GeneticHunter.fetch_ncbi_metadata(accession)

        # 4. Inject new intelligence into the record
        enriched_data['Type'] = refined_type
        enriched_data['Length'] = len(seq) # Ensure length is accurate
        enriched_data['Selection_Marker'] = ", ".join(markers) if markers else "None"
        enriched_data['Taxonomy_ID'] = tax_id
        enriched_data['PubMed_IDs'] = pubmed_refs
        enriched_data['Date_Harvested'] = pd.Timestamp.now().strftime('%Y-%m-%d')
        
        return enriched_data

def parallel_process_records(records):
    """Entry point for ThreadPool processing (Optimized for Web Requests)."""
    # Using ThreadPool instead of ProcessPool because we are doing web requests (I/O bound)
    with ThreadPoolExecutor(max_workers=5) as executor:
        results = list(executor.map(GeneticHunter.hunt_features, records))
    return results
