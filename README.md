---
title: Plaskmer
emoji: 🧬
colorFrom: blue
colorTo: green
sdk: streamlit
app_file: scripts/app.py
pinned: false
---

# Plaskmer
(You can put the rest of your normal README text down here...)
# plaskmer
Building an African Plasmid and mRNA Databank with K-mer-Based Search
# Overview
Africa remains significantly underrepresented in global microbial genomic databases, limiting effective surveillance of infectious diseases and antimicrobial resistance (AMR). This gap is particularly evident in high-burden diseases such as cholera, caused by Vibrio cholerae, which continues to drive recurrent outbreaks across the continent.
A major driver of Vibrio cholerae evolution is horizontal gene transfer via plasmids and other mobile genetic elements, enabling rapid dissemination of virulence and antimicrobial resistance determinants. However, conventional homology-based tools often fail to detect novel or highly divergent sequences, leaving substantial genomic diversity uncharacterized.
Plaskmer is an alignment-free, k-mer based bioinformatics framework designed to address this limitation by enabling reference-independent classification of plasmid and chromosomal sequences. It builds an Africa-focused plasmid and mRNA database starting with Vibrio cholerae, supporting improved genomic surveillance and AMR tracking.


---

## 🌟 Significance & Impact
Traditional genomic classification relies on homology-based tools (like BLAST), which struggle with "genomic dark matter"—novel sequences that lack close relatives in existing databases. 

Plaskmer addresses this by focusing on the **intrinsic statistical signatures** of DNA:
* **African Pathogen Surveillance:** Rapidly sort plasmids from high-burden African pathogens (e.g., *Vibrio cholerae*, *K. pneumoniae*).
* **Environmental Exploration:** Identify mobile genetic elements (MGEs) in African extreme environments, such as soda lakes and thermal pans.
* **Functional Genomics:** Track Open Reading Frames (ORFs) to understand how functional traits are shared between different bacterial species via plasmids.

---

## 🛠️ Methods & Innovation
Plaskmer utilizes a normalized $k$-mer frequency distribution model ($k=3$ to $k=6$) to generate a unique "genomic fingerprint" for any input sequence.

### The Logic
1. **Automated Harvesting:** Uses the NCBI Entrez API to fetch genomic and plasmid data based on TaxID or Species name.
2. **Signature Calculation:** Transforms raw FastA sequences into probability distributions:
   $$P(i) = \frac{c_i}{\sum_{j=1}^{4^k} c_j}$$
3. **Dimensionality Reduction:** Employs **Principal Component Analysis (PCA)** to project high-dimensional $k$-mer signatures into visual clusters.
4. **Classification:** Measures the Euclidean distance between an unknown ORF and the host chromosomal baseline to determine genomic origin.

---

## 🚀 Quick Start (Hugging Face)
Access the live dashboard here: **[Plaskmer on Hugging Face](https://jeffiq-plaskmer.hf.space/)**

1. **Harvester Tab:** Fetch data directly from NCBI to populate the local databank.
2. **Analysis Tab:** Visualize the K-mer frequency differences between chromosomal and plasmid DNA.
3. **ORF Scanner:** Predict and classify the origin of specific functional sequences.

---

## 📊 Benchmarking & Innovation
* **Dockerized Deployment:** Unlike standard Streamlit apps, Plaskmer runs in a custom Docker container with **16GB RAM**, allowing for the processing of large-scale environmental genomic datasets.
* **Alignment-Free:** Faster than BLAST/Diamond, enabling preliminary screening of uncharacterized African biodiversity.
