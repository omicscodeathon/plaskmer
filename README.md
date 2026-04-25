---
title: Plaskmer
emoji: 🧬
colorFrom: blue
colorTo: green
app_file: scripts/app.py
pinned: false
---
# plaskmer
Building an African Plasmid and mRNA Databank with K-mer-Based Search
# Overview
Plaskmer is a bioinformatics framework designed to facilitate the creation of a specialized genomic databank focused on African microbial diversity. By combining automated data harvesting from NCBI with alignment-free K-mer frequency analysis, Plaskmer allows researchers to rapidly identify and categorize genetic elements (Plasmids and mRNA/ORFs) from environmental and clinical samples across the continent.

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
Access the live dashboard here: **[Plaskmer on Hugging Face](https://huggingface.co/spaces/Jeffiq/Plaskmer)**

1. **Harvester Tab:** Fetch data directly from NCBI to populate the local databank.
2. **Analysis Tab:** Visualize the K-mer frequency differences between chromosomal and plasmid DNA.
3. **ORF Scanner:** Predict and classify the origin of specific functional sequences.

---

## 📊 Benchmarking & Innovation
* **Dockerized Deployment:** Unlike standard Streamlit apps, Plaskmer runs in a custom Docker container with **16GB RAM**, allowing for the processing of large-scale environmental genomic datasets.
* **Alignment-Free:** Faster than BLAST/Diamond, enabling preliminary screening of uncharacterized African biodiversity.
