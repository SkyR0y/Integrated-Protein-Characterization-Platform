# Integrated-Protein-Characterization-Platform

Integrated Protein Characterization Platform (IPCP) is a Streamlit-based bioinformatics web application designed for rapid in-silico protein analysis. The platform integrates multiple protein characterization utilities into a single user-friendly interface.

## Features

### Protein Analysis
Analyze protein sequences using physicochemical parameters:

- Molecular Weight
- Isoelectric Point (pI)
- Instability Index
- GRAVY (Hydrophobicity)
- Functional Interpretation

### Input Options

- Raw amino acid sequence
- UniProt ID retrieval

### Isoform Analysis

Compare multiple protein isoforms retrieved from UniProt:

- Molecular Weight comparison
- pI comparison
- Stability comparison
- Hydrophobicity comparison
- Graphical visualization

### Subcellular Localization

Integrated Phobius-based prediction for:

- Signal peptide detection
- Transmembrane topology
- Cytoplasmic / Non-cytoplasmic domain prediction

### Reports

- Downloadable protein characterization reports

---

## Technologies Used

- Python
- Streamlit
- Pandas
- Biopython
- Requests
- BeautifulSoup4

---

## Data Sources / Integrations

- UniProt REST API
- ExPASy ProtParam
- Phobius Server

---
