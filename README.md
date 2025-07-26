#  URTC-ACHE-PIPELINE

[![Python](https://img.shields.io/badge/python-3.8%2B-blue)](https://www.python.org/)  
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

An end‑to‑end **ADMET → ML → DFT** pipeline for discovering and ranking novel acetylcholinesterase (AChE) inhibitor scaffolds.

---

## Overview

1. **01_chembl_retrieval.py** – Fetch ChEMBL AChE inhibitors (IC₅₀ ≤ 100 nM)  
2. **02_admet_filter.py** – Apply strict ADMET + loose PK filters  
3. **03_machine_learning.py** – UMAP embedding & K‑means medoid clustering  
4. **04_toxicity_scoring.py** – Composite toxicity scoring (hERG, DILI, Ames, SkinSens)  
5. **05_generate_orca_inputs.py** – Generate ORCA v5.0 input files for medoids  
6. **06_parse_dft.py** – Parse ORCA outputs for HOMO–LUMO gaps & normalize  
7. **07_compute_final_ranking.py** – Combine ADMET‑distance & sDFT into final ranking  

---

## Quick Start

```bash
git clone https://github.com/abdullahchdry/urtc-ache-pipeline.git
cd urtc-ache-pipeline
pip install -r requirements.txt

```

---

## Clean (recommended)
Before you start, you can remove any prior outputs to guarantee a fresh run:
```bash
rm data/processed/* figures/*
```
## Quick start (contd.)
```bash
python scripts/01_chembl_retrieval.py
python scripts/02_admet_filter.py
python scripts/03_umap_clustering.py
python scripts/04_toxicity_scoring.py
python scripts/05_generate_orca_inputs.py
# submit ORCA jobs, place resulting .out files into stability-io/output/
python scripts/06_parse_dft_gaps.py
python scripts/07_compute_final_ranking.py
```


---
## Dependencies

- Python 3.8+
- RDKit
- pandas
- numpy
- scikit-learn
- umap-learn
- matplotlib
- chembl_webresource_client
- ORCA v5.0

