# MIT_URTC_ACHE_PIPELINE

[![Python](https://img.shields.io/badge/python-3.8%2B-blue)](https://www.python.org/)  
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

An end‑to‑end **ADMET → ML → DFT** pipeline for discovering and ranking novel acetylcholinesterase (AChE) inhibitor scaffolds.

---

## Overview

1. **01_chembl_retrieval.py** – Fetch ChEMBL AChE inhibitors (IC₅₀ ≤ 100 nM)  
2. **02_admet_filter.py** – Apply strict ADMET + loose PK filters  
3. **03_umap_clustering.py** – UMAP embedding & K‑means medoid clustering  
4. **04_toxicity_scoring.py** – Composite toxicity scoring (hERG, DILI, Ames, SkinSens)  
5. **05_generate_orca_inputs.py** – Generate ORCA v5.0 input files for medoids  
6. **06_parse_dft_gaps.py** – Parse ORCA outputs for HOMO–LUMO gaps & normalize  
7. **07_compute_final_ranking.py** – Combine ADMET‑distance & sDFT into final ranking  

---

## Quick Start

```bash
git clone https://github.com/abdullahchdry/mit_urtc_ache_pipeline.git
cd mit_urtc_ache_pipeline
pip install -r requirements.txt

mkdir -p data/raw/admet_raw data/processed figures stability-io/input stability-io/output
# extract admet_raw.zip → data/raw/admet_raw/admet_raw.csv

python scripts/01_chembl_retrieval.py
python scripts/02_admet_filter.py
python scripts/03_umap_clustering.py
python scripts/04_toxicity_scoring.py
python scripts/05_generate_orca_inputs.py
# run ORCA on HPC, place resulting .out files into stability-io/output/
python scripts/06_parse_dft_gaps.py
python scripts/07_compute_final_ranking.py

---

## Data Layout

data/
├─ raw/admet_raw/admet_raw.csv
├─ processed/
│   ├─ admet_filtered.csv
│   ├─ chembl_cleaned.csv
│   ├─ cluster_medoids.csv
│   ├─ dft_gaps.csv
│   ├─ filter_counts.csv
│   ├─ final_ranking.csv
│   └─ surviving_medoids.csv
figures/
├─ umap_clusters.png
└─ dft_gaps.png
stability-io/
├─ input/   # ORCA .inp files
└─ output/  # ORCA .out files
scripts/
├─ 01_chembl_retrieval.py
├─ 02_admet_filter.py
├─ 03_umap_clustering.py
├─ 04_toxicity_scoring.py
├─ 05_generate_orca_inputs.py
├─ 06_parse_dft_gaps.py
└─ 07_compute_final_ranking.py
