# URTC-ACHE-PIPELINE

[![Python](https://img.shields.io/badge/python-3.8%2B-blue)](https://www.python.org/)  
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

An end‑to‑end **ADMET → ML → DFT** pipeline for discovering and ranking novel acetylcholinesterase (AChE) inhibitor scaffolds.

---

## Quick Start

```bash
git clone https://github.com/abdullahchdry/urtc-ache-pipeline.git
cd urtc-ache-pipeline
pip install -r requirements.txt

# All folders and raw data (including admet_raw.zip in data/raw/admet_raw/) 
# are included in the repo—no manual setup required.

python scripts/01_chembl_retrieval.py
python scripts/02_admet_filter.py
python scripts/03_umap_clustering.py
python scripts/04_toxicity_scoring.py
python scripts/05_generate_orca_inputs.py
# submit ORCA jobs, place resulting .out files into stability-io/output/
python scripts/06_parse_dft_gaps.py
python scripts/07_compute_final_ranking.py
```
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
