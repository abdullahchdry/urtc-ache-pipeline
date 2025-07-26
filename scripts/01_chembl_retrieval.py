#!/usr/bin/env python3
"""
01_chembl_retrieval.py

Description:
    Fetch human acetylcholinesterase (AChE) inhibitors from ChEMBL (TARGET_ID CHEMBL220)
    with IC50 <= IC50_THRESH nM. Override missing SMILES for CHEMBL2448138, strip salts,
    dedupe by minimum IC50 and unique SMILES, then save results to data/processed/chembl_cleaned.csv.

Usage:
    python scripts/01_chembl_retrieval.py

Requirements:
    - chembl_webresource_client
    - RDKit

Data structure:
    data/processed/       # will contain chembl_cleaned.csv
"""

import logging
from pathlib import Path

import pandas as pd
from chembl_webresource_client.new_client import new_client
from rdkit import Chem

# Define project directories
BASE_DIR = Path(__file__).resolve().parent.parent
RAW_DIR = BASE_DIR / "data" / "raw"
PROCESSED_DIR = BASE_DIR / "data" / "processed"
PROCESSED_DIR.mkdir(parents=True, exist_ok=True)

# Configuration
TARGET_ID = "CHEMBL220"        # AChE target
IC50_THRESH = 100              # threshold in nM
OUTPUT_FILE = PROCESSED_DIR / "chembl_cleaned.csv"

# Override entry for missing SMILES
OVERRIDE_ID = "CHEMBL2448138"
OVERRIDE_SMILES = (
    "COc1cc2c(cc1OC)C(=O)C(CC1CCN(CCNc3c4c(nc5ccccc35)CCCC4)CC1)C2"
)

# Logging setup
logging.basicConfig(level=logging.INFO, format="%(message)s")
log = logging.getLogger(__name__)

def strip_salt(smiles: str) -> str:
    """Return canonical SMILES of the largest fragment by heavy-atom count."""
    if not smiles:
        return ""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return ""
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    parent = max(frags, key=lambda m: m.GetNumHeavyAtoms())
    return Chem.MolToSmiles(parent, canonical=True)

def fetch_activities() -> pd.DataFrame:
    """Retrieve IC50 activity records for the target and apply SMILES override."""
    client = new_client.activity
    query = {
        "target_chembl_id": TARGET_ID,
        "standard_type": "IC50",
        "standard_units": "nM",
        "standard_value__lte": IC50_THRESH,
        "pchembl_value__isnull": False,
    }
    activities = client.filter(**query)

    records = []
    for act in activities:
        chembl_id = act.get("molecule_chembl_id")
        raw_smiles = act.get("canonical_smiles", "")
        smiles = OVERRIDE_SMILES if chembl_id == OVERRIDE_ID else strip_salt(raw_smiles)

        try:
            ic50 = float(act["standard_value"])
            pchembl = float(act["pchembl_value"])
        except (TypeError, ValueError):
            continue

        records.append({
            "CHEMBL_ID": chembl_id,
            "SMILES": smiles,
            "IC50_nM": ic50,
            "Relation": act.get("standard_relation", ""),
            "Assay_ID": act.get("assay_chembl_id", ""),
            "Reference_ID": act.get("document_chembl_id", ""),
            "pChEMBL": pchembl,
        })
    return pd.DataFrame(records)

def dedupe_ic50(df: pd.DataFrame) -> pd.DataFrame:
    """Keep the record with lowest IC50 per CHEMBL_ID."""
    return df.sort_values("IC50_nM").groupby("CHEMBL_ID", as_index=False).first()

def dedupe_smiles(df: pd.DataFrame) -> pd.DataFrame:
    """Drop duplicate SMILES, retaining first occurrence."""
    return df.drop_duplicates(subset="SMILES", keep="first").reset_index(drop=True)

def main():
    log.info("Starting retrieval for %s (IC50 <= %d nM)", TARGET_ID, IC50_THRESH)
    df = fetch_activities()
    log.info("Fetched %d records", len(df))

    df = dedupe_ic50(df)
    log.info("After IC50 dedupe: %d unique CHEMBL IDs", len(df))

    df = dedupe_smiles(df)
    log.info("After SMILES dedupe: %d unique molecules", len(df))

    df.to_csv(OUTPUT_FILE, index=False)
    log.info("Results saved to %s", OUTPUT_FILE)
    log.info("Note: SMILES override applied for %s.", OVERRIDE_ID)

if __name__ == "__main__":
    main()
