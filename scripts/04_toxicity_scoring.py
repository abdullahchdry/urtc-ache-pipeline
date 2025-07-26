#!/usr/bin/env python3
"""
04_toxicity_scoring.py

Description:
    Compute composite toxicity scores for medoid scaffolds using weighted
    probabilities from four in silico endpoints: hERG, DILI, Ames, and SkinSens.
    Outputs a full medoid list with ToxScore (ascending by risk) and a survivors list
    (lower half) saved to data/processed/surviving_medoids.csv.

Usage:
    python scripts/04_toxicity_scoring.py

Requirements:
    - pandas
    - numpy

Data structure:
    data/processed/       # contains cluster_medoids.csv and will store toxicity outputs

Before running:
    Ensure clustering step produced data/processed/cluster_medoids.csv.
"""

import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# Define project directories
BASE_DIR = Path(__file__).resolve().parent.parent
PROCESSED_DIR = BASE_DIR / "data" / "processed"
PROCESSED_DIR.mkdir(parents=True, exist_ok=True)

# File patterns and paths
INPUT_PATTERN = "cluster_medoids*.csv"
FULL_OUTPUT = PROCESSED_DIR / "cluster_medoids_toxscore.csv"
SURVIVORS_OUTPUT = PROCESSED_DIR / "surviving_medoids.csv"

# Toxicity endpoint weights
WEIGHTS = {
    'herg': 0.35,
    'dili': 0.35,
    'ames': 0.20,
    'skin': 0.10
}

# Logging setup
logging.basicConfig(level=logging.INFO, format="%(message)s")
log = logging.getLogger(__name__)

def find_input_file() -> Path:
    """Locate the medoid CSV file in processed directory."""
    files = list(PROCESSED_DIR.glob(INPUT_PATTERN))
    if not files:
        log.error(f"No file matching '{INPUT_PATTERN}' in {PROCESSED_DIR}")
        sys.exit(1)
    return files[0]


def map_toxicity_columns(df: pd.DataFrame) -> dict:
    """Map endpoint keys to DataFrame column names."""
    mapping = {}
    for key in WEIGHTS:
        for col in df.columns:
            if key in col.lower():
                mapping[key] = col
                break
    missing = set(WEIGHTS) - set(mapping)
    if missing:
        log.error(f"Missing toxicity columns: {missing}")
        sys.exit(1)
    log.info(f"Mapped toxicity columns: {mapping}")
    return mapping


def compute_toxscore(df: pd.DataFrame, col_map: dict) -> pd.DataFrame:
    """Compute weighted composite ToxScore."""
    for col in col_map.values():
        df[col] = pd.to_numeric(df[col], errors='coerce')
    df['ToxScore'] = sum(df[col_map[k]] * w for k, w in WEIGHTS.items())
    return df


def main():
    log.info("Starting toxicity scoring pipeline...")
    input_file = find_input_file()
    log.info(f"Loading medoids from {input_file}")
    df = pd.read_csv(input_file)
    log.info(f"Loaded {len(df)} medoid records")

    # Identify identifier column
    id_col = next((c for c in df.columns if 'smiles' in c.lower()), df.columns[0])
    log.info(f"Using '{id_col}' as identifier column")

    # Map toxicity columns and compute scores
    col_map = map_toxicity_columns(df)
    df = compute_toxscore(df, col_map)
    df = df.sort_values('ToxScore').reset_index(drop=True)
    log.info("Composite ToxScore computed and medoids ranked")

    # Save full ranked list
    df.to_csv(FULL_OUTPUT, index=False)
    log.info(f"Full medoid list with ToxScore saved to {FULL_OUTPUT}")

    # Select and save survivors (lower half)
    cutoff = df['ToxScore'].median()
    survivors = df[df['ToxScore'] <= cutoff].copy()
    survivors.to_csv(SURVIVORS_OUTPUT, index=False)
    log.info(f"Survivors (ToxScore <= {cutoff:.4f}) saved to {SURVIVORS_OUTPUT}")

if __name__ == "__main__":
    main()
