#!/usr/bin/env python3
"""
02_admet_filter.py

Description:
    Apply strict ADMET and Lipinski Rule of 5 filters (with a TPSA ≥ 20 Å² floor)
    and loose pharmacokinetic property filters to the ADMETlab 3.0 dataset.
    Outputs a filtered CSV plus a summary table of counts at each stage into data/processed.

Usage:
    python scripts/02_admet_filter.py

Requirements:
    - pandas

Data structure:
    data/raw/admet_raw/      # extract admet_raw.zip here to produce admet_raw.csv
    data/processed/          # script will save filter_counts.csv and admet_filtered.csv here 

Before running:
    Ensure you have extracted admet_raw.zip into data/raw/admet_raw/ to obtain admet_raw.csv.
"""

import logging
import sys
from pathlib import Path

import pandas as pd

# Define project directories
BASE_DIR = Path(__file__).resolve().parent.parent
RAW_DIR = BASE_DIR / "data" / "raw" / "admet_raw"
PROCESSED_DIR = BASE_DIR / "data" / "processed"
PROCESSED_DIR.mkdir(parents=True, exist_ok=True)

# File paths
INPUT_FILE = RAW_DIR / "admet_raw.csv"
OUTPUT_FILE = PROCESSED_DIR / "admet_filtered.csv"
TABLE_FILE = PROCESSED_DIR / "filter_counts.csv"

# Expected initial count (adjust if needed)
INITIAL_COUNT = 1393

# Logging setup
logging.basicConfig(level=logging.INFO, format="%(message)s")
log = logging.getLogger(__name__)


def apply_strict_filters(df: pd.DataFrame) -> pd.DataFrame:
    """Filter by strict ADMET + Rule of 5 criteria."""
    strict = (
        (df['logS']   >= -4.5) &
        (df['logP']   <=  5.0) &
        (df['TPSA']   >= 20.0) &
        (df['TPSA']   <=100.0) &
        (df['nRot']   <= 10  ) &
        (df['BBB']    >=  0.70) &
        (df['MW']     <=500.0) &
        (df['nHD']    <=  5  ) &
        (df['nHA']    <= 10  )
    )
    return df[strict].reset_index(drop=True)


def apply_loose_filters(df: pd.DataFrame) -> pd.DataFrame:
    """Filter by loose pharmacokinetic criteria."""
    loose = (
        (df['logVDss']   <=   5) &
        (df['cl-plasma'] <=  50) &
        (df['t0.5']      >= 0.3) &
        (df['PPB']       <= 99.5) &
        (df['Fsp3']      >=  0.1)
    )
    return df[loose].reset_index(drop=True)


def main():
    # Check input exists
    if not INPUT_FILE.exists():
        log.error(f"Error: input not found at {INPUT_FILE}")
        sys.exit(1)

    # Load data
    df_raw = pd.read_csv(INPUT_FILE)
    log.info(f"Initial molecules: {len(df_raw)} (expected ~{INITIAL_COUNT})")

    # Apply strict filters
    df_strict = apply_strict_filters(df_raw)
    log.info(f"After strict ADMET filtering: {len(df_strict)}")

    # Apply loose filters
    df_final = apply_loose_filters(df_strict)
    log.info(f"After loose PK filtering: {len(df_final)}")

    # Save filtered data
    df_final.to_csv(OUTPUT_FILE, index=False)
    log.info(f"Filtered dataset saved to {OUTPUT_FILE}")

    # Summary counts
    summary = pd.DataFrame({
        "Stage": [
            "Initial ADMET hits",
            "After strict ADMET filter",
            "After loose PK filter"
        ],
        "Count": [
            len(df_raw),
            len(df_strict),
            len(df_final)
        ]
    })
    summary.to_csv(TABLE_FILE, index=False)
    log.info(f"Filter summary saved to {TABLE_FILE}")

if __name__ == "__main__":
    main()
