#!/usr/bin/env python3
"""
07_compute_final_ranking.py

Description:
    Compute composite rankings for medoid scaffolds by combining ADMET-based
    distance scores (d_i) with DFT stability scores (sDFT) using:

      S_i = α·(1 − d_i) + (1 − α)·sDFT_i

Usage:
    python scripts/07_compute_final_ranking.py \
        --admet data/processed/admet_filtered.csv \
        --medoids data/processed/surviving_medoids.csv \
        --dft data/processed/dft_gaps.csv \
        --alpha 0.6 \
        --out data/processed/final_ranking.csv

Requirements:
    - pandas
    - numpy
    - scikit-learn

Data structure:
    data/processed/   # contains admet_filtered.csv, surviving_medoids.csv, dft_gaps.csv;
                      # writes final_ranking.csv here
"""

import argparse
from pathlib import Path

import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler

# Define project directories
BASE_DIR = Path(__file__).resolve().parent.parent
PROCESSED_DIR = BASE_DIR / "data" / "processed"
PROCESSED_DIR.mkdir(parents=True, exist_ok=True)

# 13 ADMET features for distance calculation
FEATURES = [
    'logS','logP','TPSA','nRot','BBB','MW','nHD','nHA',
    'logVDss','cl-plasma','t0.5','PPB','Fsp3'
]


def compute_admet_distance(df: pd.DataFrame) -> pd.DataFrame:
    """Compute normalized Euclidean ADMET distance d_i."""
    X = StandardScaler().fit_transform(df[FEATURES])
    distances = np.linalg.norm(X, axis=1)
    df['d_i'] = distances / distances.max()
    return df


def main():
    parser = argparse.ArgumentParser(
        description="Compute ADMET & DFT composite ranking for medoid scaffolds"
    )
    parser.add_argument('--admet', default=str(PROCESSED_DIR / 'admet_filtered.csv'),
                        help='Path to ADMET descriptors CSV')
    parser.add_argument('--medoids', default=str(PROCESSED_DIR / 'surviving_medoids.csv'),
                        help='Path to surviving medoids CSV')
    parser.add_argument('--dft', default=str(PROCESSED_DIR / 'dft_gaps.csv'),
                        help='Path to DFT gaps with sDFT CSV')
    parser.add_argument('--alpha', type=float, default=0.6,
                        help='Weighting factor α')
    parser.add_argument('--out', default=str(PROCESSED_DIR / 'final_ranking.csv'),
                        help='Output CSV filename')
    args = parser.parse_args()

    # Load and process ADMET data
    admet_df = pd.read_csv(args.admet)
    admet_df = compute_admet_distance(admet_df)

    # Load surviving medoids
    meds_df = pd.read_csv(args.medoids, dtype={'MedoidID': int, 'smiles': str})

    # Merge ADMET distances into medoids
    merged = meds_df.merge(admet_df[['smiles', 'd_i']], on='smiles', how='left')

    # Load DFT data
    dft_df = pd.read_csv(args.dft, dtype={'MedoidID': int})

    # Merge DFT scores into merged
    merged = merged.merge(
        dft_df[['MedoidID', 'sDFT', 'Gap_eV', 'HOMO_eV', 'LUMO_eV']],
        on='MedoidID', how='left'
    )

    # Compute composite score S_i
    alpha = args.alpha
    merged['S_i'] = alpha * (1 - merged['d_i']) + (1 - alpha) * merged['sDFT']

    # Sort and save
    cols = ['MedoidID', 'smiles', 'd_i', 'sDFT', 'S_i', 'Gap_eV', 'HOMO_eV', 'LUMO_eV']
    final_df = merged.sort_values('S_i', ascending=False)
    final_df.to_csv(args.out, columns=cols, index=False)
    print(f"Wrote {len(final_df)} entries to {args.out}")

    # Print top 3 medoid scaffolds
    print("\nTop 3 medoid scaffolds by composite score:")
    print(final_df[['MedoidID', 'smiles', 'S_i']].head(3).to_string(index=False))

if __name__ == '__main__':
    main()
