#!/usr/bin/env python3
"""
06_parse_dft_gaps.py

Description:
    Parse ORCA output files for each medoid scaffold, extract the final
    "ORBITAL ENERGIES" block, determine HOMO/LUMO energies (eV), compute
    raw HOMO–LUMO gaps and a normalized DFT score, then save results.

Usage:
    python scripts/06_parse_dft_gaps.py

Requirements:
    - pandas
    - numpy
    - matplotlib

Data structure:
    stability-io/output/   # contains compound_medoid_<ID>.out files
    data/processed/        # will contain dft_gaps.csv
    figures/               # will contain dft_gaps.png

Before running:
    Copy your ORCA .out files (e.g., from HPC) into stability-io/output/.
"""

import logging
import sys
from pathlib import Path
import re

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Regex for data lines: idx, occupancy, skip, energy (eV)
DATA_LINE = re.compile(r"^\s*(\d+)\s+([\d\.]+)\s+[-\d\.E]+?\s+(-?[\d\.]+)")

# Logging setup
logging.basicConfig(level=logging.INFO, format="%(message)s")
log = logging.getLogger(__name__)

# Define project directories
BASE_DIR = Path(__file__).resolve().parent.parent
STABILITY_OUTPUT = BASE_DIR / "stability-io" / "output"
PROCESSED_DIR = BASE_DIR / "data" / "processed"
FIGURES_DIR = BASE_DIR / "figures"
# Ensure output directories exist
PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)


def parse_last_orb_block(lines):
    """Extract the last 'ORBITAL ENERGIES' block and parse orbitals."""
    headers = [i for i, l in enumerate(lines) if 'ORBITAL ENERGIES' in l.upper()]
    if not headers:
        raise ValueError("No 'ORBITAL ENERGIES' header found")
    start = next((j for j in range(headers[-1]+1, len(lines)) if DATA_LINE.match(lines[j])), None)
    if start is None:
        raise ValueError("No orbital data after header")

    orbitals = []
    for line in lines[start:]:
        m = DATA_LINE.match(line)
        if not m:
            break
        idx, occ, e_ev = int(m.group(1)), float(m.group(2)), float(m.group(3))
        orbitals.append((idx, occ, e_ev))
    if not orbitals:
        raise ValueError("Parsed zero orbitals")
    return orbitals


def find_homo_lumo(orbitals):
    """Return (HOMO, LUMO) energies from orbital list."""
    orbitals.sort(key=lambda x: x[0])
    homo = None
    for _, occ, e in orbitals:
        if occ > 0:
            homo = e
        elif homo is not None:
            return homo, e
    raise ValueError("Could not determine HOMO/LUMO")


def main():
    # Find .out files
    out_files = sorted(STABILITY_OUTPUT.glob('compound_medoid_*.out'))
    if not out_files:
        log.error("No ORCA output files found in %s", STABILITY_OUTPUT)
        sys.exit(1)

    records = []
    for path in out_files:
        match = re.match(r'compound_medoid_(\d+)\.out$', path.name)
        if not match:
            continue
        mid = int(match.group(1))
        lines = path.read_text().splitlines()
        try:
            orbitals = parse_last_orb_block(lines)
            homo, lumo = find_homo_lumo(orbitals)
            gap = lumo - homo
            records.append({'MedoidID': mid, 'HOMO_eV': homo, 'LUMO_eV': lumo, 'Gap_eV': gap})
            log.info("Parsed HOMO/LUMO for Medoid %d: %.4f/%.4f eV", mid, homo, lumo)
        except Exception as e:
            log.error("%s: %s", path.name, e)

    if not records:
        log.error("No orbital data parsed, exiting.")
        sys.exit(1)

    # Build DataFrame and compute normalized score
    df = pd.DataFrame(records).sort_values('MedoidID')
    min_gap, max_gap = df['Gap_eV'].min(), df['Gap_eV'].max()
    df['sDFT'] = (df['Gap_eV'] - min_gap) / (max_gap - min_gap)
    df = df.round({'HOMO_eV':4, 'LUMO_eV':4, 'Gap_eV':4, 'sDFT':4})

    # Save CSV
    csv_path = PROCESSED_DIR / 'dft_gaps.csv'
    df.to_csv(csv_path, index=False)
    log.info("Saved DFT gaps to %s", csv_path)

    # Plot and save
    img_path = FIGURES_DIR / 'dft_gaps.png'
    plt.figure(figsize=(8,4))
    plt.bar(df['MedoidID'].astype(str), df['Gap_eV'])
    plt.xlabel('Medoid ID')
    plt.ylabel('HOMO–LUMO Gap (eV)')
    plt.title('HOMO–LUMO Gaps of Medoid Scaffolds')
    plt.tight_layout()
    plt.savefig(img_path, dpi=300)
    plt.close()
    log.info("Saved HOMO–LUMO gap chart to %s", img_path)

if __name__ == '__main__':
    main()
