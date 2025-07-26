#!/usr/bin/env python3
"""
05_generate_orca_inputs.py

Description:
    Generate ORCA input (.inp) files for each medoid scaffold listed in
    `surviving_medoids.csv` located in `data/processed/`.
    Creates one subfolder per scaffold under `stability-io/input/compound_medoid_<MedoidID>/` containing:
      - compound_medoid_<MedoidID>.inp  : B3LYP/6-31G(d) Opt using 16 cores

    Multiplicity is set to 1 for even-electron systems, 2 for odd.
    Charge is assumed zero.

Usage:
    python scripts/05_generate_orca_inputs.py

Requirements:
    - pandas
    - rdkit

Data structure:
    data/processed/surviving_medoids.csv   # input CSV with MedoidID & smiles
    stability-io/
      input/                              # will receive compound_medoid_<ID>/ folders

Before running:
    Ensure `data/processed/surviving_medoids.csv` exists from toxicity step.
"""

import logging
import sys
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

# ORCA input template
ORCA_INP = '''! B3LYP 6-31G(d) Opt TightSCF
%pal
   nprocs 16
end

* xyz 0 {mult}
{xyz}
*'''

# Logging setup
logging.basicConfig(level=logging.INFO, format="%(message)s")
log = logging.getLogger(__name__)

# Define project directories
BASE_DIR = Path(__file__).resolve().parent.parent
PROCESSED_DIR = BASE_DIR / "data" / "processed"
STABILITY_IO_DIR = BASE_DIR / "stability-io"
INPUT_DIR = STABILITY_IO_DIR / "input"
# Ensure input directory exists
INPUT_DIR.mkdir(parents=True, exist_ok=True)


def embed_molecule(smiles: str) -> Chem.Mol:
    """Embed SMILES to 3D and run UFF optimization."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    if AllChem.EmbedMolecule(mol, params) != 0:
        raise RuntimeError(f"Embedding failed for SMILES: {smiles}")
    AllChem.UFFOptimizeMolecule(mol)
    return mol


def mol_to_xyz(mol: Chem.Mol) -> str:
    """Convert molecule to XYZ coordinate block."""
    conf = mol.GetConformer()
    lines = []
    for atom in mol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        lines.append(f"{atom.GetSymbol():<2} {pos.x:10.6f} {pos.y:10.6f} {pos.z:10.6f}")
    return "\n".join(lines)


def determine_multiplicity(mol: Chem.Mol) -> int:
    """Multiplicity: 1 if even-electron, else 2 (charge=0)."""
    total_protons = sum(a.GetAtomicNum() for a in mol.GetAtoms())
    return 1 if (total_protons % 2 == 0) else 2


def main():
    # Input CSV path
    input_csv = PROCESSED_DIR / "surviving_medoids.csv"

    # Validate input CSV
    if not input_csv.exists():
        log.error(f"CSV not found at {input_csv}")
        sys.exit(1)
    df = pd.read_csv(input_csv, dtype=str)
    log.info(f"Loaded {len(df)} medoids from {input_csv}")

    # Identify required columns
    cols_lower = {c.lower(): c for c in df.columns}
    if 'medoidid' not in cols_lower or 'smiles' not in cols_lower:
        log.error("Input CSV must contain 'MedoidID' and 'smiles' columns.")
        sys.exit(1)
    id_col = cols_lower['medoidid']
    smi_col = cols_lower['smiles']

    # Generate ORCA input for each medoid
    for _, row in df.iterrows():
        mid = row[id_col]
        smiles = row[smi_col]
        log.info(f"Processing MedoidID={mid}")

        # Embed and optimize geometry
        mol3d = embed_molecule(smiles)
        mult = determine_multiplicity(mol3d)
        log.info(f"  Multiplicity set to {mult}")

        # Prepare output folder under stability-io/input
        medoid_folder = INPUT_DIR / f"compound_medoid_{mid}"
        medoid_folder.mkdir(exist_ok=True)

        # Write ORCA input file
        xyz_block = mol_to_xyz(mol3d)
        inp_content = ORCA_INP.format(mult=mult, xyz=xyz_block)
        inp_path = medoid_folder / f"compound_medoid_{mid}.inp"
        inp_path.write_text(inp_content)
        log.info(f"  Wrote ORCA input to {inp_path}")

if __name__ == "__main__":
    main()
