#!/usr/bin/env python3
"""
Desalt a CSV file containing SMILES strings.

For multi-fragment SMILES (salts), keeps only the largest fragment
by heavy atom count. Automatically detects the SMILES column.

Usage:
    python3 desalt.py input.csv                     # writes to input_desalted.csv
    python3 desalt.py input.csv -o output.csv       # writes to output.csv
    python3 desalt.py input.csv --smiles-col smi    # specify SMILES column name
    python3 desalt.py input.csv --drop-salts        # drop salt entries instead of desalting
"""

######################## Libraries to import ########################
import argparse
import csv
import os
import sys

from rdkit import Chem, RDLogger

######################## Suppress RDKit warnings ########################
RDLogger.DisableLog("rdApp.*")

COMMON_SMILES_HEADERS = {"smiles", "smi", "canonical_smiles", "cansmiles", "isosmiles", "structure"}


def find_smiles_column(headers):
    """Auto-detect the SMILES column from headers."""
    for h in headers:
        if h.lower() in COMMON_SMILES_HEADERS:
            return h
    # Fallback: check if any column contains 'smi' as substring
    for h in headers:
        if "smi" in h.lower():
            return h
    return None


def desalt_smiles(smi):
    """
    Return the largest fragment (by heavy atom count) from a SMILES string.
    Returns None if no valid fragment is found.
    """
    frags = smi.split(".")
    best_smi = None
    best_size = -1
    for frag in frags:
        mol = Chem.MolFromSmiles(frag)
        if mol is not None:
            n = mol.GetNumHeavyAtoms()
            if n > best_size:
                best_size = n
                best_smi = frag
    return best_smi


def main():
    parser = argparse.ArgumentParser(description="Remove salt fragments from SMILES in a CSV file.")
    parser.add_argument("input", help="Input CSV file")
    parser.add_argument("-o", "--output", help="Output CSV file (default: <input>_desalted.csv)")
    parser.add_argument("--smiles-col", help="Name of the SMILES column (auto-detected if not given)")
    parser.add_argument("--drop-salts", action="store_true",
                        help="Drop salt entries entirely instead of keeping the largest fragment")
    parser.add_argument("--delimiter", default=",", help="CSV delimiter (default: comma)")
    args = parser.parse_args()

    if not os.path.isfile(args.input):
        print(f"Error: {args.input} not found", file=sys.stderr)
        sys.exit(1)

    if args.output is None:
        base, ext = os.path.splitext(args.input)
        args.output = f"{base}_desalted{ext}"

    with open(args.input, "r") as fin:
        reader = csv.DictReader(fin, delimiter=args.delimiter)
        headers = reader.fieldnames

        # Detect SMILES column
        smiles_col = args.smiles_col or find_smiles_column(headers)
        if smiles_col is None or smiles_col not in headers:
            print(f"Error: Could not detect SMILES column. Headers: {headers}", file=sys.stderr)
            print("Use --smiles-col to specify it.", file=sys.stderr)
            sys.exit(1)

        print(f"Input:        {args.input}")
        print(f"Output:       {args.output}")
        print(f"SMILES col:   {smiles_col}")
        print(f"Mode:         {'drop salts' if args.drop_salts else 'keep largest fragment'}")
        print()

        total = 0
        salts = 0
        dropped = 0
        written = 0

        with open(args.output, "w", newline="") as fout:
            writer = csv.DictWriter(fout, fieldnames=headers, delimiter=args.delimiter)
            writer.writeheader()

            for row in reader:
                total += 1
                smi = row[smiles_col].strip()

                if "." in smi:
                    salts += 1
                    if args.drop_salts:
                        dropped += 1
                        continue
                    desalted = desalt_smiles(smi)
                    if desalted is None:
                        dropped += 1
                        continue
                    row[smiles_col] = desalted
                else:
                    mol = Chem.MolFromSmiles(smi)
                    if mol is None:
                        dropped += 1
                        continue

                writer.writerow(row)
                written += 1

    print(f"Total entries:      {total}")
    print(f"Salt forms found:   {salts}")
    print(f"Dropped:            {dropped}")
    print(f"Written:            {written}")


if __name__ == "__main__":
    main()
