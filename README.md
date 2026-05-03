# SMILES Desalter

This script is designed to remove salts from SMILES datasets coming from any chemical databases. It automatically removes the smaller fragments from the dataset, ensuring you are left with only the primary molecule of interest.

For multi-fragment SMILES (separated by a `.`), this script defaults to keeping only the largest fragment (determined by heavy atom count) and discarding the smaller salt/solvent fragments. Alternatively, you can choose to drop the entire row if a salt is detected.

## Features
* **Auto-detection**: Automatically finds the SMILES column by looking for common headers (e.g., `smiles`, `smi`, `canonical_smiles`).
* **Flexible Output**: Automatically generates a new `_desalted.csv` file or lets you specify an output path.
* **Modes**: 
  * *Default*: Replaces the salt SMILES with the largest single fragment.
  * *Drop*: Removes the entire entry if it contains a salt.
* **Customizable**: Supports custom column names and custom CSV delimiters.

## Requirements

This script requires [RDKit](https://www.rdkit.org/). 

You can install the required dependencies using pip:

```bash
pip install -r requirements.txt
