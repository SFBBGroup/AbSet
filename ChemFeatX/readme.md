# ChemFeatX
![Python 3.13](https://img.shields.io/static/v1?label=python&message=3.13&color=blue&style=flat-square)

This repository contains a script for calculating molecular descriptors to capture the characteristics of amino acid residues. The calculated descriptors include:
* Solvent Accessible Surface Area 
* Relative Accessible Surface Area
* Atomic depth
* Potrusion index
* Hydrophobicity
* Sequence
* Half-sphere exposure calculations
* Cα coordinates
* ϕ and ψ dihedral angles
* Secondary structure of the protein

Additionally, the script supports multiprocessing to efficiently compute multiple `.pdb` files.

## Setup

We provide a way to run the program using a Conda environment:

### Conda environment
```bash
$ conda env create -f environment.yml
$ conda activate ChemFeatX
```

## Prerequisites

Before running, you need to install the following tools on your machine:
* [DSSP](https://github.com/PDB-REDO/dssp)
* [MGLtools](https://ccsb.scripps.edu/mgltools/)

Make sure DSSP is functional by running a test command:
```bash
mkdssp {pdb_file} --output-format dssp > {filename}.tbl
```

Additionally, make sure that the PDB files do not contain heteroatoms.

## Running the Script

Place one or more `.pdb` files in the script's execution directory. To execute, run:

```bash
python molecular_descriptors_calculator.py
```

### Output format

At the end of the execution, a .csv file will be generated at:
```bash
./output/pdb_name.csv
```
