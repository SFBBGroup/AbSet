
# ChemFeatX
![Python 3.13](https://img.shields.io/static/v1?label=python&message=3.13&color=blue&style=flat-square)

This repository contains a script to filter structures based on their resolution, using resolution data available from the Protein Data Bank (PDB)

## Setup

The filter_resolution.py script can be executed directly from the terminal. Alternatively, a Google Colab notebook (filter_resolution.ipynb) is provided for users who prefer a browser-based environment.


## Prerequisites

Biopython

## Running the Script

The resolution threshold can be modified.

Place one or more `.pdb` files in the script's execution directory. To execute, run:

```bash
python resolution_filter.py
```

### Output format

At the end of the execution, wo new folders are created: one containing the PDB files with a resolution up to the desired threshold, and another with those exceeding the threshold:
```bash
./pdb_up_to_threshold/pdb_file.pdb
./pdb_above_threshold/pdb_file.pdb
```
