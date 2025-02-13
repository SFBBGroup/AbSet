import numpy as np
import pandas as pd
from Bio.PDB import PDBParser, PDBIO, NeighborSearch, Selection, DSSP
from Bio.PDB.SASA import ShrakeRupley
from Bio.PDB.ResidueDepth import get_surface, min_dist as bio_min_dist
from collections import defaultdict
from glob import glob
from Bio import PDB
import multiprocessing
from Bio.PDB.ResidueDepth import ResidueDepth
import subprocess
import os

OUT_PATH = './output'

# Function to calculate Solvent Accessible Surface Area (SASA) and Relative Accessible Surface Area (RASA)
def calculate_SASA_RASA(model):
    # Get all chain IDs in the model
    chains = [chain.get_id() for chain in model.get_chains()]
    sr = ShrakeRupley()          # Initialize SASA calculator
    sr.compute(model, level="R") # Compute SASA at the residue level
    data_ASA = []

    for chain_id in chains:
        chain = model[chain_id]
        for residue in chain:
            # Get residue number and store SASA data
            residue_num = str(residue.get_id()[1]) + residue.get_id()[2]
            residue_num = residue_num.strip()
            data_ASA.append((chain_id, residue_num, residue.get_resname(), round(residue.sasa, 2)))

    # Create DataFrame with SASA data
    df_SASA = pd.DataFrame(data_ASA, columns=['chain', 'residue_num', 'residue_name', 'SASA'])

    # Define maximum ASA values for amino acids
    asa_max = {"ALA": 107.24, "ARG": 233.01, "ASN": 150.85, "ASP": 144.06,
               "CYS": 131.46, "GLN": 177.99, "GLU": 171.53, "GLY": 80.54,
               "HIS": 180.93, "ILE": 173.40, "LEU": 177.87, "LYS": 196.14,
               "MET": 186.80, "PHE": 200.93, "PRO": 133.78, "SER": 115.30,
               "THR": 136.59, "TRP": 240.12, "TYR": 213.21, "VAL": 149.34}

    # Calculate RASA and add it to the DataFrame
    df_SASA_RASA = df_SASA.copy()
    df_SASA_RASA['RASA'] = round((df_SASA_RASA['SASA'] / df_SASA_RASA['residue_name'].str[:3].map(asa_max).fillna(1)) * 100, 3)
    
    return df_SASA_RASA

# Function to calculate Half-Sphere Exposure (HSE) values
def calculate_HSE(model):
    chains = [chain.get_id() for chain in model.get_chains()]
    data_hse = []
    RADIUS = 12.0 # Define radius for HSE calculation
    hse_CA = PDB.HSExposure.HSExposureCA(model, RADIUS)
    hse_CB = PDB.HSExposure.HSExposureCB(model, RADIUS)


    for chain_id in chains:
        chain = model[chain_id]
        for residue in chain:
            residue_num = str(residue.get_id()[1]) + residue.get_id()[2]
            residue_num = residue_num.strip()
            xse = residue.xtra # Retrieve HSE data from residue
            entry = {
                'chain': chain_id,
                'residue_num': residue_num ,
                'EXP_HSE_B_U': xse.get('EXP_HSE_B_U', None),
                'EXP_HSE_B_D': xse.get('EXP_HSE_B_D', None),
                'EXP_HSE_A_U': xse.get('EXP_HSE_A_U', None),
                'EXP_HSE_A_D': xse.get('EXP_HSE_A_D', None)
            }
            data_hse.append(entry)
    
    # Convert HSE data to a DataFrame
    return pd.DataFrame(data_hse)

# Function to calculate hydrophobicity based on a predefined scale
def calculate_HDR(model):
    chains = [chain.get_id() for chain in model.get_chains()]
    data_HDR = []

    hdr_scale = {'ILE': 4.5, 'VAL': 4.2, 'LEU': 3.8, 'PHE': 2.8, 'CYS': 2.5,
                 'MET': 1.9, 'ALA': 1.8, 'GLY': -0.4, 'THR': -0.7, 'SER': -0.8,
                 'TRP': -0.9, 'TYR': -1.3, 'PRO': -1.6, 'HIS': -3.2, 'GLN': -3.5,
                 'ASN': -3.5, 'GLU': -3.5, 'ASP': -3.5, 'LYS': -3.9, 'ARG': -4.0}

    for chain_id in chains:
        chain = model[chain_id]
        for residue in chain:
            residue_num = str(residue.get_id()[1]) + residue.get_id()[2]
            residue_num = residue_num.strip()
            data_HDR.append([chain_id, residue_num, residue.get_resname()])

    df_HDR = pd.DataFrame(data_HDR, columns=['chain', 'residue_num', 'residue_name'])
    df_HDR['hydrophobicity'] = df_HDR['residue_name'].map(hdr_scale)
    df_HDR = df_HDR.drop(columns=['residue_name'])

    return df_HDR

# Function to calculate residue depth (DPX) for each residue in the structure
def calculate_dpx(structure):
    data_dpx = []
    rd = ResidueDepth(structure[0]) # Calculate residue depth for the first model

    for chain in structure[0]:
        for residue in chain:
            residue_num = str(residue.get_id()[1]) + residue.get_id()[2]
            residue_num = residue_num.strip()
            residue_depth = rd[chain.id, residue.id]
            residue_depth = round(residue_depth[0], 3)

            data_dpx.append({
                'chain': chain.get_id(),
                'residue_num': residue_num,
                'mean_dpx': residue_depth,
            })

    return pd.DataFrame(data_dpx)

# Function to calculate the protrusion index (IP) for residues
def calculate_IP(model):
    data_IP = []
    results = defaultdict(list)
    atoms = Selection.unfold_entities(model, 'A') # Get all atoms in the model
    ns = NeighborSearch(atoms) # Initialize neighbor search for proximity analysis

    for atom in atoms:
        close_atoms = ns.search(atom.coord, 10) # Search atoms within 10 Å radius
        parent = atom.get_parent()              # Get residue containing the atom
        residue_name = parent.get_resname()
        residue_df = str(parent.get_id()[1]) + parent.get_id()[2]
        residue_id = residue_df.strip()
        chain_id = parent.get_full_id()[2]
        key = (chain_id, residue_id, residue_name)
        results[key].append(len(close_atoms))

    for key, value in results.items():
        data_residues = (list(key) + value)
        chain = data_residues[0]
        residue_num = data_residues[1]
        residue_name = data_residues[2]
        values = data_residues[3:]
        mean_contacts = sum(values) / len(values)
        max_contacts = max(values)
        min_contacts = min(values)
        mean_value = indice_protrusion(mean_contacts)
        value_max = indice_protrusion(max_contacts)
        value_min = indice_protrusion(min_contacts)
        data_IP.append([chain, residue_num, mean_value, value_max, value_min])

    return pd.DataFrame(data_IP, columns=['chain', 'residue_num', 'mean_IP', 'max_IP', 'min_IP'])

# Function to calculate the Cα (alpha carbon) coordinates for each residue in a PDB model
def calculate_CA(model, pdb_file):
    chains = [chain.get_id() for chain in model.get_chains()]
    data_DSSP = []
    data_CA = []
    parser = PDBParser()

    for chain_id in chains:
        chain = model[chain_id]
        for residue in chain:
            residue_num = str(residue.get_id()[1]) + residue.get_id()[2]
            residue_num = residue_num.strip()
            res_id = residue.get_id()

            # Check if the residue contains a Cα atom
            if residue.has_id("CA"):
                ca_atom = residue["CA"]
                data_CA.append({
                    'chain': chain.id,
                    'residue_num': residue_num,
                    'CA_x': ca_atom.coord[0],
                    'CA_y': ca_atom.coord[1],
                    'CA_z': ca_atom.coord[2]
                })

    df_CA = pd.DataFrame(data_CA)

    # Return the Cα DataFrame
    return df_CA

# Helper function to calculate the protrusion index (CX)
def indice_protrusion(num_contatos):
    radius = 10.0
    mean_atomic_volume = 20.1
    volume_int = num_contatos * mean_atomic_volume
    volume_sphere = (4 / 3) * np.pi * (radius ** 3)
    volume_ext = volume_sphere - volume_int
    cx_value = round(volume_ext / volume_int, 3)

    # Return the protrusion index
    return cx_value

# Function to calculate DSSP values for secondary structure
def calculate_DSSP(pdb_file):
    filename = f'{pdb_file[:-4]}'
    subprocess.run(f"mkdssp {pdb_file} --output-format dssp > {filename}.tbl", shell=True)
    dssp_file = f'{filename}.tbl'

    # Parse the DSSP file and write relevant data to a CSV
    with open(f'{filename}_dssp.csv', "w") as output_file:
        subprocess.run(
            f"grep -A500000000000 'RESIDUE AA' {dssp_file} | "
            f"awk -F '' '{{print $6$7$8$9$10$11\",\"$12\",\"$17\",\"$104$105$106$107$108$109\",\"$110$111$112$113$114$115$116}}' | "
            f"tr -d ' '",
            shell=True,
            stdout=output_file
        )
    
    # Update the header of the CSV file
    with open(f'{filename}_dssp.csv', 'r+') as file:
        lines = file.readlines()
        lines[0] = 'residue_num,chain,secondary_structure,phi,psi\n'
        file.seek(0) 
        file.writelines(lines)
    
    df_DSSP = pd.read_csv(f'{filename}_dssp.csv')
    df_DSSP['secondary_structure'].fillna('L', inplace=True) # Fill missing values with 'L' (loop)

    # Return the DSSP DataFrame
    return df_DSSP

# Function to process a single PDB file
def process_pdb_file(pdb_file):
    # Get the structure pdb
    parser = PDBParser()
    structure = parser.get_structure("estrutura", pdb_file)
    model = structure[0]

    # Calculate various descriptors
    df_SASA_RASA = calculate_SASA_RASA(model)
    df_HSE = calculate_HSE(model)
    df_HDR = calculate_HDR(model)
    df_IP = calculate_IP(model)
    df_dpx = calculate_dpx(structure)
    df_CA = calculate_CA(model, pdb_file)
    df_DSSP = calculate_DSSP(pdb_file)

    # Merge all calculated DataFrames
    df_final = df_SASA_RASA.merge(df_HSE, on=['chain', 'residue_num'], how='left')
    df_final = df_final.merge(df_HDR[['chain', 'residue_num', 'hydrophobicity']], on=['chain', 'residue_num'], how='left')
    df_final = df_final.merge(df_IP, on=['chain', 'residue_num'], how='left')
    df_final = df_final.merge(df_dpx, on=['chain', 'residue_num'], how='left')
    df_final = df_final.merge(df_CA, on=['chain', 'residue_num'], how='left')
    df_final['residue_num'] = df_final['residue_num'].astype(str)
    df_DSSP['residue_num'] = df_DSSP['residue_num'].astype(str)
    df_final = df_final.merge(df_DSSP, on=['chain', 'residue_num'], how='left')
    

    # Save the final DataFrame to a CSV file
    df_final.to_csv(f"./{OUT_PATH}/{pdb_file[:-4]}.csv", index=False)

# Main function to process all PDB files
def process_all_pdb_files():
    # Get a list of all PDB files in the current directory
    pdb_files = glob("*.pdb")
                
    # Use multiprocessing to process PDB files in parallel
    pool = multiprocessing.Pool(processes=4) # Create a pool with 12 processes
    pool.map(process_pdb_file, pdb_files)     # Map the processing function to all files
    pool.close()                              # Close the pool
    pool.join()                               # Wait for all processes to complete                                       
   

# Run the main function
if __name__ == "__main__":
    # Create output folder
    if not os.path.exists(OUT_PATH):
        os.mkdir(OUT_PATH)
    
    process_all_pdb_files()
    
    print("""
        ________________________
       |                        |
       |    May the Feature     |
       |      be with you.      |
       |________________________|
                 |  |
                 |  |
                 |  |
                 |  |
             ____|  |____
            |____________|
    """)

    print("Thank you for using ChemFeatX, please cite us.")

