import pickle
import os
from rdkit import Chem
from ase import Atoms
from ase.io import write
import shutil

# Load the data
cache_path = "/mnt/ligandpro/data/dfrolova/flowdock_data/cache/moad_None_None_hf_esm_12_embeddings_pt_MOAD_PDBBind_txt"
filepath = os.path.join(cache_path, 'complexes.pkl')

with open(filepath, 'rb') as f:
    complexes = pickle.load(f)

print(f"Data successfully loaded from {filepath}!")
print("len(complexes)", len(complexes))

count = 0

file_name = '/mnt/ligandpro/data/dfrolova/flowdock_data/data/splits/MOAD_PDBBind.txt'
complex_names_all = []
with open(file_name, 'r') as file:
    for line in file:
        stripped_line = line.strip()
        complex_names_all.append(stripped_line)

print(complex_names_all[:10])
print(len(complex_names_all))   #40772
print(len(complexes))           #42713

print(b.shape)
# Iterate over each complex
for i in range(len(complexes)):
    print(complexes[i].name)
    #print(complexes[i].name.split('_mol')[0])
    
    # Access the ligand structure
    ligand = complexes[i].ligand.orig_mol
    #print(ligand.GetNumAtoms())
    
    # Prepare the directory name based on the complex name
    directory_name = f'train/{complexes[i].name}'

    # Create the directory if it does not exist
    os.makedirs(directory_name, exist_ok=True)  # This creates the directory and does nothing if it already exists
    
    # Define the SDF file path
    sdf_file_path = f'sdf_mols/{complexes[i].name}.sdf'
    pdb_file_path = f"closest_chains_pdbs/{complexes[i].name.split('_mol')[0]}.pdb"
    
    # Copy the SDF file to the new directory
    if os.path.exists(sdf_file_path) and os.path.exists(pdb_file_path):
        new_file_path_sdf = os.path.join(directory_name, f'{complexes[i].name}.sdf')
        new_file_path_pdb = os.path.join(directory_name, f"{complexes[i].name.split('_mol')[0]}.pdb")
        # os.rename(sdf_file_path, new_file_path_sdf)  # or use shutil.copy to copy instead if you want to keep the original
        # os.rename(pdb_file_path, new_file_path_pdb)  # or use shutil.copy to copy instead if you want to keep the original
        shutil.copy(sdf_file_path, new_file_path_sdf)  
        shutil.copy(pdb_file_path, new_file_path_pdb)  
        #print(f"Copied {sdf_file_path} to {new_file_path_sdf}")
        #print(f"Copied {pdb_file_path} to {new_file_path_pdb}")
    else:
        print(f"Warning: {sdf_file_path} or {pdb_file_path} does not exist!")

    # count += 1
    # if count ==10:
    #     break
