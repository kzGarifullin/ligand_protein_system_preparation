import pickle
import os
from rdkit import Chem
from ase import Atoms
from ase.io import write
cache_path = "/mnt/ligandpro/data/dfrolova/flowdock_data/cache/moad_None_None_hf_esm_12_embeddings_pt_MOAD_PDBBind_txt"

filepath = os.path.join(cache_path, 'complexes.pkl')
with open(filepath, 'rb') as f:
    complexes = pickle.load(f)
print(f"Data successfully loaded from {filepath}!")
print("len(complexes)", len(complexes))
count = 0
for i in range(len(complexes)):
    print(complexes[i].name)
    ligand = complexes[i].ligand.orig_mol
    print(ligand.GetNumAtoms())
    #print(complexes[i].ligand.pos)
    mol2_file_path = 'sdf_mols/'+complexes[i].name+'.sdf'
    w = Chem.SDWriter(mol2_file_path)  
    w.write(ligand)
    w.close()
    count+=1
    # if count==3:
    #     break