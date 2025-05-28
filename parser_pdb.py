from collections import defaultdict
import copy
import os
import warnings

import numpy as np
import torch
from Bio.PDB import PDBParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from rdkit import Chem
from rdkit.Chem.rdchem import BondType as BT
from rdkit.Chem import AllChem, GetPeriodicTable, RemoveHs
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdMolTransforms
from rdkit.Geometry import Point3D
from scipy import spatial

import torch.nn.functional as F
import networkx as nx
import prody

def mol_to_graph(mol):
    # Initialize graph
    G = nx.Graph()
    for i, atom in enumerate(mol.GetAtoms()):
        G.add_node(i)
    for bond in mol.GetBonds():
        start, end = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        G.add_edge(start, end)
    return G

def split_molecule(mol, min_lig_size=7):
    G = mol_to_graph(mol)

    molecule_parts = []
    for atom_indices in nx.connected_components(G):

        # take the connected component
        atoms_to_remove = list(set(G.nodes) - set(atom_indices))
        atoms_to_remove.sort(reverse=True)

        em1 = Chem.EditableMol(copy.deepcopy(mol))
        for atom in atoms_to_remove:
            em1.RemoveAtom(atom)
            
        mol_part = em1.GetMol()
        try:
            Chem.SanitizeMol(mol_part)
        except:
            print('mol_part sanitization failed')
        if mol_part.GetNumAtoms() >= min_lig_size:
            molecule_parts.append(mol_part)
    return molecule_parts

def parse_receptor(pdbid, pdbbind_dir):
    rec = parsePDB(pdbid, pdbbind_dir)
    return rec

def parsePDB(pdbid, pdbbind_dir):
    rec_path = os.path.join(pdbbind_dir, f'{pdbid}_protein.pdb')
    #rec_path = os.path.join(pdbbind_dir, 'pdb_protein', f'{pdbid.split("_superlig")[0]}_protein.pdb')
    return parse_pdb_from_path(rec_path)


def parse_pdb_from_path(path):
    pdb = prody.parsePDB(path)
    return pdb

def read_molecule(molecule_file, sanitize=False, calc_charges=False, remove_hs=False):
    """
    Read a molecular structure from a file and optionally process it.

    This function reads a molecular structure from various file formats and provides options to sanitize the molecule,
    calculate Gasteiger charges, and remove hydrogen atoms.

    Parameters:
    molecule_file (str): Path to the molecular structure file. Supported formats are .mol2, .sdf, .pdbqt, and .pdb.
    sanitize (bool): If True, sanitize the molecule (default: False).
    calc_charges (bool): If True, calculate Gasteiger charges for the molecule (default: False).
    remove_hs (bool): If True, remove hydrogen atoms from the molecule (default: False).

    Returns:
    RDKit.Chem.Mol or None: The RDKit molecule object if the molecule is successfully read and processed, None otherwise.

    Raises:
    ValueError: If the file format is not supported.

    Notes:
    - Sanitization ensures the molecule's valence states are correct and that the structure is reasonable.
    - Gasteiger charges are partial charges used for computational chemistry methods.
    - Removing hydrogen atoms can be useful for simplifying the molecule, though it may lose information.

    Example:
    >>> from rdkit import Chem
    >>> mol = read_molecule('molecule.mol2', sanitize=True, calc_charges=True, remove_hs=True)
    >>> if mol:
    >>>     print(Chem.MolToSmiles(mol))
    """
    if molecule_file.endswith('.mol2'):
        mol = Chem.MolFromMol2File(molecule_file, sanitize=False, removeHs=False)
    elif molecule_file.endswith('.sdf'):
        supplier = Chem.SDMolSupplier(molecule_file, sanitize=False, removeHs=False)
        mol = supplier[0]
    elif molecule_file.endswith('.pdbqt'):
        with open(molecule_file) as file:
            pdbqt_data = file.readlines()
        pdb_block = ''
        for line in pdbqt_data:
            pdb_block += '{}\n'.format(line[:66])
        mol = Chem.MolFromPDBBlock(pdb_block, sanitize=False, removeHs=False)
    elif molecule_file.endswith('.pdb'):
        mol = Chem.MolFromPDBFile(molecule_file, sanitize=False, removeHs=False)
    else:
        raise ValueError('Expect the format of the molecule_file to be '
                         'one of .mol2, .sdf, .pdbqt and .pdb, got {}'.format(molecule_file))

    try:
        if sanitize or calc_charges:
            Chem.SanitizeMol(mol)

        if calc_charges:
            # Compute Gasteiger charges on the molecule.
            try:
                AllChem.ComputeGasteigerCharges(mol)
            except:
                warnings.warn('Unable to compute charges for the molecule.')

        if remove_hs:
            mol = Chem.RemoveHs(mol, sanitize=sanitize)
    except Exception as e:
        print(e)
        print("RDKit was unable to read the molecule.")
        return None

    return mol

def extract_receptor_structure_prody(rec, lig):
    """
    Extract and process the structure of a receptor in the context of its interaction with a ligand.

    This function extracts the atomic coordinates of amino acids in the receptor, particularly focusing on
    backbone atoms (C-alpha, N, and C). It filters out non-amino acid residues and identifies the chains
    that are valid (contain amino acids) and those that are in close proximity to the ligand.

    Parameters:
    rec (Bio.PDB.Structure.Structure): The receptor structure, typically a Bio.PDB structure object.
    lig (rdkit.Chem.Mol): The ligand molecule, typically an RDKit molecule object.
    lm_embedding_chains (list of np.ndarray, optional): Optional embeddings for each chain from a language model.
        If provided, it should have the same number of chains as the receptor structure.

    Returns:
    tuple:
        - rec (Bio.PDB.Structure.Structure): The modified receptor structure with invalid chains removed.
        - c_alpha_coords (np.ndarray): A numpy array of shape (n_residues, 3) containing the C-alpha atom coordinates of
          valid residues.
        - lm_embeddings (np.ndarray or None): A concatenated numpy array of the valid language model embeddings for the chains,
          if lm_embedding_chains is provided. Otherwise, None.
    """
    if lig is not None:
        conf = lig.GetConformer()
        lig_coords = conf.GetPositions()
    seq = rec.ca.getSequence()
    coords = get_coords(rec)

    res_chain_ids = rec.ca.getChids()
    res_seg_ids = rec.ca.getSegnames()
    res_chain_ids = np.asarray([s + c for s, c in zip(res_seg_ids, res_chain_ids)])
    chain_ids = np.unique(res_chain_ids)
    seq = np.array([s for s in seq])

    sequences = []
    valid_chain_names = []
    lm_embeddings = []
    c_alpha_coords = []
    full_coords = []
    min_distances_to_lig = []
    chain_distances = {}
    for i, chain_id in enumerate(chain_ids):
        chain_mask = res_chain_ids == chain_id
        chain_coords = coords[chain_mask]
        nonempty_coords = chain_coords.reshape(-1, 3)
        nonempty_coords = nonempty_coords[np.isnan(nonempty_coords).sum(axis=1) == 0]

        min_dist_to_lig = 0
        if lig is not None:
            distances = np.linalg.norm(lig_coords[None] - nonempty_coords[:, None], axis=-1)
            min_dist_arr = distances.min(axis=0)
            min_dist_to_lig = distances.min()
            chain_distances[chain_id] = min_dist_to_lig

        if min_dist_to_lig < 4.5:
            sequences.append(tokenized_seq)
            lm_embeddings.append(embeddings)
            valid_chain_names.append(chain_id)
            c_alpha_coords.append(chain_coords[:, 1].astype(np.float32))
            full_coords.append(nonempty_coords)
            if lig is not None:
                min_distances_to_lig.append(min_dist_arr)

    if len(c_alpha_coords) == 0:
        print('NO VALID CHAIN!!!')
        print(chain_distances)
        return None, None, None, None, None

    chain_lengths = [len(seq) for seq in sequences]
    c_alpha_coords = np.concatenate(c_alpha_coords, axis=0)  # [n_residues, 3]
    full_coords = np.concatenate(full_coords, axis=0) # [n_protein_atoms, 3]
    lm_embeddings = np.concatenate(lm_embeddings, axis=0)
    sequences = np.concatenate(sequences, axis=0)

    if lig is not None:
        min_distances_to_lig = np.stack(min_distances_to_lig)
        min_distances_to_lig = min_distances_to_lig.min(axis=0)

        distance_cutoff = 5.
        is_buried_threshold = 0.3 # -100
        buried_atoms_mask = min_distances_to_lig <= distance_cutoff
        fraction_buried = buried_atoms_mask.mean()

        if fraction_buried < is_buried_threshold:
            print(f'Ligand is not buried (fraction_buried = {fraction_buried})')
            return None, None, None, None, None

    return c_alpha_coords, lm_embeddings, sequences, chain_lengths, full_coords, valid_chain_names


# complex_names_all read from  /mnt/ligandpro/data/dfrolova/flowdock_data/data/splits/MOAD_PDBBind.txt


file_name = '/mnt/ligandpro/data/dfrolova/flowdock_data/data/splits/MOAD_PDBBind.txt'
complex_names_all = []
with open(file_name, 'r') as file:
    for line in file:
        stripped_line = line.strip()
        complex_names_all.append(stripped_line)


protein_to_complex_names = defaultdict(list)
for name in complex_names_all:
    protein_to_complex_names[name.split('_superlig')[0]].append(name)

print(protein_to_complex_names['6t4c_1'])
print(protein_to_complex_names['4fch_1'])


for protein_name, protein_complex_names in protein_to_complex_names.items():
    print("protein_name:", protein_name)
    print("protein_to_complex_names:", protein_complex_names)

    rec_model = parse_receptor(protein_name, "/mnt/ligandpro/data/BindingMOAD_2020_processed/pdb_protein")
                #parse_receptor(pdbid, pdbbind_dir)
    print("rec_model:", rec_model)
    for name in protein_complex_names:
        
        print('complex', name)

        ligs = [read_molecule(os.path.join('/mnt/ligandpro/data/BindingMOAD_2020_processed', 'pdb_superligand', f'{name}.pdb'), remove_hs=False, sanitize=True)]


        ligs = [split_molecule(lig_mol, min_lig_size=7) for lig_mol in ligs]
        ligs = [lig_mol for lig_mol_list in ligs for lig_mol in lig_mol_list if lig_mol is not None]
        print(ligs)
        max_lig_size = 200
        for lig_idx, lig_mol in enumerate(ligs):
            if max_lig_size is not None and lig_mol.GetNumHeavyAtoms() > max_lig_size:
                print(f'Ligand with {lig_mol.GetNumHeavyAtoms()} heavy atoms is larger than max_lig_size {self.max_lig_size}. Not including {name} in preprocessed data.')
                continue

            #try:
            c_alpha_coords_list, lm_embeddings_list, sequences_list, chain_lengths, full_coords, valid_chain_names = extract_receptor_structure_prody(
                    copy.deepcopy(rec_model), lig_mol)
                
            # except Exception as e:
            #     print(f"An unexpected error occurred: {e}")
                
            # TODO extract chains valid_chain_names and save to pdb
            # lig_mol to sdf

            #final_name = f'{name}_mol{lig_idx}'
   
            # smiles = names2smiles[final_name]
            # frcmod_name = smiles2name[smiles][0]
            

