





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






protein_to_complex_names = defaultdict(list)
for name in complex_names_all:
    protein_to_complex_names[name.split('_superlig')[0]].append(name)





for protein_name, protein_complex_names in protein_to_complex_names.items():

    rec_model = parse_receptor(complex_names[0], self.data_dir, self.dataset_type)

    for name in protein_complex_names:
        
        print('complex', name)

        if self.dataset_type == 'pdbbind' or self.dataset_type == 'lpce':
            ligs = read_mols(self.data_dir, name, remove_hs=False)
        elif self.dataset_type == 'moad':
            ligs = [read_molecule(os.path.join(self.data_dir, 'pdb_superligand', f'{name}.pdb'), remove_hs=False, sanitize=True)]
        elif self.dataset_type == 'dockgen' or self.dataset_type == 'dockgen_full':
            ligs = [read_molecule(os.path.join(self.data_dir, name, f'{name}_ligand.pdb'), remove_hs=False, sanitize=True)]
        elif self.dataset_type == 'astex' or self.dataset_type == 'posebusters':
            ligs = [read_molecule(os.path.join(self.data_dir, name, f'{name}_ligand.sdf'), remove_hs=False, sanitize=True)]
        elif self.dataset_type == 'astex_conf' or self.dataset_type == 'posebusters_conf':
            ligs = [read_molecule(os.path.join(self.data_dir, name, f'{name}_ligand_start_conf.sdf'), remove_hs=False, sanitize=True)]
        else:
            raise ValueError(f'Unknown dataset type: {self.dataset_type}')

        ligs = [split_molecule(lig_mol, min_lig_size=7) for lig_mol in ligs]
        ligs = [lig_mol for lig_mol_list in ligs for lig_mol in lig_mol_list if lig_mol is not None]

        for lig_idx, lig_mol in enumerate(ligs):
            if self.max_lig_size is not None and lig_mol.GetNumHeavyAtoms() > self.max_lig_size:
                print(f'Ligand with {lig_mol.GetNumHeavyAtoms()} heavy atoms is larger than max_lig_size {self.max_lig_size}. Not including {name} in preprocessed data.')
                continue

            try:
                c_alpha_coords_list, lm_embeddings_list, sequences_list, chain_lengths, full_coords, valid_chain_names = extract_receptor_structure_prody(
                        copy.deepcopy(rec_model), lig_mol)
                

            
                
            # TODO extract chains valid_chain_names and save to pdb
            # lig_mol to sdf

            final_name = f'{name}_mol{lig_idx}'

            smiles = names2smiles[final_name]
            frcmod_name = smiles2name[smiles][0]
            

