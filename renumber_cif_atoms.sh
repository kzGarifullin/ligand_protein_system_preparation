#!/bin/bash

# Read each protein name from the file and copy its folder
while IFS= read -r protein_name; do
    if [ -d "train/$protein_name" ]; then
        echo "Modifying $protein_name..."
        #rm train/$protein_name/$protein_name"_ligand_modified.mol2"
        python renumber_cif_atoms.py /mnt/ligandpro/data/garifullin/system_preparation_train_right/train/$protein_name/$protein_name"_ligand_reordered.cif" /mnt/ligandpro/data/garifullin/system_preparation_train_right/train/$protein_name/$protein_name"_ligand_reordered_renumbered.cif"
    else
        echo "Warning: Folder for $protein_name not found in train"
    fi
done < "PDBBind_train"

echo "cif file creation completed!" 