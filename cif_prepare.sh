#!/bin/bash


# Read each protein name from the file and copy its folder
while IFS= read -r protein_name; do
    if [ -d "train/$protein_name" ]; then
        echo "Modifying $protein_name..."
        python cif_prepare.py train/$protein_name/$protein_name"_ligand_reordered_renumbered.cif" train/$protein_name/$protein_name"_ligand_reordered_renumbered_prepared.cif"
    else
        echo "Warning: Folder for $protein_name not found in train"
    fi
done < "PDBBind_train"

echo "Copy operation completed!" 