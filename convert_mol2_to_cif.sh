#!/bin/bash

# Read each protein name from the file and copy its folder
while IFS= read -r protein_name; do
    if [ -d "train/$protein_name" ]; then
        echo "Modifying $protein_name..."
        #rm train/$protein_name/$protein_name"_ligand_modified.mol2"
        python convert_mol2_to_cif.py train/$protein_name/$protein_name"_ligand.mol2" train/$protein_name/$protein_name"_ligand_modified.cif"
    else
        echo "Warning: Folder for $protein_name not found in train"
    fi
done < "PDBBind_train"

echo "cif file creation completed!" 