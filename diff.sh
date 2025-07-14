#!/bin/bash

# Read each protein name from the file and copy its folder
while IFS= read -r protein_name; do
    if [ -d "train/$protein_name" ]; then
        #echo "Modifying $protein_name..."
        #diff train/$protein_name/$protein_name"_ligand_modified.cif" train/$protein_name/$protein_name"_ligand_reordered.cif"
        if ! diff train/$protein_name/$protein_name"_ligand_modified.cif" train/$protein_name/$protein_name"_ligand_reordered.cif" > /dev/null; then
            echo "Files differ for: $protein_name"
            diff train/$protein_name/$protein_name"_ligand_modified.cif" train/$protein_name/$protein_name"_ligand_reordered.cif"
        fi
    else
        echo "Warning: Folder for $protein_name not found in train"
    fi
done < "MOAD_train"

echo "checked!" 