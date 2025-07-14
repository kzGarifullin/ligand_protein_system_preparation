#!/bin/bash

# Read each protein name from the file and copy its folder
while IFS= read -r protein_name; do
    if [ -d "train/$protein_name" ]; then
        echo "Modifying $protein_name..."
        #python reorder_cif_by_mol2.py /mnt/ligandpro/data/garifullin/system_preparation_train_right/train/$protein_name/$protein_name"_ligand.mol2" /mnt/ligandpro/data/garifullin/system_preparation_train_right/train/$protein_name/$protein_name"_ligand_modified.cif" /mnt/ligandpro/data/garifullin/system_preparation_train_right/train/$protein_name/$protein_name"_ligand_reordered.cif"
        #python reorder_cif_by_mol2.py /mnt/ligandpro/data/garifullin/ligand_protein_system_preparation/output4fem_1_superlig_0_mol0.sdf /mnt/ligandpro/data/garifullin/ligand_protein_system_preparation/output4fem_1_superlig_0_mol0_modified.cif output4fem_1_superlig_0_mol0_modified_reordered.cif
        python reorder_cif_by_mol2.py train/$protein_name/$protein_name".sdf" train/$protein_name/$protein_name"_ligand_modified.cif" train/$protein_name/$protein_name"_ligand_reordered.cif"
    
    else
        echo "Warning: Folder for $protein_name not found in train"
    fi
done < "MOAD_train"

echo "cif file creation completed!" 