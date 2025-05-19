#!/bin/bash

mkdir test_combined
cd test_combined
cp ../PDBBind_test .

while IFS= read -r protein_name; do
    if [ -d "/mnt/ligandpro/data/garifullin/system_preparation_test/test/$protein_name" ]; then
        echo "Copying $protein_name..."
        mkdir $protein_name
        cd $protein_name
        cp "/mnt/ligandpro/data/garifullin/system_preparation_test/test/$protein_name/${protein_name}_ligand_modified_new_new.cif"  "${protein_name}_ligand.cif" 
        cp "/mnt/ligandpro/data/garifullin/system_preparation_test/test/$protein_name/${protein_name}_protein_processed.pdb" .
        cp "/mnt/ligandpro/data/garifullin/system_preparation_test/test/$protein_name/${protein_name}_ligand.sdf" .
        cp "/mnt/ligandpro/data/garifullin/system_preparation_test/test/$protein_name/${protein_name}_ligand.mol2" .
        cp "/mnt/ligandpro/data/garifullin/system_preparation_test/build_amber/$protein_name/structure.pdb" .
        cp "/mnt/ligandpro/data/garifullin/system_preparation_test/build_amber/$protein_name/structure.prmtop" .
        cd ..
    else
        echo "Warning: Folder for $protein_name not found in $SOURCE_DIR"
    fi
done < "PDBBind_test"

echo "Copy operation completed!" 