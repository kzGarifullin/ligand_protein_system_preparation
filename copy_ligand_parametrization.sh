#!/bin/bash

# Source directory
SOURCE_DIR="/mnt/ligandpro/data/dfrolova/flowdock_data/torsion_histograms/intermediate_files"
# Destination directory (current directory)
DEST_DIR="train"

# Read each protein name from the file and copy its folder
while IFS= read -r protein_name; do
    if [ -d "$SOURCE_DIR/${protein_name}_mol0.acpype" ]; then
        echo "Copying parametrization of $protein_name..."
        cp "$SOURCE_DIR/${protein_name}_mol0.acpype/"*.frcmod "$DEST_DIR/$protein_name/"
    else
        echo "Warning: Folder for $protein_name not found in $SOURCE_DIR"
    fi
done < "PDBBind_train"

echo "Copy operation completed!" 