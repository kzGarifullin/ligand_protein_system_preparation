#!/bin/bash

# Source directory
mkdir train
SOURCE_DIR="/mnt/ligandpro/data/PDBBind_processed"
# Destination directory (current directory)
DEST_DIR="train"

# Read each protein name from the file and copy its folder
while IFS= read -r protein_name; do
    if [ -d "$SOURCE_DIR/$protein_name" ]; then
        echo "Copying $protein_name..."
        cp -r "$SOURCE_DIR/$protein_name" "$DEST_DIR/"
    else
        echo "Warning: Folder for $protein_name not found in $SOURCE_DIR"
    fi
done < "PDBBind_train"

echo "Copy operation completed!" 