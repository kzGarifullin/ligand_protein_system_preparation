#!/bin/bash

./copy_protein_folders.sh

conda activate pymol_env
./convert_mol2_to_cif.sh
conda_deactivate

./reorder_cif_by_mol2.sh
./renumber_cif_atoms.sh
./cif_prepare.sh
./cif_prepare_mol.sh
./copy_ligand_parametrization.sh

conda activate htmd
nohup python system_preparation.py > output.log 2>&1 &
conda_deactivate