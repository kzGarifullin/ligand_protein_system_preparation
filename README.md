# Ligand-protein system preparation and parametrization

## Setup

* Clone this repo:
```bash
git clone https://github.com/kzGarifullin/ligand_protein_system_preparation.git
cd MaterialFusion
```

* Setup the pymol environment. Conda environment `pymol_env` will be created and you can use it.
```bash
conda create --name pymol_env
conda activate pymol_env
conda install conda-forge::pymol-open-source
conda deactivate
```

* Setup the htmd environment. Conda environment `htmd` will be created and you can use it.
```bash
conda create --name hemd
conda activate htmd
conda install acellera::htmd
conda deactivate
```

## Quickstart

To prepare and parametrize ligand-protein systems run:
```bash
./run.sh
```

## Scripts description

- **copyproteinfolders.sh**: Copies source files of ligand-protein systems from PDBBind.
- **convertmol2tocif.sh**: Converts .mol2 files to .cif format.
- **reordercifbymol2.sh**: Rearranges the order of atoms in the generated SIF files to match the order in the MOL2 files, ensuring consistency in atom arrangement.
- **renumbercifatoms.sh**: Renumbers atoms in the generated SIF files.
- **cifprepare.sh**: Prepares correct atom IDs for all atoms from the generated SIF file.
- **cifpreparemol.sh**: Renames the residues to "MOL".
- **copyligandparametrization.sh**: Copies .frcmod files if they exist.
- **systempreparation.py**: Prepares the parametrization of the ligand-protein system.
