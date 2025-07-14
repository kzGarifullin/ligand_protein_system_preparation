# Ligand-protein system preparation and parametrization

The main branch of this repository is for parameterizing the PDBBind dataset, while the moad_parameterization branch is for parameterizing the MOAD dataset.

## Setup

* Clone this repo:
```bash
git clone https://github.com/kzGarifullin/ligand_protein_system_preparation.git
cd ligand_protein_system_preparation
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

## Quickstart for PDBBind dataset

To prepare and parametrize ligand-protein systems run:
```bash
./run.sh
```

## Scripts description for PDBBind dataset parametrization

- **copy_protein_folders.sh**: Copies source files of ligand-protein systems from PDBBind.
- **convert_mol2_to_cif.sh**: Converts .mol2 files to .cif format.
- **reorder_cif_by_mol2.sh**: Rearranges the order of atoms in the generated SIF files to match the order in the MOL2 files, ensuring consistency in atom arrangement.
- **renumber_cif_atoms.sh**: Renumbers atoms in the generated SIF files.
- **cif_prepare.sh**: Prepares correct atom IDs for all atoms from the generated SIF file.
- **cif_prepare_mol.sh**: Renames the residues to "MOL".
- **copy_ligand_parametrization.sh**: Copies .frcmod files if they exist.
- **system_preparation.py**: Prepares the parametrization of the ligand-protein system.


## MOAD dataset parametrization

For the MOAD dataset, the sequence of actions is slightly different from that used for PDBBind dataset.

```bash
conda activate flowdock
python moad_ligand_process.py
python parser_pdb.py
python create_folders_moad.py
conda deactivate 
make_systems_file.sh 
conda activate pymol_env
./convert_sdf_to_cif.sh
conda deativate
./reorder_cif_by_mol2.sh
./renumber_cif_atoms.sh
./cif_prepare.sh
./cif_prepare_mol.sh
./copy_ligand_parametrization.sh
conda activate htmd
nohup python system_preparation.py > output.log 2>&1 &
conda_deactivate
```

## Scripts description for MOAD parametrization

- **moad_ligand_process.py** - сохрание молекул в виде сдф из кэша
- **parser_pdb** - сохраняет в отдельные pdb файлы только те части белка, которые лежат недалеко от лиганда
- **create_folders_moad.py** - создаем папку train/ в которой буду храниться наши комплексы. Сразу туда копируем наши sdf и pdb (только ближайшие цепи)
- **make_systems_file.sh**  - создаст файл с названиями систем, которые возникли в папке train
- **convert_sdf_to_cif.sh** - перевод из sdf в cif
- **reorder_cif_by_mol2.sh** - теперь в cif нужно поменять атомы местами так, чтобы их порядок был такой же как и в sdf, это важно потому что далее при расчете энергии параметризация создается из cif! сам же пайплайн генерации использует sdf файлы, поэтому порядок в sdf и cif файлах должен совпадать