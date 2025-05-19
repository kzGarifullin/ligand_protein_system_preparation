from Bio import PDB
import htmd
from moleculekit.molecule import Molecule
from htmd.ui import *
from htmd.home import home
from os.path import join
from moleculekit.config import config
import json


def system_prepare(protein_path, ligand_folder, ligand_path, build_system_outdir, pH, ligand_frcmod=None, resname='MOL'):
    prot = Molecule(protein_path)
    print(prot)
    prot = systemPrepare(prot, pH=pH)
    prot = autoSegment(prot, sel='protein')
    ligand = Molecule(join(ligand_folder, ligand_path))
    print(ligand)
    ligand.set('segid', 'L')  # The segment ID of each atom.
    # 'L' is the new value being assigned to the 'segid' property. It changes the segment identifier of the ligand to 'L'.
    ligand.set('resname', resname)
    # короче теперь segid это L, а resname теперь BEN, как я понимаю, это для удобства
    mol = Molecule(name='combo')
    mol.append(prot)
    mol.append(ligand)
    mol.reps.add(sel='protein', style='NewCartoon', color='Secondary Structure')
    mol.reps.add(sel='resname '+resname, style='Licorice')
    print(mol.coords.shape)
    print("join(ligand_folder, ligand_path)", join(ligand_folder, ligand_path))
    topos_amber = amber.defaultTopo() + [join(ligand_folder, ligand_path)]
    if ligand_frcmod==None:
        print("join(ligand_folder, 'MOL.frcmod')", join(ligand_folder, 'MOL.frcmod'))
        frcmods_amber = amber.defaultParam() + [join(ligand_folder, 'MOL.frcmod')]
    else:
        print("join(ligand_folder,"+ ligand_frcmod+")", join(ligand_folder, ligand_frcmod))
        frcmods_amber = amber.defaultParam() + [join(ligand_folder, ligand_frcmod)]
    bmol_amber = amber.build(mol,  topo=topos_amber, param=frcmods_amber, outdir=build_system_outdir , ionize = False)
    print(ligand.coords.shape)


def translation_vector(ligand_folder, ligand_path, structure_file, resname='MOL'):
    ligand = Molecule(join(ligand_folder, ligand_path))
    residue_name = resname

    # Get coordinates
    mol_coords = get_coordinates_of_residue(structure_file, residue_name)
    diffs = []
    coord_error_flag = False
    for i in range(len(ligand.coords)):
        diff = mol_coords[i] - ligand.coords[i].T
        diffs.append(diff)
        #print(diffs[0] - diffs[-1])
        #assert np.allclose(diffs[0], diffs[-1], atol=1e-3), f"Coordinates at index {i} are not equal!"
        if not np.allclose(diffs[0], diffs[-1], atol=1e-3):
            print(f"Coordinates at index {i} are not equal!")
            coord_error_flag = True
            print("diffs:", diffs)
    return diffs[0], coord_error_flag

def get_coordinates_of_residue(pdb_file, residue_name):
    '''
    Function to get atom coordinates for a specific residue name
    '''
    parser = PDB.PDBParser()
    structure = parser.get_structure('Structure', pdb_file)

    coords = []

    # Iterating over all residues in all models
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.get_resname() == residue_name:
                    for atom in residue:
                        coords.append(atom.get_coord())

    return coords

def resname_parse(cif_file):
    resnames = set()
    with open(cif_file, "r") as f:
        for line in f:
            if line.startswith("HETATM") or line.startswith("ATOM"):
                parts = line.split()
                if len(parts) > 5:
                    resname = parts[5]
                    resnames.add(resname)
    return list(resnames)[0]

if __name__ == "__main__":
    #protein = "6dz3"
    input_file = "/mnt/ligandpro/data/garifullin/system_preparation_train/PDBBind_train"
    #input_file = "/mnt/ligandpro/data/garifullin/system_preparation/re_train.txt"
    translation_vectors_dict = {}
    trans_vector_error_protein_list = []
    error_protein_list = []
    with open(input_file, "r") as f:
        for line in f:
            protein = line.strip()
            if protein:  # skip empty lines
                #print(pdb_id)
                print(protein)
                resname = resname_parse("train/"+protein+"/"+protein+"_ligand_final.cif")
                print(resname)

                #print(a.shape)
                protein_path = 'train/' + protein + '/' + protein + '_protein_processed.pdb'
                ligand_folder = 'train/' + protein + '/'
                build_system_outdir='build_amber/' + protein
                structure_file='build_amber/' + protein + '/structure.pdb'
                try:
                    system_prepare(protein_path=protein_path, ligand_folder=ligand_folder, ligand_path=protein + '_ligand_final.cif', ligand_frcmod=protein+'_mol0_AC.frcmod', build_system_outdir=build_system_outdir, pH=7.4, resname=resname)
                    trans_vector, coord_error_flag = translation_vector(ligand_folder=ligand_folder, ligand_path=protein + '_ligand_final.cif', structure_file=structure_file, resname=resname)
                    print(trans_vector)
                    if coord_error_flag == False:
                        translation_vectors_dict[protein] = list(trans_vector[0])
                    else:
                        trans_vector_error_protein_list.append(protein)
                except Exception as e:
                    print(f"Error processing {protein}: {e}")
                    error_protein_list.append(protein)
                #break
    print(translation_vectors_dict)
    translation_vectors_dict = {k: [float(x) for x in v] for k, v in translation_vectors_dict.items()}

    with open('translation_vectors.json', 'w') as f:
        json.dump(translation_vectors_dict, f, indent=4)
        
    with open('trans_vector_error_protein_list.txt', 'w') as f:
        for protein in trans_vector_error_protein_list:
            f.write(protein + '\n')
        
    with open('error_protein_list.txt', 'w') as f:
        for protein in error_protein_list:
            f.write(protein + '\n')


