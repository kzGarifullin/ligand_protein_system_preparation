import Bio
print("Biopython v" + Bio.__version__)

from Bio.PDB import PDBParser
from Bio.PDB import PDBIO

# Parse and get basic information
parser=PDBParser()
path = '/mnt/ligandpro/data/BindingMOAD_2020_processed/pdb_protein/1p44_1_protein.pdb'
protein_1p49 = parser.get_structure('STS', path)
protein_1p49_resolution = protein_1p49.header["resolution"]
#protein_1p49_keywords = protein_1p49.header["keywords"]

print("Sample name: " + str(protein_1p49))
print("Resolution: " + str(protein_1p49_resolution))
#print("Keywords: " + str(protein_1p49_keywords))
print("Model: " + str(protein_1p49[0]))
print(protein_1p49[0].ca)

#initialize IO 
io=PDBIO()

#custom select
class Select():
    def accept_model(self, model):
        return True
    def accept_chain(self, chain):
        return True
    def accept_residue(self, residue):
        return True       
    def accept_atom(self, atom):
        print("atom id:" + atom.get_id())
        print("atom name:" + atom.get_name())
        if atom.get_name() == 'CA':  
            print("True") 
            return True
        else:
            return False

#write to output file
io.set_structure(protein_1p49)
io.save("1p49_out.pdb", Select())
