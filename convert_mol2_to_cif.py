import pymol
import argparse


def main():
    #python mol2_prepare.py input.mol2 output.mol2
    parser = argparse.ArgumentParser(description='Modify atom names in MOL2 file based on element types')
    parser.add_argument('input_file', help='Input MOL2 file path')
    parser.add_argument('output_file', help='Output MOL2 file path')
    
    args = parser.parse_args()
    
    pymol.cmd.load(args.input_file, 'molecule')
    pymol.cmd.save(args.output_file, 'molecule')
    print("Conversion completed: "+args.input_file+" -> "+args.output_file) 

if __name__ == "__main__":
    main()
