import argparse

def renumber_cif_atoms(input_cif, output_cif):
    with open(input_cif) as fin, open(output_cif, "w") as fout:
        atom_counter = 1
        for line in fin:
            if line.startswith("HETATM"):
                parts = line.split()
                # Replace the atom serial number (second column) with the new number
                parts[1] = str(atom_counter)
                # Rebuild the line, preserving spacing as much as possible
                new_line = "{:<6}{:>4} ".format(parts[0], parts[1]) + " ".join(parts[2:]) + "\n"
                fout.write(new_line)
                atom_counter += 1
            else:
                fout.write(line)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Renumber atom serials in CIF HETATM records.")
    parser.add_argument("input_cif", help="Input CIF file")
    parser.add_argument("output_cif", help="Output CIF file")
    args = parser.parse_args()
    renumber_cif_atoms(args.input_cif, args.output_cif) 