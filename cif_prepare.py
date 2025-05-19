from collections import defaultdict
import argparse

def modify_cif_atom_names(input_file, output_file):
    with open(input_file, 'r') as infile:
        lines = infile.readlines()

    element_counters = defaultdict(int)
    modified_lines = []

    for line in lines:
        if line.startswith('HETATM') or line.startswith('ATOM'):
            parts = line.split()
            if len(parts) >= 4:
                element = parts[2]  # 3rd column: type_symbol
                element_counters[element] += 1
                parts[3] = f"{element}{element_counters[element]}"  # 4th column: label_atom_id
                # Reconstruct the line with single spaces between columns
                new_line = ' '.join(parts) + '\n'
                modified_lines.append(new_line)
            else:
                modified_lines.append(line)
        else:
            modified_lines.append(line)

    with open(output_file, 'w') as outfile:
        outfile.writelines(modified_lines)

def main():
    parser = argparse.ArgumentParser(description='Modify atom names in CIF file based on element types')
    parser.add_argument('input_file', help='Input CIF file path')
    parser.add_argument('output_file', help='Output CIF file path')
    args = parser.parse_args()
    modify_cif_atom_names(args.input_file, args.output_file)

if __name__ == "__main__":
    main()


