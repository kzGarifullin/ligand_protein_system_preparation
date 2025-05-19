import re
import argparse

def parse_mol2_atoms(mol2_path):
    atoms = []
    with open(mol2_path) as f:
        in_atom_section = False
        for line in f:
            if line.startswith("@<TRIPOS>ATOM"):
                in_atom_section = True
                continue
            if line.startswith("@<TRIPOS>"):
                in_atom_section = False
            if in_atom_section:
                parts = line.split()
                if len(parts) < 6:
                    continue
                x, y, z = float(parts[2]), float(parts[3]), float(parts[4])
                atoms.append((x, y, z))
    return atoms

def parse_cif_atoms_and_header(cif_path):
    header = []
    atoms = []
    footer = []
    in_atom_block = False
    with open(cif_path) as f:
        for line in f:
            if line.startswith("HETATM"):
                in_atom_block = True
                parts = line.split()
                # Координаты в CIF: Cartn_x (10), Cartn_y (11), Cartn_z (12) (индексы 10,11,12 или 9,10,11 если считать с 0)
                x, y, z = float(parts[10]), float(parts[11]), float(parts[12])
                atoms.append({'x': x, 'y': y, 'z': z, 'line': line, 'used': False})
            elif in_atom_block and not line.startswith("HETATM"):
                in_atom_block = False
                footer.append(line)
            if not in_atom_block:
                header.append(line)
    # Удаляем строки атомов из header
    header = [l for l in header if not l.startswith("HETATM")]
    return header, atoms, footer

def atoms_match(a, b, tol=0.01):
    return abs(a[0] - b['x']) <= tol and abs(a[1] - b['y']) <= tol and abs(a[2] - b['z']) <= tol

def find_and_mark_atom(mol2_atom, cif_atoms, tol=0.01):
    for atom in cif_atoms:
        if not atom['used'] and atoms_match(mol2_atom, atom, tol):
            atom['used'] = True
            return atom['line']
    return None

def main():
    #python mol2_prepare.py input.mol2 output.mol2
    parser = argparse.ArgumentParser(description='Modify atom names in MOL2 file based on element types')
    parser.add_argument('mol2_path', help='Input MOL2 file path')
    parser.add_argument('cif_path', help='Input CIF file path')
    parser.add_argument('output_cif_path', help='Output CIF file path')
    args = parser.parse_args()
    print(args.mol2_path)

    TOL = 0.01  # допуск по координатам
    print(TOL)
    # 1. Получаем координаты атомов из mol2
    mol2_atoms = parse_mol2_atoms(args.mol2_path)

    # 2. Получаем атомы и header из CIF
    cif_header, cif_atoms, cif_footer = parse_cif_atoms_and_header(args.cif_path)

    # 3. Записываем новый CIF
    with open(args.output_cif_path, "w") as out:
        for line in cif_header:
            out.write(line)
        for mol2_atom in mol2_atoms:
            cif_line = find_and_mark_atom(mol2_atom, cif_atoms, TOL)
            if cif_line:
                out.write(cif_line)
            else:
                print(f"Warning: Atom with coords {mol2_atom} not found in CIF!\n")
        for line in cif_footer:
            out.write(line) 

    print("Conversion completed: "+args.cif_path+" -> "+args.output_cif_path) 

if __name__ == "__main__":
    main()