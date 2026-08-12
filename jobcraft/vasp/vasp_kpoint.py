import numpy as np
import ase.io
import argparse



def k_point_adjuster(ase_atoms, kspacing=0.11, slab=False):
    '''
    This function is stolen from
    https://github.com/Felixrccs/Cu-111-oxide/blob/main/wrapper/vasp_dft.py
    '''
    b = ase_atoms.get_reciprocal_cell()
    k = np.ceil(np.linalg.norm(b, axis=1) * 2 * np.pi / kspacing).astype(int)
    k = [1 if i == 0 else i for i in k]
    if slab:
        k[2] = 1
    return np.array(k)

def write_kpoints(filename="KPOINTS", kpoints=[1,1,1], slab=True, gamma=True):
    # kpts = k_point_adjuster(atoms, kspacing=kspacing, slab=slab)
    with open(filename, "w") as f:
        f.write("KPOINTS generated from kspacing\n")
        f.write("0\n")  # automatic grid
        f.write("Gamma\n" if gamma else "Monkhorst-Pack\n")
        f.write(f"{kpoints[0]} {kpoints[1]} {kpoints[2]}\n")
        f.write("0 0 0\n")
    print(f"KPOINTS written: {kpoints}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_file",help="file for which calculate k-grid",type=str)
    parser.add_argument("--output_points",help="k grid file name",default="KPOINTS",type=str)
    parser.add_argument("--kspacing",help="k spacing in reciprocal A",default=0.11,type=float)
    # parser.add_argument("--k_grid",help="k grid, if kspacing is selected this is ignored",nargs="+"default=0.11,type=float)
    parser.add_argument("--gamma_not_centered",help="if gamma should be center",action="store_false")
    parser.add_argument("--is_slab",help="If true, k grid in z direction is 1",action="store_true")
    args = parser.parse_args()
    try:
        # this is here in case you want to read extxyz, but your file is xyz
        mols = ase.io.read(args.input_file,format="extxyz")
    except:
        mols = ase.io.read(args.input_file)
    k_points = k_point_adjuster(ase_atoms=mols,kspacing=args.kspacing,slab=args.is_slab)
    write_kpoints(mols, filename="KPOINTS", kpoints=k_points, gamma=args.gamma_not_centered)
