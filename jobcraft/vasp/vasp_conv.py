from ase.io import read, write
import numpy as np
import argparse

def readEF(folders_list):
    len_mols = []
    conv_forces = []
    conv_energies = []
    for folder in folders_list:
        mol = read(f"{folder}/OUTCAR@-1")
        len_mols.append(len(mol))
        conv_forces.append(mol.get_forces())
        conv_energies.append(mol.get_potential_energy())
    return len_mols,conv_energies,conv_forces
        # do your work here

def check_conv(len_mols,cE, cF,folders):
    N_steps = len(cE)
    
    for i in reversed(range(1, N_steps)):
        dE = cE[i]/len_mols[i] - cE[i-1]/len_mols[i-1]
        dF = np.max(np.abs(cF[i] - cF[i-1]))        
        print(f"{folders[i]}-{folders[i-1]}","dE meV/atom",dE*1000)
        print(f"{folders[i]}-{folders[i-1]}","dF meV/A",dF*1000)

        
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--conv_folders",
        nargs="+",
        help="List of folders to iterate through",
        required=True,
    )
    
    args = parser.parse_args()
    folders = args.conv_folders
    len_mols,cE,cF = readEF(folders)
    check_conv(len_mols,cE,cF,folders)
    

    
    
