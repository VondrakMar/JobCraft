import numpy as np
import ase.io
import os


def read_aims_output(mol=None,mol_file_name=None,mol_file_format="extxyz",output_name="aims.out",properties=["energy"]):
    '''
    currently implemented readings. I was (heavily) inspired by the ASE, but adding reading of hirshfeld charges was easir this way,
    becaues I have a small brain
    '''
    if mol_file_name == None:
        assert type(mol) != list
        natoms = len(mol)
    elif mol == None:
        mol = ase.io.read(mol_file_name,format=mol_file_format)
        natoms = len(mol)
    elif mol == None and mol_file_name == None:
        print("specify name of the file to load or ase.Atoms object")
        exit()
    elif mol != None and mol_file_name != None:
        print("You set up both, mol and mol_file_name, assuming that mol is what to use")
        assert type(mol) != list
        natoms = len(mol)
    
    output_file = open(output_name,"r").readlines()
    results = {}
    if "energy" in properties:
        for line in output_file:
            if line.rfind('Total energy corrected') > -1:
                E0 = float(line.split()[5])
        results["energy"] = E0

    if "forces" in properties:
        forces = np.zeros([natoms, 3])
        for n, line in enumerate(output_file):
            if line.rfind('Total atomic forces') > -1:
                for iatom in range(natoms):
                    data = output_file[n + iatom + 1].split()
                    for iforce in range(3):
                        forces[iatom, iforce] = float(data[2 + iforce])
        results['forces'] = forces
    
    if "hirshfeld" in properties:
        hirshfeld = []
        for n,line in enumerate(output_file):
            if (line.rfind("Performing Hirshfeld analysis of fragment charges and moments.")) >-1:
                count = 0
                for iatom in range(natoms):
                    data = output_file[n + iatom*10 + 3].split()
                    hirshfeld.append(float(data[-1]))
        results['hirshfeld'] = np.array(hirshfeld)

    return results
