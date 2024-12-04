import numpy as np
import ase.io
import os
from typing import Any, Callable, Dict, List
from .hartree_multipoles import read_multipoles_from_output_file,change_conventions,reverse_search_for

def read_aims_output(mol: str =None,
                     mol_file_name : str = None,
                     mol_file_format : str ="extxyz",
                     output_name : str="aims.out",
                     properties : List[str]=["energy"],
                     xc : str = None,
                     spin : bool =False,
                     max_l : int = 2):
    '''
    currently implemented readings. I was (heavily) inspired by the ASE, but adding reading of hirshfeld charges was easir this way,
    becaues I have a small brain

    calculation_time, hartree_multipoles, fermi_level, VBM, CBM and number_of_scf_cylces are stole from Will Baldwin
    '''
    implemented_properties = ["energy","forces","hirshfeld","calculation_time","number_of_scf_cylces","fermi_level","fermi_level_up","fermi_level_down","VBM","CBM","hartee_multipoles"]
    for property in properties:
        assert property in implemented_properties, f"{property} not implemented, please use one of the {implemented_properties}"
    if mol_file_name == None:
        assert type(mol) != list, "mol should be a ase.atoms object not a list"
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
                if xc == "pbesol":
                    for iatom in range(natoms):
                        if spin:
                            data = output_file[n + iatom*11 + 7].split()
                        else:
                            data = output_file[n + iatom*10 + 7].split()
                        hirshfeld.append(float(data[-1]))
                else:
                    for iatom in range(natoms):
                        if spin:
                            data = output_file[n + iatom*11 + 3].split()
                        else:
                            data = output_file[n + iatom*10 + 3].split()
                        hirshfeld.append(float(data[-1]))
        results['hirshfeld'] = np.array(hirshfeld)
    if "calculation_time" in properties:
        line_start = reverse_search_for(output_file, ["Detailed time accounting"]) + 1
        print("time",line_start)
        tot_time = float(output_file[line_start].split(":")[-1].strip().split()[0])
        results['calculation_time'] = tot_time
    if "number_of_scf_cylces" in properties:
        line_start = reverse_search_for(output_file, ["| Number of self-consistency cycles"])
        num_scf = float(output_file[line_start].split(":")[-1].strip().split()[0])
        results['number_of_scf_cylces'] = num_scf
    if "VBM" in properties:
        line_start = reverse_search_for(output_file, ["Highest occupied state (VBM) at"])
        e_vbm = float(output_file[line_start].split()[-6])
        results['VBM'] = e_vbm
    if "CBM" in properties:
        line_start = reverse_search_for(output_file, ["Lowest unoccupied state (CBM) at"])
        e_cbm = float(output_file[line_start].split()[-6])
        results['CBM'] = e_cbm
    if "fermi_level" in properties:
        line_start = reverse_search_for(output_file, ["| Chemical potential (Fermi level):"])
        e_fermi = float(output_file[line_start].split(":")[-1].strip().split()[0])
        results['fermi_level'] = e_fermi
    if "fermi_level_up" in properties:
        line_start = reverse_search_for(output_file, ["| Chemical Potential, spin up"])
        e_fermi = float(output_file[line_start].split(":")[-1].strip().split()[0])
        results['fermi_level_up'] = e_fermi
    if "fermi_level_down" in properties:
        line_start = reverse_search_for(output_file, ["| Chemical Potential, spin down"])
        e_fermi = float(output_file[line_start].split(":")[-1].strip().split()[0])
        results['fermi_level_down'] = e_fermi
    if "hartee_multipoles" in properties:
        multipoles = read_multipoles_from_output_file(output_name, mol.get_global_number_of_atoms(), max_l)
        results['atomic_multipoles'] = change_conventions(multipoles, max_l)
    return results

