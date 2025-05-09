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
                     xc : str = "pbe",
                     spin : bool =False,
                     max_l : int = 2,
                     save_hirshfeld_volumes=False):
    '''
    currently implemented readings. I was (heavily) inspired by the ASE, but adding reading of hirshfeld charges was easir this way,
    becaues I have a small brain

    calculation_time, hartree_multipoles, fermi_level, VBM, CBM and number_of_scf_cylces are stole from Will Baldwin
    '''
    implemented_properties = ["energy","forces","hirshfeld","hirshfeld_spin","calculation_time","number_of_scf_cylces","fermi_level","fermi_level_up","fermi_level_down","VBM","CBM","hartee_multipoles"]
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
    if "hirshfeld_spin" in properties:
        hirshfeld_q = []
        hirshfeld_spin=[]
        hirshfeld_dipole=[]
        hirshfeld_quadruple=[]
        free_atom_volume = []
        hirshfeld_volume = []
        for n,line in enumerate(output_file):
            if (line.rfind("Performing Hirshfeld analysis of fragment charges and moments.")) >-1:
                # this is tested only for PBE
                count = 0
                if spin:
                    for iatom in range(natoms):
                        hirshfeld_q.append(float(output_file[n + iatom*11 + 3].split()[-1]))
                        free_atom_volume.append(float(output_file[n + iatom*11 + 4].split()[-1]))
                        hirshfeld_volume.append(float(output_file[n + iatom*11 + 5].split()[-1]))
                        hirshfeld_spin.append(float(output_file[n + iatom*11 + 6].split()[-1]))
                        hirshfeld_dipole.append(np.array([float(tmp) for tmp in output_file[n + iatom*11 + 7].split()[-3:]]))
                        tmp_quadr = []
                        tmp_quadr.extend([float(tmp) for tmp in output_file[n + iatom*11 + 9].split()[-3:]])
                        tmp_quadr.extend([float(tmp) for tmp in output_file[n + iatom*11 + 10].split()[-3:]])
                        tmp_quadr.extend([float(tmp) for tmp in output_file[n + iatom*11 + 11].split()[-3:]])
                        hirshfeld_quadruple.append(tmp_quadr)
                elif not spin:
                    print("Hello")
                    for iatom in range(natoms):
                        hirshfeld_q.append(float(output_file[n + iatom*10 + 3].split()[-1]))
                        free_atom_volume.append(float(output_file[n + iatom*10 + 4].split()[-1]))
                        hirshfeld_volume.append(float(output_file[n + iatom*10 + 5].split()[-1]))
                        hirshfeld_dipole.append(np.array([float(tmp) for tmp in output_file[n + iatom*10 + 6].split()[-3:]]))
                        tmp_quadr = []
                        tmp_quadr.extend([float(tmp) for tmp in output_file[n + iatom*10 + 8].split()[-3:]])
                        tmp_quadr.extend([float(tmp) for tmp in output_file[n + iatom*10 + 9].split()[-3:]])
                        tmp_quadr.extend([float(tmp) for tmp in output_file[n + iatom*10 + 10].split()[-3:]])
                        hirshfeld_quadruple.append(tmp_quadr)
        if save_hirshfeld_volumes:
            results["free_atom_volume"] = np.array(free_atom_volume)
            results["hirshfeld_volume"] = np.array(hirshfeld_volume)
        if spin:
            results["hirshfeld_spin"] = np.array(hirshfeld_spin)
        results["hirshfeld_q"] = np.array(hirshfeld_q)
        results["hirshfeld_dipole"] = np.array(hirshfeld_dipole)
        results["hirshfeld_quadruple"] = np.array(hirshfeld_quadruple)
    if "hirshfeld" in properties and "hirshfeld_spin" not in properties:
        hirshfeld = []
        for n,line in enumerate(output_file):
            if (line.rfind("Performing Hirshfeld analysis of fragment charges and moments.")) >-1:
                count = 0
                if xc == "pbe":
                    for iatom in range(natoms):
                        if spin:
                            data = output_file[n + iatom*11 + 3].split()
                        else:
                            data = output_file[n + iatom*10 + 3].split()
                        hirshfeld.append(float(data[-1]))                    
                elif xc == "pbesol":
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

