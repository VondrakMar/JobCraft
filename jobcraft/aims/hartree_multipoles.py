from scipy.special import sph_harm
import numpy as np
import ase.io.cube
import ase.io
from ase.atoms import Atoms
import os
import csv
import ase.io.cube
from copy import deepcopy
'''
courtesy of Will
'''


ATOM_HEADER = ['|', 'Atom']

def reverse_search_for(lines_obj, keys, line_start=0):
    for ll, line in enumerate(lines_obj[line_start:][::-1]):
        if any([key in line for key in keys]):
            return len(lines_obj) - ll - 1        

def search_for(lines_obj, keys, line_start=0):
    for ll, line in enumerate(lines_obj[line_start:]):
        if any([key in line for key in keys]):
            return line_start + ll
        

def change_conventions(imoments, max_l):
    # phase
    moments = deepcopy(imoments)
    for ll in range(max_l+1):
        for m in range(2*ll+1):
            m_prime = m - ll
            if m_prime > 0:
                moments[:,(ll**2) + m] *= (-1)**(m_prime)

    # units and l's
    for ll in range(max_l+1):
        moments[:,(ll**2):((ll+1)**2)] *= (0.529177210)**ll * (2*ll + 1)**0.5

    return moments


        
def read_an_atom(lines, start_index, max_l):
    moments = []
    index=start_index+1
    for offset in range(121):
        if lines[index+offset] == [] or lines[index+offset].split()[:-2] == ATOM_HEADER:
            raise ValueError('unexpected end of multipole block, possibly too high l requested')
        
        l = int(lines[index+offset].split()[1])
        if l == max_l+1:
            break
    
        moments.append(float(lines[index+offset].split()[3]))
    return np.asarray(moments)


def read_multipoles_from_output_file(rundir, num_atoms, max_l):
    with open(rundir) as f:
        lines = f.readlines()

        line_start = reverse_search_for(lines, ["What we write is -moment*sqrt(pi/4) , which gives charges in e etc."])

        all_moments = []
        for ii in range(num_atoms):
            line_start = search_for(lines, ["Atom"], line_start=line_start)
            assert int(lines[line_start].split()[2].replace(':','')) == ii+1
            all_moments.append(read_an_atom(lines, line_start, max_l))
            line_start+=1
        
        all_moments = np.asarray(all_moments)
    return all_moments
