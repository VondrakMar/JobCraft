import h5py
import numpy as np
from ase import Atoms
from ase.io import write
import sys

import argparse

def load_h5md(path: str) -> list[Atoms]:
    with h5py.File(path, "r") as f:
        atomic_numbers = f["data/atomic_numbers"][0]          # (N,)
        positions_all  = f["data/positions"][:]               # (n_pos, N, 3)
        cell_all       = f["data/cell"][:]                    # (n_pos, 3, 3)
        pbc            = f["data/pbc"][:].astype(bool)        # (3,)
        energy_all     = f["data/potential_energy"][:, 0]     # (n_e,)
        temp_all       = f["data/temperature"][:, 0]          # (n_e,)
        pos_steps      = f["steps/positions"][:]              # (n_pos,)
        # print(positions_all)
        '''
        energy_steps   = f["steps/potential_energy"][:]       # (n_e,)
        ke_all         = f["data/kinetic_energy"][:, 0]       # (n_e,)
    energy_map = dict(zip(energy_steps, energy_all))
    ke_map     = dict(zip(energy_steps, ke_all))
    temp_map   = dict(zip(energy_steps, temp_all))
    '''
    frames = []
    for i, pos in enumerate(pos_steps):
        atoms = Atoms(
            numbers=atomic_numbers,
            positions=positions_all[i],
            cell=cell_all[i],
            pbc=pbc,
        )
        # atoms.info["energy"]         = float(energy_all[i])
        # atoms.info["kinetic_energy"] = float(ke_all[i])
        # atoms.info["temperature"]    = float(temp_all[i])
        frames.append(atoms)

    return frames

def convert_h5_xyz(h5md_path: str, out_path: str | None = None) -> str:
    frames = load_h5md(h5md_path)
    write(out_path, frames, format="extxyz")
    return out_path

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_file",help="Where is the h5 file",type=str)
    parser.add_argument("--output_file",help="Name of the xyz file",type=str)
    args = parser.parse_args()
    convert_h5_xyz(h5md_path=args.input_file,out_path=args.output_file)
