import h5py
import numpy as np
from ase import Atoms
from ase.io import write
import sys

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
        print(f["data/potential_energy"].shape)
        '''
        energy_steps   = f["steps/potential_energy"][:]       # (n_e,)
        ke_all         = f["data/kinetic_energy"][:, 0]       # (n_e,)
    energy_map = dict(zip(energy_steps, energy_all))
    ke_map     = dict(zip(energy_steps, ke_all))
    temp_map   = dict(zip(energy_steps, temp_all))
    '''
    print(energy_all)
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




# This was originally written only for ethanol - benzene mixed system, if you reading this I was just lazy to ever do it for a general system
N_at_ethanol = 9
benzene_start = 1305+9
N_at_benzene = 12

def get_separation(mol,supercell_size=(3, 3, 3)):
    mults = supercell_size[0]*supercell_size[1]*supercell_size[2]
    nAts = int(len(mol)/mults)
    ethanol_nAts = benzene_start
    benzene_nAts = nAts - ethanol_nAts
    ethanol_part = []
    benzene_part = []
    for mult in range(mults):
        ethanol_part.extend(range(mult*(nAts),mult*(nAts)+ethanol_nAts))
        benzene_part.extend(range(mult*(nAts)+ethanol_nAts,mult*(nAts)+nAts))
    ethanol = mol[ethanol_part]
    benzene = mol[benzene_part]
    return ethanol,benzene


def get_molecule_list(mol,supercell_size=[3,3,3]):
    mults = supercell_size[0]*supercell_size[1]*supercell_size[2]
    nAts = int(len(mol)/mults)
    moldy_ethen, benzene = get_separation(mol,supercell_size)
    ethanol_nAts = int(benzene_start*mults)
    benzene_nAts = int((nAts - ethanol_nAts)*mults)
    ethanol_list = []
    benzene_list = []
    for id_ethanol in range(0,ethanol_nAts,9):
        ethanol_list.append(moldy_ethen[id_ethanol:id_ethanol+9])
    for id_benzene in range(0,benzene_nAts,12):
        benzene_list.append(benzene[id_benzene:id_benzene+12])
        # ethanol_list[-1].cell = [False,False,False]
        # ethanol_list[-1].center()
    return benzene_list,ethanol_list

def extract_clusters(mol,r_cut = 8.0):
    print("extracting clusters")
    assert mol.cell, "Works only for PBC systems"
    # cluster_list = []
    center = [mol.cell[0][0]/2,mol.cell[1][1]/2,mol.cell[2][2]/2]
    benzene_list,ethanol_list = get_molecule_list(mol,[1,1,1])

    for id_benzene,benzene in enumerate(benzene_list):
        positions = benzene.get_positions()
        distances = np.linalg.norm(positions - center,axis=1)
        is_in_question_mark = distances <= r_cut
        if True in is_in_question_mark:
            try:
                final_cluster += benzene
            except:
                final_cluster = benzene
    for id_ethanol,ethanol in enumerate(ethanol_list):
        positions = ethanol.get_positions()
        distances = np.linalg.norm(positions - center,axis=1)
        is_in_question_mark = distances <= r_cut
        if True in is_in_question_mark:
            try:
                final_cluster += ethanol
            except:
                final_cluster = ethanol
    return final_cluster


if __name__ == "__main__":
    convert_h5_xyz(f"{sys.argv[1]}",out_path="skibidi.xyz")
