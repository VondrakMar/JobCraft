import ase
import ase.io
from ase.data import covalent_radii,atomic_numbers 
import numpy as np
from ase.visualize import view
from copy import deepcopy

def align_longest_axis(atoms, target_axis=np.array([1.0, 0.0, 0.0])):
    """
    Rotates the molecule so that its longest spatial extent aligns with target_axis.
    
    Parameters
    ----------
    atoms : ase.Atoms
        Input structure (modified in-place)
    target_axis : array-like
        Axis to align to (default: x-axis)
    """
    tmp_atoms = deepcopy(atoms)
    positions = tmp_atoms.get_positions()
    center = positions.mean(axis=0)
    positions -= center
    cov = np.dot(positions.T, positions)
    eigvals, eigvecs = np.linalg.eigh(cov)
    principal_axis = eigvecs[:, np.argmax(eigvals)]
    principal_axis /= np.linalg.norm(principal_axis)
    target_axis = target_axis / np.linalg.norm(target_axis)
    v = np.cross(principal_axis, target_axis)
    s = np.linalg.norm(v)
    c = np.dot(principal_axis, target_axis)
    if np.isclose(s, 0):
        
        return tmp_atoms
    vx = np.array([
        [0, -v[2], v[1]],
        [v[2], 0, -v[0]],
        [-v[1], v[0], 0]
    ])
    R = np.eye(3) + vx + np.dot(vx, vx) * ((1 - c) / (s**2))
    rotated_positions = positions @ R.T
    rotated_positions += center
    tmp_atoms.set_positions(rotated_positions)

    return tmp_atoms

def get_sphere(solvent,factor=1.5) -> float:
    axis_vector = np.array([1.0,0.0,0.0]) 
    rotated_solv_base = align_longest_axis(solvent,axis_vector)
    pos = rotated_solv_base.positions
    x_min_id,y_min_id,z_min_id,x_max_id,y_max_id,z_max_id = pos[:, 0].argmin(),pos[:, 1].argmin(),pos[:, 2].argmin(),pos[:, 0].argmax(),pos[:, 1].argmax(),pos[:, 2].argmax()
    x_min,y_min,z_min,x_max,y_max,z_max = pos[:, 0].min(),pos[:, 1].min(),pos[:, 2].min(),pos[:, 0].max(),pos[:, 1].max(),pos[:, 2].max()
    dist_x = abs(x_max - x_min)
    dist_y = abs(y_max - y_min)
    dist_z = abs(z_max - z_min)
    sphere_rad = max([dist_x,dist_y,dist_z])
    
    return factor*sphere_rad

def random_rotate_atoms(atoms: ase.Atoms, rng=None) -> ase.Atoms:
    if rng is None:
        rng = np.random.default_rng()

    mol = atoms.copy()
    center = mol.get_center_of_mass()
    angles = rng.uniform(0.0, 360.0, size=3)
    mol.rotate(angles[0], "x", center=center)
    mol.rotate(angles[1], "y", center=center)
    mol.rotate(angles[2], "z", center=center)

    return mol

def is_inside_cell(atoms):
    scaled = atoms.get_scaled_positions(wrap=False)
    return ((scaled >= 0.0) & (scaled < 1.0)).all()

def solvate_slab(slab, solvent, solv_axis=2, slab_element=None, N_mols=20,sphere_factor=1.5):
    sphere_rad = get_sphere(solvent,sphere_factor)

    if slab_element is None:
        el_max_slab = slab.positions[:, solv_axis].argmax()
    else:
        el_max_slab = slab_element
    
    elnum_max_slab = atomic_numbers[el_max_slab]
    dir_max_slab = slab.positions[:, solv_axis].max()
    z_offset_base = dir_max_slab + sphere_rad + 0.7 * covalent_radii[elnum_max_slab]

    cell = slab.cell
    x_max = cell[0, 0]
    y_max = cell[1, 1]

    slab_solv = deepcopy(slab)
    placed = 0
    max_attempts_per_mol = 200
    
    x_pos = 0.5
    y_pos = 1.0
    z_layer = 0  

    for admolecule in range(N_mols):
        placed_this_mol = False
        attempts = 0

        while not placed_this_mol and attempts < max_attempts_per_mol:
            tmp_solv = deepcopy(solvent)
            tmp_solv = random_rotate_atoms(tmp_solv)

            z_trans = z_offset_base + z_layer * sphere_rad
            tmp_solv.translate([x_pos, y_pos, z_trans])

            tmp_slab_solv = slab_solv + tmp_solv

            if is_inside_cell(tmp_slab_solv):
                slab_solv = tmp_slab_solv
                placed_this_mol = True
                
                x_pos += sphere_rad
                if x_pos + sphere_rad > x_max:
                    x_pos = 0.5
                    y_pos += sphere_rad
                    if y_pos + sphere_rad > y_max:
                        y_pos = 1.0
                        z_layer += 1
            else:
                
                x_pos += sphere_rad
                if x_pos + sphere_rad > x_max:
                    x_pos = 0.5
                    y_pos += sphere_rad
                if y_pos + sphere_rad > y_max:
                    y_pos = 1.0
                    z_layer += 1
                attempts += 1

        if not placed_this_mol:
            print(f"Warning: could not place molecule {admolecule}")
            break

    return slab_solv
