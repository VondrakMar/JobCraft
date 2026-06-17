from ase import Atoms
from ase.build import molecule
import numpy as np

import numpy as np
from ase import Atoms
from ase.geometry import get_distances


def expand_box_rigid_molecules(atoms, scale, center=None):
    atoms = atoms.copy()
    atoms.center()
    if center is None:
        center = atoms.get_center_of_mass()
    center = np.asarray(center, dtype=float)

    # --- identify molecules via bond connectivity ---
    mol_ids = _connected_components(atoms)

    new_positions = atoms.get_positions().copy()
    masses = atoms.get_masses()

    for ids in mol_ids:
        idx = np.array(ids)
        m = masses[idx]
        com = (new_positions[idx] * m[:, None]).sum(0) / m.sum()
        # rigid-body displacement of the molecular COM relative to center
        new_com = center + scale * (com - center)
        shift = new_com - com
        new_positions[idx] += shift  # same shift for all atoms => bonds preserved

    out = atoms.copy()
    out.set_positions(new_positions)

    # scale the cell about the chosen center, keeping fractional layout sane
    cell = atoms.get_cell()
    out.set_cell(cell * scale, scale_atoms=False)
    out.center()
    return out


def _connected_components(atoms, mult=1.25):
    """Group atoms into molecules using covalent-radius bond cutoffs (PBC-aware)."""
    from ase.data import covalent_radii

    n = len(atoms)
    pos = atoms.get_positions()
    cell = atoms.get_cell()
    pbc = atoms.get_pbc()
    radii = covalent_radii[atoms.get_atomic_numbers()]

    # all pairwise minimum-image distances
    _, D = get_distances(pos, pos, cell=cell, pbc=pbc)

    adj = [[] for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            cutoff = mult * (radii[i] + radii[j])
            if D[i, j] < cutoff:
                adj[i].append(j)
                adj[j].append(i)

    seen = [False] * n
    comps = []
    for i in range(n):
        if seen[i]:
            continue
        stack, comp = [i], []
        seen[i] = True
        while stack:
            a = stack.pop()
            comp.append(a)
            for b in adj[a]:
                if not seen[b]:
                    seen[b] = True
                    stack.append(b)
        comps.append(comp)
    return comps

def add_molecule(host, mol, box, min_dist=2.5, max_tries=500, rng=None):
    rng = rng or np.random.default_rng()
    mol = mol.copy()
    mol.translate(-mol.get_center_of_mass())
    for _ in range(max_tries):
        trial = mol.copy()
        # random orientation
        trial.rotate(rng.uniform(0, 360), 'x', center=(0, 0, 0))
        trial.rotate(rng.uniform(0, 360), 'y', center=(0, 0, 0))
        trial.rotate(rng.uniform(0, 360), 'z', center=(0, 0, 0))
        # random position with margin from walls
        trial.translate(rng.uniform(2.0, box - 2.0, size=3))
        if len(host) == 0:
            return host + trial
        combined = host + trial
        d = combined.get_all_distances(mic=True)
        # distances between new atoms and existing atoms
        sub = d[len(host):, :len(host)]
        if sub.min() > min_dist:
            return combined
    raise RuntimeError("Could not place molecule without overlap")