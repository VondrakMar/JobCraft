"""
Pull some per-frame scalars into atoms.info and a custom per-atom array,
keep the weight_vector -> forces mapping on (it's on by default):
    python h5md_to_xyz.py --input_file traj.h5 --output_file traj.xyz \\
        --info energy=data/potential_energy \\
        --info temperature=data/temperature \\
        --info kinetic_energy=data/kinetic_energy \\
        --array custom_field=data/some_per_atom_dataset

Turn the weight_vector -> forces mapping off:
    python h5md_to_xyz.py --input_file traj.h5 --output_file traj.xyz \\
        --no-weights-as-forces

Unwrap positions so atoms don't jump across periodic boundaries:
    python h5md_to_xyz.py --input_file traj.h5 --output_file traj.xyz \\
        --unwrap-positions

Or call it directly from Python:
    from h5md_to_xyz import convert_h5_xyz
    convert_h5_xyz(
        "traj.h5", "traj.xyz",
        info_fields={"energy": "data/potential_energy"},
        array_fields={},
        weights_as_forces=True,
    )
"""

import sys
import argparse
from typing import Optional

import h5py
import numpy as np
from ase import Atoms
from ase.io import write

def unwrap_positions(positions: np.ndarray, cell: np.ndarray) -> np.ndarray:
    '''
    stolen from torch_SSPD code, since I dont want that as a dependency
    https://gitlab.tuwien.ac.at/e165-03-1_theoretische_materialchemie/madsen-s-research-group/StochasticSaddlePointDynamics
    '''

    positions = np.asarray(positions, dtype=float)
    cell = np.asarray(cell, dtype=float)
    frac = np.matmul(positions, np.linalg.inv(cell))
    delta = frac[1:] - frac[:-1]
    delta -= np.round(delta)
    unwrapped_frac = np.concatenate(
        [frac[:1], frac[:1] + np.cumsum(delta, axis=0)], axis=0
    )
    return np.matmul(unwrapped_frac, cell)


def _load_step_map(f: h5py.File, h5_path: str, steps_path: Optional[str] = None):
    if h5_path not in f:
        print(f"warning: {h5_path!r} not found in file, skipping", file=sys.stderr)
        return None
    
def _load_step_map_arrays(f: h5py.File, h5_path: str, steps_path: Optional[str] = None,n_atoms = None):
    if h5_path not in f:
        print(f"warning: {h5_path!r} not found in file, skipping", file=sys.stderr)
        return None

    name = h5_path.rsplit("/", 1)[-1]
    steps_path = steps_path or f"steps/{name}"
    if steps_path not in f:
        print(f"warning: {steps_path!r} not found, skipping {h5_path!r}", file=sys.stderr)
        return None

    values = f[h5_path][:]
    if n_atoms is not None:
        values = values.reshape(values.shape[0],n_atoms,(int(values.shape[-1]/n_atoms)))
    steps = f[steps_path][:]
    return dict(zip(steps.tolist(), values))


def load_h5md(
    path: str,
    info_fields: Optional[dict[str, str]] = None,
    array_fields: Optional[dict[str, str]] = None,
    weights_as_forces: bool = True,
    weight_vector_path: str = "data/weight_vector",
    unwrap: bool = False,
) -> list[Atoms]:
    info_fields = info_fields or {}
    array_fields = array_fields or {}

    with h5py.File(path, "r") as f:
        atomic_numbers = f["data/atomic_numbers"][0]          # (N,)
        positions_all  = f["data/positions"][:]               # (n_pos, N, 3)
        cell_all       = f["data/cell"][:]                    # (n_pos, 3, 3)
        pbc            = f["data/pbc"][:].astype(bool)        # (3,)
        pos_steps      = f["steps/positions"][:]              # (n_pos,)
        n_atoms = atomic_numbers.shape[0]

        if unwrap:
            positions_all = unwrap_positions(positions_all, cell_all)

        info_maps = {}
        for key, h5_path in info_fields.items():
            m = _load_step_map(f, h5_path)
            if m is not None:
                info_maps[key] = m

        array_maps = {}
        for key, h5_path in array_fields.items():
            m = _load_step_map_arrays(f, h5_path,n_atoms=atomic_numbers.shape[0])
            if m is not None:
                array_maps[key] = m

        weight_map = None
        if weights_as_forces and weight_vector_path in f:
            weight_vector_all = f[weight_vector_path][:]      # (n_w, 1, N*3)
            weight_vector_all = weight_vector_all.reshape(
                weight_vector_all.shape[0], n_atoms, 3
            )
            weight_name = weight_vector_path.rsplit("/", 1)[-1]
            weight_steps = f[f"steps/{weight_name}"][:]       # (n_w,)
            weight_map = dict(zip(weight_steps.tolist(), weight_vector_all))
        elif weights_as_forces:
            print(
                f"warning: weights_as_forces=True but {weight_vector_path!r} "
                "not found, skipping",
                file=sys.stderr,
            )

    frames = []
    for i, step in enumerate(pos_steps):
        atoms = Atoms(
            numbers=atomic_numbers,
            positions=positions_all[i],
            cell=cell_all[i],
            pbc=pbc,
        )

        for key, value_map in info_maps.items():
            value = value_map.get(int(step))
            if value is not None:
                arr = np.asarray(value)
                atoms.info[key] = arr.item() if arr.size == 1 else arr

        for key, value_map in array_maps.items():
            value = value_map.get(int(step))
            if value is not None:
                atoms.arrays[key] = np.asarray(value)

        if weight_map is not None:
            weights = weight_map.get(int(step))
            if weights is not None:
                norm = np.linalg.norm(weights, axis=1).max()
                if norm > 0:
                    weights = weights / norm
                atoms.arrays["forces"] = weights

        frames.append(atoms)

    return frames


def convert_h5_xyz(
    h5md_path: str,
    out_path: Optional[str] = None,
    info_fields: Optional[dict[str, str]] = None,
    array_fields: Optional[dict[str, str]] = None,
    weights_as_forces: bool = True,
    unwrap: bool = False,
) -> str:
    frames = load_h5md(
        h5md_path,
        info_fields=info_fields,
        array_fields=array_fields,
        weights_as_forces=weights_as_forces,
        unwrap=unwrap,
    )
    write(out_path, frames, format="extxyz")
    return out_path


def _parse_field_args(pairs: list[str]) -> dict[str, str]:
    """Turn ['key=data/path', ...] (as given repeatedly on the CLI) into a dict."""
    result: dict[str, str] = {}
    for pair in pairs:
        if "=" not in pair:
            raise argparse.ArgumentTypeError(f"expected KEY=H5_PATH, got {pair!r}")
        key, h5_path = pair.split("=", 1)
        result[key] = h5_path
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--input_file", help="Where is the h5 file", type=str, required=True)
    parser.add_argument("--output_file", help="Name of the xyz file", type=str, required=True)
    parser.add_argument(
        "--info",
        action="append",
        default=[],
        metavar="KEY=H5_PATH",
        help="Extra per-frame value to store in atoms.info. Repeatable. "
             "Example: --info energy=data/potential_energy",
    )
    parser.add_argument(
        "--array",
        action="append",
        default=[],
        metavar="KEY=H5_PATH",
        help="Extra per-atom array (already shaped n_steps x n_atoms x ...) "
             "to store in atoms.arrays. Repeatable. "
             "Example: --array custom_field=data/some_dataset",
    )
    parser.add_argument(
        "--weights-as-forces",
        dest="weights_as_forces",
        action="store_true",
        default=True,
        help="Load data/weight_vector, reshape + normalize it, and store as "
             "atoms.arrays['forces'] (default: on)",
    )
    parser.add_argument(
        "--no-weights-as-forces",
        dest="weights_as_forces",
        action="store_false",
        help="Disable the weight_vector -> forces mapping",
    )
    parser.add_argument(
        "--unwrap-positions",
        dest="unwrap",
        action="store_true",
        default=False,
        help="Unwrap positions across frames so atoms don't jump when they "
             "cross a periodic boundary (off by default -- positions are "
             "written as stored)",
    )
    args = parser.parse_args()

    convert_h5_xyz(
        h5md_path=args.input_file,
        out_path=args.output_file,
        info_fields=_parse_field_args(args.info),
        array_fields=_parse_field_args(args.array),
        weights_as_forces=args.weights_as_forces,
        unwrap=args.unwrap,
    )



'''
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
        # energy_all     = f["data/potential_energy"][:, 0]     # (n_e,)
        # temp_all       = f["data/temperature"][:, 0]          # (n_e,)
        pos_steps      = f["steps/positions"][:]              # (n_pos,)
        # print(positions_all)
        n_atoms = atomic_numbers.shape[0]

        weight_vector_all = None
        weight_map = None
        if "weight_vector" in f["data"]:
            # saved as (n_steps, 1, n_atoms * 3) -> reshape per-frame to (n_atoms, 3)
            weight_vector_all = f["data/weight_vector"][:]        # (n_w, 1, N*3)
            weight_vector_all = weight_vector_all.reshape(
                weight_vector_all.shape[0], n_atoms, 3
            )
            weight_steps = f["steps/weight_vector"][:]            # (n_w,)
            weight_map = dict(zip(weight_steps.tolist(), weight_vector_all))

    frames = []
    for i, step in enumerate(pos_steps):
        atoms = Atoms(
            numbers=atomic_numbers,
            positions=positions_all[i],
            cell=cell_all[i],
            pbc=pbc,
        )
        # atoms.info["energy"]         = float(energy_all[i])
        # atoms.info["kinetic_energy"] = float(ke_all[i])
        # atoms.info["temperature"]    = float(temp_all[i])
        if weight_map is not None:
            weights = weight_map.get(int(step))
            if weights is not None:
                norm = np.linalg.norm(weights, axis=1).max()
                if norm > 0:
                    weights = weights / norm
                # store as "forces" so ase gui / extxyz round-trip it as a
                # per-atom vector field and lets you view it with the
                # built-in force-vector visualization
                atoms.arrays["forces"] = weights
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
'''
