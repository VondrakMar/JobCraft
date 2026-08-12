import ase
from ase import Atoms
from ase.io import read, write
from jobcraft.vasp.vasp_kpoint import write_kpoints, k_point_adjuster
import numpy as np
from pathlib import Path
import shutil
from pymatgen.io.vasp.inputs import Incar
from copy import deepcopy

# ---------------------------------------------------------------------------
# User options: pick the functional and whether to add DFT-D3(BJ) dispersion
# ---------------------------------------------------------------------------
# FUNCTIONAL: "pbe", "rpbe" or "r2scan"
FUNCTIONAL = "rpbe"
# VDW_D3BJ: True -> add DFT-D3 with Becke-Johnson damping (IVDW = 12)
VDW_D3BJ = False#True


# parameters from https://raw.githubusercontent.com/dftd3/simple-dftd3/main/assets/parameters.toml
D3BJ_PARAMS = {
    # PBE: VASP built-in defaults (listed here for transparency / reproducibility)
    "pbe":    {"vdw_s6": 1.0000, "vdw_s8": 0.7875, "vdw_a1": 0.4289, "vdw_a2": 4.4407},
    # RPBE-D3(BJ)
    "rpbe":   {"vdw_s6": 1.0000, "vdw_s8": 0.8318, "vdw_a1": 0.1820, "vdw_a2": 4.0094},
    # r2SCAN-D3(BJ)
    "r2scan": {"vdw_s6": 1.0000, "vdw_s8": 0.6019, "vdw_a1": 0.5156, "vdw_a2": 5.7734},
}


def apply_functional(kwargs: dict, functional: str, vdw_d3bj: bool) -> dict:
    """
    Set INCAR tags for the requested exchange-correlation functional and,
    optionally, DFT-D3(BJ) dispersion.

    - "pbe":    GGA = PE
    - "rpbe":   GGA = RP
    - "r2scan": METAGGA = R2SCAN (GGA tag removed; needs LASPH = .TRUE.)
    - vdw_d3bj: IVDW = 12 (Grimme D3 with Becke-Johnson damping); the
                functional-specific damping parameters from D3BJ_PARAMS are
                written explicitly so r2SCAN does not silently use PBE values.
    """
    kw = deepcopy(kwargs)
    functional = functional.lower()

    # Clear any XC / vdW tags that might leak in from the base settings
    for key in ("gga", "metagga", "ivdw",
                "vdw_s6", "vdw_s8", "vdw_a1", "vdw_a2"):
        kw.pop(key, None)

    if functional == "pbe":
        kw["gga"] = "PE"
    elif functional == "rpbe":
        kw["gga"] = "RP"
    elif functional == "r2scan":
        kw["metagga"] = "R2SCAN"
        kw["lasph"] = True          # required for meta-GGA
        kw["algo"] = "All"          # robust for meta-GGA SCF
    else:
        raise ValueError(
            f"Unknown functional: {functional!r}. Use 'pbe', 'rpbe' or 'r2scan'."
        )

    if vdw_d3bj:
        kw["ivdw"] = 12             # DFT-D3 with Becke-Johnson damping
        kw.update(D3BJ_PARAMS[functional])

    return kw


def write_incar(kwargs, path="INCAR"):
    with open(path, "w") as f:
        for k, v in kwargs.items():
            if isinstance(v, bool):
                v = ".TRUE." if v else ".FALSE."
            elif isinstance(v, (list, tuple)):
                v = " ".join(map(str, v))
            f.write(f"{k.upper()} = {v}\n")


def create_POSCAR(atoms: ase.atoms.Atoms, file_name="POSCAR"):
    write(f"{file_name}", atoms, format="vasp", sort=True)


def check_element_lines(file_name="POSCAR"):
    wanted = [1, 6]
    max_num = max(wanted)
    lines = {}
    with open(f"{file_name}") as f:
        for i, line in enumerate(f, start=1):
            if i in wanted:
                lines[i] = line.split()
            if i >= max_num:
                break
    return lines[wanted[0]] == lines[wanted[1]], lines[wanted[0]]


def create_POTCAR(element_list, pot_path, dict_paws, potcar_name="POTCAR"):
    for el in element_list:
        paw_name = dict_paws.get(el, el)
        tmp_pot_file = f"{pot_path}{paw_name}/POTCAR"
        with open(tmp_pot_file, "r") as src, open(potcar_name, "a") as dst:
            dst.write(src.read())


def move_to_folder(file_list_move, file_list_copy=[], base_name="struc",
                   id_name=None, num_digits=5, verbose=False):
    if id_name is not None:
        folder_name = f"{base_name}{id_name:0{num_digits}d}"
    else:
        folder_name = base_name
    folder = Path(folder_name)
    folder.mkdir(exist_ok=True)
    for file_name in file_list_move:
        file_path = Path(file_name)
        if file_path.exists():
            shutil.move(str(file_path), str(folder / file_name))
            if verbose:
                print(f"Moved {file_name} -> {folder_name}/")
        else:
            print(f"{file_name} not found")
    for file_name in file_list_copy:
        file_path = Path(file_name)
        if file_path.exists():
            shutil.copy(str(file_path), str(folder / file_name))
            if verbose:
                print(f"Copied {file_name} -> {folder_name}/")
        else:
            print(f"{file_name} not found")

def compute_dipol(atoms: Atoms, idipol: int = 3, mass_weighted: bool = True):
    if mass_weighted:
        center_cart = atoms.get_center_of_mass()
    else:
        center_cart = atoms.get_positions().mean(axis=0)
    frac = np.linalg.solve(atoms.cell.array.T, center_cart)
    frac = frac % 1.0
    if idipol in (1, 2, 3):
        axis = idipol - 1
        dipol = np.array([0.5, 0.5, 0.5])
        dipol[axis] = frac[axis]
    else:
        dipol = frac
    return dipol.tolist()


def detect_config_type(atoms: Atoms,config_string: str) -> str:
    """
    Decide system type. Priority:
      1. Fallback: no cell or pbc all False -> molecule.
      2. Explicit atoms.info['config_type'] if present and valid.
      3. Otherwise raise (caller must set config_type for slab/bulk).
    """
    # Fallback detection for molecules with no cell/pbc
    cell_lengths = atoms.cell.lengths() if atoms.cell is not None else np.zeros(3)
    no_cell = np.allclose(cell_lengths, 0.0)
    pbc = np.asarray(atoms.pbc, dtype=bool)
    no_pbc = not pbc.any()
    if no_cell or no_pbc:
        return "molecule"

    valid = {"molecule", "slab", "bulk"}
    ct = atoms.info.get(config_string, None)
    if ct is not None:
        ct = ct.lower()
        if ct in valid:
            return ct

    raise ValueError(
        "Could not determine config_type. Set atoms.info['config_type'] "
        "to one of {'molecule', 'slab', 'bulk'}."
    )


def mol_in_the_box(atoms: Atoms, vacuum: float = 8.0, idipol: int = 3) -> Atoms:
    """
    Place a molecule inside a cubic box with `vacuum` Å padding on each side,
    center it, and turn PBC on (VASP requires a cell).
    Returns a new Atoms object.
    """
    new_atoms = atoms.copy()
    new_atoms.set_pbc([True, True, True])

    # Build a box large enough to hold the molecule + vacuum on each side
    positions = new_atoms.get_positions()
    extents = positions.max(axis=0) - positions.min(axis=0)
    box = extents + 2.0 * vacuum
    # Make it cubic-ish (use the largest dimension) to keep things simple
    L = float(np.max(box))
    new_atoms.set_cell([L, L, L])
    new_atoms.center()
    return new_atoms


def prepare_molecule_kwargs(atoms: Atoms, base_kwargs: dict) -> dict:
    """
    INCAR settings for an isolated molecule in a box:
      - Gaussian smearing, small sigma
      - Gamma-point only typically (handled by k_point_adjuster outside)
      - Dipole correction in all three directions (IDIPOL=4)
    """
    kw = deepcopy(base_kwargs)
    kw["ldipol"] = True
    kw["idipol"] = 4  # full 3D dipole correction for isolated molecules
    kw["dipol"] = compute_dipol(atoms, idipol=4)
    # Molecules are typically open-shell-safe with these
    kw["ismear"] = 0
    kw["sigma"] = 0.01
    return kw


def prepare_slab_kwargs(atoms: Atoms, base_kwargs: dict, idipol: int = 3) -> dict:
    """
    INCAR settings for a slab: dipole correction along the surface normal.
    """
    kw = deepcopy(base_kwargs)
    kw["ldipol"] = True
    kw["idipol"] = idipol
    kw["dipol"] = compute_dipol(atoms, idipol=idipol)
    return kw


def prepare_bulk_kwargs(atoms: Atoms, base_kwargs: dict) -> dict:
    """
    INCAR settings for bulk: no dipole correction, Methfessel-Paxton is
    fine for metals but we keep Gaussian for safety. Tune as needed.
    """
    kw = deepcopy(base_kwargs)
    # Make sure no dipole keys leak in
    for key in ("ldipol", "idipol", "dipol"):
        kw.pop(key, None)
    return kw


# ---------------------------------------------------------------------------
# Base INCAR (no XC/dipole keys by default — XC is set via apply_functional)
# ---------------------------------------------------------------------------
vasp_kwargs = {
    "system": "skibidi",
    "encut": 900,
    "prec": "Accurate",
    "ediff": 1e-6,
    "algo": "All",
    "nelm": 250,
    "lreal": "Auto",
    "ismear": 0,
    "sigma": 0.05,
    "ibrion": -1,
    "nsw": 0,
    "isif": 0,
    "lwave": False,
    "lcharg": False,
    "ncore": 10,
    "lasph": True,
    "kpar": 1,
}


# ---------------------------------------------------------------------------
# Main driver
# ---------------------------------------------------------------------------
import sys
if __name__ == "__main__":
    PBE_POTCHAR_PATH = "/home/mvondrak/software/vasp.6.4.2paws/PBE/"
    dict_paws = {
        "Cu": "Cu",
        "H": "H",
        "O": "O",
    }
    mol_name = str(sys.argv[1])

    # Optional CLI overrides:  python prepPOTCAR_functional.py <file> [functional] [vdw_d3bj]
    functional = FUNCTIONAL
    vdw_d3bj = VDW_D3BJ
    if len(sys.argv) > 2:
        functional = sys.argv[2]
    if len(sys.argv) > 3:
        vdw_d3bj = sys.argv[3].lower() in ("1", "true", "yes", "on")

    print(f"Functional: {functional} | DFT-D3(BJ): {vdw_d3bj}")

    #mols = read(f"{mol_name}@:", format="extxyz")
    mols = read(f"{mol_name}@:")
    vasp_kwargs["encut"] = 900
    for id_mol, mol in enumerate(mols):
        config_type = "slab"

        if config_type == "molecule":
            atoms_out = mol_in_the_box(mol, vacuum=8.0, idipol=3)
            cur_vasp_kwargs = prepare_molecule_kwargs(atoms_out, vasp_kwargs)
            kgrid = [1, 1, 1]  # Gamma-only for isolated molecules
        elif config_type == "slab":
            atoms_out = mol
            cur_vasp_kwargs = prepare_slab_kwargs(atoms_out, vasp_kwargs, idipol=3)
            kgrid = k_point_adjuster(atoms_out, kspacing=0.15, slab=True)
        elif config_type == "bulk":
            atoms_out = mol
            cur_vasp_kwargs = prepare_bulk_kwargs(atoms_out, vasp_kwargs)
            kgrid = k_point_adjuster(atoms_out, kspacing=0.15, slab=False)
        else:
            raise RuntimeError(f"Unhandled config_type: {config_type}")

        # Apply exchange-correlation functional + optional DFT-D3(BJ)
        cur_vasp_kwargs = apply_functional(cur_vasp_kwargs, functional, vdw_d3bj)

        cur_vasp_kwargs["system"] = f"{config_type}_{id_mol}"

        # Write inputs
        write_kpoints(kpoints=kgrid)
        create_POSCAR(atoms_out, file_name="POSCAR")
        write_incar(cur_vasp_kwargs)

        is_fine, el_list = check_element_lines()
        if is_fine:
            create_POTCAR(
                element_list=el_list,
                pot_path=PBE_POTCHAR_PATH,
                dict_paws=dict_paws,
            )
            move_to_folder(
                file_list_move=["POSCAR", "POTCAR", "KPOINTS", "INCAR"],
                file_list_copy=[],
                # base_name=f"{config_type}_struc",
                base_name=f"struc",
                id_name=id_mol,
                num_digits=5,
            )
        else:
            print(f"[{id_mol}] element line mismatch in POSCAR — skipping")
