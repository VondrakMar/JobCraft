# run_nvt_parallel.py
import torch
import torch_sim as ts
import os
from ase.build import bulk
from ase.io import read,write
from mace.calculators.foundations_models import mace_mp
from torch_sim.models.mace import MaceModel
from pathlib import Path
from torch_sim.models.dispersion import D3DispersionModel
from torch_sim.models.interface import SumModel
import time
from nvalchemiops.torch.interactions.dispersion import (
    D3Parameters,
    dftd3,
)

def load_d3_parameters(
    param_file: Path | None = None) -> D3Parameters:
    if param_file is None:
        cache_dir = Path(os.path.expanduser("~")) / ".cache" / "nvalchemiops"
        param_file = cache_dir / "dftd3_parameters.pt"

    if not param_file.exists():
        raise FileNotFoundError(
            f"DFT-D3 parameter file not found: {param_file}\n"
            "Please run one of the example scripts to generate the parameter file.\n"
            "Example: python examples/dispersion/01_dftd3_molecule.py"
        )

    state_dict = torch.load(param_file, map_location="cpu", weights_only=True)

    return D3Parameters(
        rcov=state_dict["rcov"],
        r4r2=state_dict["r4r2"],
        c6ab=state_dict["c6ab"],
        cn_ref=state_dict["cn_ref"],
    )


def runMDmaceD3():
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    mace_raw = mace_mp(model="mace-mh-1.model",device="cuda",head="omat_r2scan",return_raw_model=True)
    mace_model = MaceModel(model=mace_raw, device=device)

    base_atoms = read("struc.xyz",format="extxyz")

    n_replicas = 2
    systems = [base_atoms.copy() for _ in range(n_replicas)]

    traj_files = [f"nvt_CuO_replica_{i:04d}.h5md" for i in range(n_replicas)]
    params = load_d3_parameters()
    d3 = D3DispersionModel(
        a1=0.49484001,
        s8=0.78981345,
        a2=5.73083694,
        s6=1.0,
        cutoff=20.0,
        d3_params = params,
        dtype=mace_model.dtype,
        device=mace_model.device,
    )

    model = SumModel(mace_model, d3)
    t0 = time.time() 
    relaxed_state = ts.optimize(
        system=systems,
        model=model,
        optimizer=ts.Optimizer.fire,
        autobatcher=False,
        max_steps=50,
        init_kwargs=dict(cell_filter=ts.CellFilter.frechet),
    )
    t1 = time.time()
    final_state = ts.integrate(
        system=relaxed_state,
        model=model,
        n_steps=100,
        timestep=0.001,          # ps; 1 fs
        temperature=300,         # K
        integrator=ts.Integrator.nvt_langevin,
        trajectory_reporter=dict(
            filenames=traj_files,
            prop_calculators={
                10:{
                    "potential_energy":lambda state:state.energy,
                    "forces": lambda state:state.forces,
                    "temperature":lambda state:state.calc_temperature()
                }
            },
            state_frequency=10,
        ),
    )
    final_atoms = final_state.to_atoms()


if __name__ == "__main__":
    runMDmaceD3()
