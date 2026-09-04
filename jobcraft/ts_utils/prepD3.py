import os
from pathlib import Path

import torch

from D3utils import (
    extract_dftd3_parameters,
    save_dftd3_parameters,
)

from nvalchemiops.torch.neighbors import neighbor_list

param_file = (
    Path("dftd3_parameters.pt")
)

print(Path.cwd())
if not param_file.exists():
    print("Downloading DFT-D3 parameters...")
    params = extract_dftd3_parameters()
    save_dftd3_parameters(params)
else:
    params = torch.load(param_file, weights_only=True)
    print("Loaded cached DFT-D3 parameters")
