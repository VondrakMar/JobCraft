#!/bin/bash
rm *resubmit*
script_path=/home/mvondrak/development/JobCraft/jobcraft/restart_script.py 
python ${script_path} raven \
    --usedN 4 \
    -N 2 \
    --prep_submit \
    --strucs test_file.xyz \
    --method="aims" \
    --aims_basis="tight" \
    --per_submit 2 \
    --aims_kwargs '{"xc": "hse_06 0.11", "hse_unit":"bohr-1", "relativistic": ["atomic_zora", "scalar"], "compute_forces":true, "spin":"collinear","charge":1,"fixed_spin_moment":1}' \
    --wall_time="3:00:00" \
    --aims_species_path="/home/mvondrak/software/fhi-aims.231212_1/species_defaults/defaults_2020/" \
    
    # "initial_momen 0.1" "initial_momen 2.5" "initial_momen 0.3" \
    # --aims_kwargs '{"xc": "hse_06 0.11", "hse_unit":"bohr-1", "relativistic": ["atomic_zora", "scalar"], "compute_forces":true, "k_grid":[3,3,3], "spin":"collinear","charge":1,"fixed_spin_moment":1}' \ 
