#!/bin/bash
script_path=/home/mvondrak/development/JobCraft/jobcraft/preparation_script.py 
python ${script_path} raven \
    --usedN 4 \
    -N 1 \
    -n 16 \
    --prep_submit \
    --strucs bulk_sample.xyz \
    --method="aims" \
    --aims_basis="tight" \
    --per_submit 2 \
    --aims_kwargs '{"xc": "hse_06 0.11", "k_grid_density":3,"hse_unit":"bohr-1", "relativistic": ["atomic_zora", "scalar"], "compute_forces":true, "spin":"collinear","charge":1,"fixed_spin_moment":1}' \
    --aims_outputs hirshfeld hartree_multipoles dipole \
    --wall_time="1:00:00" \
    --aims_species_path="/home/mvondrak/software/fhi-aims.231212_1/species_defaults/defaults_2020/" \
    --aims_geometry_lines "homogeneous_field 0.25 0 0" "another test" "sea creatures are living on the ground" \
    --aims_atoms_lines "initial_moment 0.5" \
    --aims_atoms_indeces 1 7 17
    
    # "initial_momen 0.1" "initial_momen 2.5" "initial_momen 0.3" \
    # --aims_kwargs '{"xc": "hse_06 0.11", "hse_unit":"bohr-1", "relativistic": ["atomic_zora", "scalar"], "compute_forces":true, "k_grid":[3,3,3], "spin":"collinear","charge":1,"fixed_spin_moment":1}' \ 
rm header_file.temp control.in temp.in
     
