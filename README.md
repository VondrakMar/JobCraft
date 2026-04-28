# Job Operations and Batch Creation, Resource Allocation, File Templating
Originally this was only for creating simple DFT jobs and their processing. Now it is more, but I doubt it is usefull for anyone except me. 

If you want to use it, I advice you againt it, but feel free. There are 0 tests, so be aware

## Examples of generating aims files from xyz. This is prolly depricated, sorry
```Bash
python /home/mvondrak/work/JobCraft/jobcraft/preparation_script.py raven --usedN 4 -N 1 -n 4 --prep_submit --strucs waters.xyz --method="aims" --aims_basis="tight" --wall_time="1:00:00" --aims_species_path="/home/mvondrak/software/fhi-aims.231212_1/species_defaults/defaults_2020/" --aims_geometry_lines "homogeneous_field 0.25 0 0" "another test" "sea creatures are living on the ground"
```
