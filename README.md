# Job Operations and Batch Creation, Resource Allocation, File Templating
JobCraft is simple software allowing to create submit files and jobs for my PhD.

It is now usable only for me (and people at the FHI Berlin, since all the settings are for our HPCs), so you have to modify `my_presets.py` to let if work for your system also.

Final submission scripts are using `slurm`, and GNU `parllalel`.

Parts of the code are also scripts for working with output files of softwares I am using. This is what ASE is doing, but I tried to change it, however it was causing more problems than I was willing to try to solve, so it was faster to write my simple code.

Yes, the name was created with chatGPT, but I really like starcraft, so deal with it

## Examples of usage, these are now random just for me
```Bash
python /home/mvondrak/work/JobCraft/jobcraft/preparation_script.py raven --usedN 4 -N 1 -n 4 --prep_submit --strucs waters.xyz --method="aims" --aims_basis="tight" --wall_time="1:00:00" --aims_species_path="/home/mvondrak/software/fhi-aims.231212_1/species_defaults/defaults_2020/" --aims_geometry_lines "homogeneous_field 0.25 0 0" "another test" "sea creatures are living on the ground"
```