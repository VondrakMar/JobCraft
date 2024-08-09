# Job Operations and Batch Creation, Resource Allocation, File Templating
JobCraft is simple software allowing to create submit files and jobs for my PhD.

It is now usable only for me (and people at the FHI Berlin, since all the settings are for our HPCs), so you have to modify `my_presets.py` to let if work for your system also.

Final submission scripts are using `slurm`, and GNU `parllalel`.

Parts of the code are also scripts for working with output files of softwares I am using. This is what ASE is doing, but I tried to change it, however it was causing more problems than I was willing to try to solve, so it was faster to write my simple code.

Yes, the name was created with chatGPT, but I really like starcraft, so deal with it
