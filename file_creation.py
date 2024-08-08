def prep_aims_ase_file(base_file="aims_ase_template.py",final_name="runAims.py",aims_command="raven",aims_species="raven/species/"):
    with open(base_file, 'r') as base_file:
        base_lines = base_file.read()
    with open(f'{final_name}', 'w') as final_file:
        final_file.write(f'aims_command="{aims_command}"\n')
        final_file.write(f'aims_species="{aims_species}"\n')
        final_file.write("\n"+base_lines)

# prep_aims_ase_file()
def prep_submit_header(submitted_nodes,NTASK_PER_NODE,PRESETS_FOR_HEADER,wall_time="1:00:00"):
    header_file = open("header_file.temp","w")
    header_file.write("#!/bin/bash -l\n")
    header_file.write("#SBATCH -o ./tjob.out.%j\n") # stardart output file
    header_file.write("#SBATCH -e ./tjob.err.%j\n") # error output file
    header_file.write("#SBATCH -D ./\n") # working directory 
    # header_file.write(f"#SBATCH -J {job_name}\n") # name of the file
    header_file.write("\n")
    header_file.write(f"#SBATCH --nodes={submitted_nodes}\n")
    header_file.write(f"#SBATCH --ntasks-per-node={NTASK_PER_NODE}\n")
    header_file.write(f"#SBATCH --ntasks-per-node=1\n")
    header_file.write(f"#SBATCH --time={wall_time}\n")
    header_file.write("\n")
    header_file.write(f"{PRESETS_FOR_HEADER}")
    #header_file.write("#SBATCH --cpus-per-task=1\n")
    header_file.close()