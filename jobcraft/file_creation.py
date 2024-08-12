import os
import ase.io



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
    header_file.write(f"#SBATCH --cpus-per-task=1\n")
    header_file.write(f"#SBATCH --time={wall_time}\n")
    header_file.write("\n")
    header_file.write(f"{PRESETS_FOR_HEADER}")
    #header_file.write("#SBATCH --cpus-per-task=1\n")
    header_file.close()

def save_results_to_xyz(mol,results,prepositon="dft_",saved_name="struc_saved.xyz"):
    for res in results:
        if type(results[res]) == float:
            mol.info[f"{prepositon}{res}"] = results[res]
        elif len(results[res]) == len(mol):
            mol.arrays[f"{prepositon}{res}"] = results[res]
        else: 
            mol.info[f"{prepositon}{res}"] = results[res] # this is not tested, it should load dipole vector for example
    ase.io.write(saved_name,mol)