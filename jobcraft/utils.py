# this is just garbage of stuff I wrote and did not want to loose them
def prep_restart_DFTs(per_file = 4,at_the_same_time=2):
    import subprocess
    command = "for a in */; do echo $a >> temp; tail -2 ${a}aims.out | head -1 >> temp ; done"
    subprocess.run(command, shell=True, capture_output=True, text=True) # this is way faster than open each file in python
    lines = open("temp","r").readlines()
    prev_line = lines[0]
    to_restart = []
    for line in range(1,len(lines)):
        if "struc" in lines[line] and "struc" in prev_line:
            to_restart.append(prev_line) 
        if "struc" not in lines[line] and ("Have a nice day" not in lines[line] or "scf_solver: SCF cycle not converged" not in lines[line]) and "struc" in prev_line:
            to_restart.append(prev_line)
        prev_line = lines[line]
    #paral_line = "srun -N 4 -n 288 /u/mvondrak/software/fhi-aims.231212_1/build/aims.231212_1.scalapack.mpi.x >> aims.out"
    temp = open("submit_file1.sl").readlines()

    # parallel --delay 0.2 --joblog task.log --progress -j 2 < paral_file1
    head_file = ["#!/bin/bash -l\n"]
    for line in temp:
        if "SBATCH" in line or "module" in line or "export" in line:
            head_file.append(line)

    paral_line = open("paral_file1").readlines()[0].split(";")[1]
    count = -1
    
    with open("newCMD","w") as newCMD:
        for id_line,restart_line in enumerate(to_restart):
            if id_line%per_file == 0:
                if id_line != 0:
                    slurm_file.close()
                    paral_file.close()
                count+=1
                slurm_file = open(f"submit_restart{count}.sl","w")
                for head_line in head_file:
                    slurm_file.write(head_line)
                slurm_file.write(f"parallel --delay 0.2 --joblog task.log --progress -j {at_the_same_time} < restart_file{count}")
                paral_file = open(f"restart_file{count}","w")
            current_com = f"cd {restart_line[:-1]} ; {paral_line}\n"
            paral_file.write(f"{current_com}")#f"cd {restart_line}; srun -N {self.node_per_job} -n {self.cpu_per_job_to_srun} {aims_command} > aims.out")


            #newCMD.write(current_com)
# prep_restart_DFTs(per_file=2,at_the_same_time=1)



def add_initial_moment_water():
    from ase.io import read
    import sys
    from matscipy.neighbours import neighbour_list
    import shutil
    mol = read(f"geometry.in",format="aims")
    mol.center(vacuum=10)
    i,j = neighbour_list("ij",mol,cutoff=1.2)



    elements = mol.get_chemical_symbols()
    neigh1 = i
    neigh2 = j

    # Find indices of oxygen atoms
    oxygen_indices = [i for i, elem in enumerate(elements) if elem == 'O']


    neigh1_count = {}

    for number in neigh1:
        if number in neigh1_count:
            neigh1_count[number] += 1
        else:
            neigh1_count[number] = 1

    for i in neigh1_count:
        if i in oxygen_indices:
            if neigh1_count[i] != 2:
                single_n = i
    shutil.move("geometry.in","temp.in")
    f = open("temp.in","r").readlines()
    geometry = open("geometry.in","w")
    for id_line in range(0,5+single_n):
        geometry.write(f[id_line])
    geometry.write("initial_moment 1\n")
    for id_line in range(5+single_n,len(f)):
        geometry.write(f[id_line])
