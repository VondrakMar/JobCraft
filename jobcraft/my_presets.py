#aims for raven
aims_for_raven = 'export OMP_NUM_THREADS=1\nexport MKL_DYNAMIC=FALSE\nexport MKL_NUM_THREADS=1\nmodule load parallel/201807\nmodule load anaconda/3/2020.02\nmodule load intel/21.4.0\nmodule load mkl/2021.4\nmodule load impi/2021.4\nexport LD_LIBRARY_PATH="$LD_LIBRARY_PATH:${INTEL_HOME}/compiler/2021.4.0/linux/compiler/lib/intel64_lin/"\n'
aims_for_viper = 'export OMP_NUM_THREADS=1\nulimit -s unlimited\nmodule load anaconda/3/2023.03\nmodule load gcc/13 openmpi/4.1 mkl/2024.0\nmodule load parallel/201807\n'
aims_path_raven = '/u/mvondrak/software/fhi-aims.231212_1/'
aims_exec_raven="build/aims.231212_1.scalapack.mpi.x"
aims_species_raven='species_defaults/defaults_2020/'

aims_path_viper = "trololo"

