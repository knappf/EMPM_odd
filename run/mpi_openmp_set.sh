#export MPIRUN_PATH="/software/openmpi-1.8.2/intel/bin/"      # Metacentrum urga, ursa, uruk
#export MPIRUN_PATH="/opt/intel/compilers_and_libraries_2020.1.217/linux/mpi/intel64/bin/"    # ipnp16
# number of MPI procs
#alias mpirun="/software/openmpi-1.8.2/intel/bin/mpirun"
export NUM_MPI_PROCS=16
echo "# MPI procs"
echo $NUM_MPI_PROCS 
# number of OMP threads 
export OMP_NUM_THREADS=16
echo "# OMP threads"
echo $OMP_NUM_THREADS
