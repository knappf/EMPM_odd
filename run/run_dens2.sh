echo "Calculation of 2-phonon densities"
#cat $PBS_NODEFILE > nodes.txt
#source mpi_openmp_set.sh
/opt/intel/oneapi/mpi/2021.8.0/bin/mpirun -np $NUM_MPI_PROCS ./phon_dens2_MPI > log_dens2 2>error_dens2
#/software/openmpi-1.8.2/intel/bin/mpirun -np $NUM_MPI_PROCS ./phon_dens2_MPI > log_dens2 2>error_dens2
