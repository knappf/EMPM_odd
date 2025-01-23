echo "Calculation of 1-phonon densities"
#cat $PBS_NODEFILE > nodes.txt
mpirun -np $NUM_MPI_PROCS ./phon_dens1_MPI > log_phon_dens1 2>error_phon_dens1

#/software/openmpi-1.8.2/intel/bin/mpirun -np $NUM_MPI_PROCS ./eqm_even_admat > log_ad 2>error_ad
