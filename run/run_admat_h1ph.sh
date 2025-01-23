echo "Calculation of 2-phonon AD matrices"
#cat $PBS_NODEFILE > nodes.txt
mpirun -np $NUM_MPI_PROCS ./eqm_hole_1phon_admat > log_ad 2>error_ad

#/software/openmpi-1.8.2/intel/bin/mpirun -np $NUM_MPI_PROCS ./eqm_even_admat > log_ad 2>error_ad
