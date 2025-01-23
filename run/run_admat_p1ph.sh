echo "Calculation of 2-phonon AD matrices"
#cat $PBS_NODEFILE > nodes.txt
/opt/intel/oneapi/mpi/2021.10.0/bin/mpirun -np $NUM_MPI_PROCS ./eqm_part_1phon_admat > log_ad 2>error_ad

#/software/openmpi-1.8.2/intel/bin/mpirun -np $NUM_MPI_PROCS ./eqm_even_admat > log_ad 2>error_ad
