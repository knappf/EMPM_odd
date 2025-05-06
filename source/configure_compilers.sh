export XLF="ifx -O3" #  -qopenmp" 
export XLFMPI="mpiifx -O3 -mcmodel=large"  # this should work on Metacentrum with appropiate modules  
#export XLFMPI="mpiifx -O3 -mcmodel=large"  # newer version of Intel Fotran compiller, use @ipnp36 
export mpirun="/afs/ics.muni.cz/software/intel-oneapi/intel-oneapi-hpc-toolkit-2025.1.0/oneapi/mpi/2021.15/bin/mpirun"  # alias for Metacentrum with complete path 
export EMPM_RUN_DIR="../run/"
export EMPM_BIN_DIR="../bin/"
mkdir -p ../bin/


