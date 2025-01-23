export XLF="ifx -O3" #  -qopenmp" 
export XLFMPI="mpiifort -O3 -mcmodel=large"  # this should work on Metacentrum with appropiate modules  
#export XLFMPI="mpiifx -O3 -mcmodel=large"  # newer version of Intel Fotran compiller, use @ipnp36 
export mpirun="/opt/intel/oneapi/mpi/2021.14/bin/mpirun"  # alias for Metacentrum with complete path 
export EMPM_RUN_DIR="../run/"
export EMPM_BIN_DIR="../bin/"
mkdir -p ../bin/


