# compiler 
#XLFMPI= ifort -O3 -parallel 
# options
#OPT= -L$MKLPATH -I$MKLINCLUDE -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread
# exe file
EXF= phon_dens1_MPI
#
EXE_PATH=../${EMPM_BIN_DIR}
all: input_sp.o anglib.o dens_MPI.o phon_dens1_MPI.o
	$(XLFMPI) input_sp.o anglib.o dens_MPI.o phon_dens1_MPI.o -o ${EXE_PATH}${EXF}

input_sp.o: input_sp.f90
	${XLFMPI} -c input_sp.f90

anglib.o: anglib.f90
	${XLFMPI} -c anglib.f90

dens_MPI.o: dens_MPI.f90
	${XLFMPI} -c dens_MPI.f90

phon_dens1_MPI.o: phon_dens1_MPI.f90
	${XLFMPI} -c phon_dens1_MPI.f90

clean:
	rm -f *.o *.mod
	rm -f ${EXF}$     

