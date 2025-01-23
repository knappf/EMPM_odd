# compiler 
#XLF= ifort -g #-O3 -parallel #-O3 -parallel
#XLFMPI= mpiifort -g  -mcmodel=large #-g -debug full #-O3 -parallel 
# options
#OPT= -L$MKLPATH -I$MKLINCLUDE -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread

EXE_PATH=../${EMPM_BIN_DIR}
# exe file
EXF= phon_dens2_MPI
#
all: anglib.o dens2_MPI.o input_sp.o phon_dens2_MPI.o
	$(XLFMPI) anglib.o dens2_MPI.o input_sp.o phon_dens2_MPI.o ${OPT} -o ${EXE_PATH}${EXF}

anglib.o: anglib.f90
	${XLFMPI} -c anglib.f90

dens2_MPI.o: dens2_MPI.f90
	${XLFMPI} -c dens2_MPI.f90

input_sp.o: input_sp.f90
	${XLFMPI} -c input_sp.f90

phon_dens2_MPI.o: phon_dens2_MPI.f90
	${XLFMPI} -c phon_dens2_MPI.f90

clean:
	rm -f *.o *.mod
	rm -f ${EXF}$     

