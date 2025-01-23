# compiler 
XLF= ifort -O3 -parallel 
# options
#OPT= -L$MKLPATH -I$MKLINCLUDE -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread
# exe file
EXF= phon_dens1
#
EXE_PATH=../${EMPM_BIN_DIR}
all: input_sp.o anglib.o dens_pn.o phon_dens1.o
	$(XLF) input_sp.o anglib.o dens_pn.o phon_dens1.o -o ${EXE_PATH}${EXF}

input_sp.o: input_sp.f90
	${XLF} -c input_sp.f90

anglib.o: anglib.f90
	${XLF} -c anglib.f90

dens_pn.o: dens_pn.f90
	${XLF} -c dens_pn.f90

phon_dens1.o: phon_dens1.f90
	${XLF} -c phon_dens1.f90

clean:
	rm -f *.o *.mod
	rm -f ${EXF}$     

