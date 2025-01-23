# compiler 
XLF= ifort -O3 -parallel 
# options
#OPT= -L$MKLPATH -I$MKLINCLUDE -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread

EXE_PATH=../${EMPM_BIN_DIR}
# exe file
EXF= phon_int
#
all: input_sp.o anglib.o read_dens_and_int.o interac_iso.o phon_int.o
	$(XLF) input_sp.o anglib.o read_dens_and_int.o interac_iso.o phon_int.o ${OPT} -o ${EXE_PATH}${EXF}

input_sp.o: input_sp.f90
	${XLF} -c input_sp.f90

anglib.o: anglib.f90
	${XLF} -c anglib.f90

read_dens_and_int.o: read_dens_and_int.f90
	${XLF} -c read_dens_and_int.f90

interac_iso.o: interac_iso.f90 
	${XLF} -c interac_iso.f90

phon_int.o: phon_int.f90
	${XLF}$ -c phon_int.f90

clean:
	rm -f *.o *.mod
	rm -f ${EXF}$     



