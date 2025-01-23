# compiler 
#XLF= ifort -O3 -parallel
XLF= ifort -g
# options
OPT= -L$MKLPATH -I$MKLINCLUDE -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread

EXE_PATH=../${EMPM_BIN_DIR}
# exe file
EXF= eqm_part_phon
#
all: input_sp.o anglib.o base_p_phon.o read_dens_and_int.o dmat.o eqm_part.o
	$(XLF) input_sp.o anglib.o base_p_phon.o read_dens_and_int.o dmat.o eqm_part.o ${OPT} -o ${EXE_PATH}${EXF}

input_sp.o: input_sp.f90
	${XLF} -c input_sp.f90

anglib.o: anglib.f90
	${XLF} -c anglib.f90

base_p_phon.o: base_p_phon.f90
	${XLF} -c base_p_phon.f90

read_dens_and_int.o: read_dens_and_int.f90
	${XLF} -c read_dens_and_int.f90

dmat.o: dmat.f90
	${XLF} -c dmat.f90

eqm_part.o: eqm_part.f90 
	${XLF} -c eqm_part.f90

clean:
	rm -f *.o *.mod
	rm -f ${EXF}




	
