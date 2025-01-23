# compiler 
XLF= ifort -O3 -parallel
#LF= ifort -g
# options
OPT= -L$MKLPATH -I$MKLINCLUDE -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread
#OPT= -mkl=sequential
EXE_PATH=../${EMPM_BIN_DIR}
# exe file
EXF= eqm_hole_phon_svd
#
all: input_sp.o anglib.o base_h_phon.o read_dens_and_int.o choles.o dmat.o eqm_hole.o
	$(XLF) input_sp.o anglib.o base_h_phon.o read_dens_and_int.o choles.o dmat.o eqm_hole.o ${OPT} -o ${EXE_PATH}${EXF}

input_sp.o: input_sp.f90
	${XLF} -c input_sp.f90

anglib.o: anglib.f90
	${XLF} -c anglib.f90

base_h_phon.o: base_h_phon.f90
	${XLF} -c base_h_phon.f90

read_dens_and_int.o: read_dens_and_int.f90
	${XLF} -c read_dens_and_int.f90

choles.o: choles.f
	${XLF} -c choles.f

dmat.o: dmat.f90
	${XLF} -c dmat.f90

eqm_hole.o: eqm_hole.f90
	${XLF} -c eqm_hole.f90

clean:
	rm -f *.o *.mod
	rm -f ${EXF}




	
