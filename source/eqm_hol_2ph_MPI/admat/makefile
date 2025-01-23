# compiler 
#XLF= ifort -O3 -parallel
#XLF= ifort -g
# options
OPT= -L$MKLPATH -I$MKLINCLUDE -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread
EXE_PATH=../../${EMPM_BIN_DIR}
# exe file
EXF= eqm_hole_2phon_admat
#
all: input_sp.o anglib.o base_h_phon.o read_dens_and_int.o choles.o admat.o dens_list.o eqm_hol_2ph.o
	$(XLFMPI) input_sp.o anglib.o base_h_phon.o read_dens_and_int.o choles.o admat.o dens_list.o eqm_hol_2ph.o ${OPT} -o ${EXE_PATH}${EXF}

input_sp.o: input_sp.f90
	${XLFMPI} -c input_sp.f90

anglib.o: anglib.f90
	${XLFMPI} -c anglib.f90

base_h_phon.o: base_h_phon.f90
	${XLFMPI} -c base_h_phon.f90

read_dens_and_int.o: read_dens_and_int.f90
	${XLFMPI} -c read_dens_and_int.f90

choles.o: choles.f
	${XLFMPI} -c choles.f

admat.o: admat.f90
	${XLFMPI} -c admat.f90

eqm_hol_2ph.o: eqm_hol_2ph.f90
	${XLFMPI} -c eqm_hol_2ph.f90

dens_list.o: dens_list.f90
	${XLFMPI} -c dens_list.f90

clean:
	rm -f *.o *.mod
	rm -f ${EXFMPI}




	
