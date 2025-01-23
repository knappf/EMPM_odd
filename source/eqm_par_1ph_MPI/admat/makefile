# compiler 
#XLF= ifort -O3 -parallel
#XLFMPI= mpiifort -g
# options
OPT= -L$MKLPATH -I$MKLINCLUDE -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread
EXE_PATH=../../${EMPM_BIN_DIR}
# exe file
EXF= eqm_part_1phon_admat
#
all: types_eqm.o input_sp.o anglib.o base_p_phon.o read_dens_and_int.o choles.o admat.o dens_list.o eqm_par_1ph.o
	$(XLFMPI) types_eqm.o input_sp.o anglib.o base_p_phon.o read_dens_and_int.o choles.o admat.o dens_list.o eqm_par_1ph.o ${OPT} -o ${EXE_PATH}${EXF}

types_eqm.o: types_eqm.f90
	${XLFMPI} -c types_eqm.f90

input_sp.o: input_sp.f90
	${XLFMPI} -c input_sp.f90

anglib.o: anglib.f90
	${XLFMPI} -c anglib.f90

base_p_phon.o: base_p_phon.f90
	${XLFMPI} -c base_p_phon.f90

read_dens_and_int.o: read_dens_and_int.f90
	${XLFMPI} -c read_dens_and_int.f90

choles.o: choles.f
	${XLFMPI} -c choles.f

admat.o: admat.f90
	${XLFMPI} -c admat.f90

eqm_par_1ph.o: eqm_par_1ph.f90
	${XLFMPI} -c eqm_par_1ph.f90

dens_list.o: dens_list.f90
	${XLFMPI} -c dens_list.f90

clean:
	rm -f *.o *.mod
	rm -f ${EXFMPI}




	
