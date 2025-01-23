# compiler 
#XLF= ifort -O3 -parallel
#XLF= ifort -g
# options
OPT= -L$MKLPATH -I$MKLINCLUDE -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread
EXE_PATH=../../${EMPM_BIN_DIR}
# exe file
EXF= eqm_hole_1phon_MPI
#
all: types_eqm.o input_sp.o anglib.o base_h_phon.o read_dens_and_int.o read_admat.o choles.o cm_ort_svd.o dens_list.o ham.o eqm_hol_1ph.o
	$(XLF) types_eqm.o input_sp.o anglib.o base_h_phon.o read_dens_and_int.o read_admat.o choles.o cm_ort_svd.o dens_list.o ham.o eqm_hol_1ph.o ${OPT} -o ${EXE_PATH}${EXF}

types_eqm.o: types_eqm.f90
	${XLF} -c types_eqm.f90
	
input_sp.o: input_sp.f90
	${XLF} -c input_sp.f90

anglib.o: anglib.f90
	${XLF} -c anglib.f90

base_h_phon.o: base_h_phon.f90
	${XLF} -c base_h_phon.f90

read_dens_and_int.o: read_dens_and_int.f90
	${XLF} -c read_dens_and_int.f90

eqm_hol_1ph.o: eqm_hol_1ph.f90
	${XLF} -c eqm_hol_1ph.f90

dens_list.o: dens_list.f90
	${XLF} -c dens_list.f90

read_admat.o: read_admat.f90
	${XLF} -c read_admat.f90

choles.o: choles.f90
	${XLF} -c choles.f90

cm_ort_svd.o: cm_ort_svd.f90
	${XLF} -c cm_ort_svd.f90

ham.o: ham.f90
	${XLF} -c ham.f90

clean:
	rm -f *.o *.mod
	rm -f ${EXF}




	
