# compiler 
#XLF= ifort -O3 -parallel #-O3 -parallel#
#XLF= ifort -g
# options
OPT= -L$MKLPATH -I$MKLINCLUDE -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread

EXE_PATH=../../${EMPM_BIN_DIR}
# exe file
EXF= eqm_even_MPI
#
all: types_eqm.o input_sp.o anglib.o base.o read_admat.o choles.o cm_ort_svd.o ham.o dmat_test.o eqm_even.o 
	$(XLF) types_eqm.o input_sp.o anglib.o base.o read_admat.o choles.o cm_ort_svd.o ham.o dmat_test.o eqm_even.o ${OPT} -o ${EXE_PATH}${EXF}

input_sp.o: input_sp.f90
	${XLF} -c input_sp.f90

anglib.o: anglib.f90
	${XLF} -c anglib.f90

base.o: base.f90
	${XLF} -c base.f90

choles.o: choles.f90
	${XLF} -c choles.f90

ham.o: ham.f90
	${XLF} -c ham.f90

dmat_test.o: dmat_test.f
	${XLF} -c dmat_test.f

read_admat.o: read_admat.f90
	${XLF} -c read_admat.f90

types_eqm.o: types_eqm.f90
	${XLF} -c types_eqm.f90

cm_ort_svd.o: cm_ort_svd.f90
	${XLF} -c cm_ort_svd.f90

eqm_even.o: eqm_even.f90
	${XLF} -c eqm_even.f90

clean:
	rm -f *.o *.mod
	rm -f ${EXF}




	
