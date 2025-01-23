# compiler 
#XLF= ifort -O3 -parallel #-O3 -parallel
#XLFMPI= mpiifort -g
# options
OPT= -L$MKLPATH -I$MKLINCLUDE -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread

EXE_PATH=../../${EMPM_BIN_DIR}
# exe file
EXF= eqm_even_admat
#
all: input_sp.o anglib.o base.o admat_cp_par.o eqm_even.o 
	$(XLFMPI) input_sp.o anglib.o base.o admat_cp_par.o eqm_even.o ${OPT} -o ${EXE_PATH}${EXF}

input_sp.o: input_sp.f90
	${XLFMPI} -c input_sp.f90

anglib.o: anglib.f90
	${XLFMPI} -c anglib.f90

base.o: base.f90
	${XLFMPI} -c base.f90

admat_cp_par.o: admat_cp_par.f90
	${XLFMPI} -c admat_cp_par.f90

eqm_even.o: eqm_even.f90
	${XLFMPI} -c eqm_even.f90

clean:
	rm -f *.o *.mod
	rm -f ${EXF}




	
