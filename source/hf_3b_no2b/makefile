FC = ifx
FCOPTS = -O3 -qopenmp
LD = ifx
LDOPTS = -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread
#LDOPTS = -L$MKLPATH -I$MKLINCLUDE -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread

EXE_PATH=../${EMPM_BIN_DIR}
EXF = HF_run
OBJS =HFB.o read_input.o dimm.o deallocate_all.o initialize_basis.o occupations.o kinetic.o v2b_bin.o v3b_no2b_bin.o interaction_binary.o v3body.o transition.o hfb_iteration.o transf_interaction.o HFB_energy.o hp_field.o hn_field.o Dp_field.o Dn_field.o diagonalization.o make_sp_levels.o make_densities.o check_densities.o R_val.o gauss_int.o TDA.o FNNtz.o VNNtz.o TDAtz.o make_qsp_levels.o MBPT_energy.o transit_calc.o transit_calc_spur.o phonon_density_calc.o diag_subspace.o spur_vec.o spur_vectz.o geom.o technical.o math.o Tcm_2b.o

$(EXF): $(OBJS)
	$(LD) $(LDOPTS) -o ${EXE_PATH}${EXF} $(OBJS)

HFB.o: HFB.f technical.o math.o
	$(FC) $(FCOPTS) -c HFB.f

read_input.o: read_input.f
	$(FC) $(FCOPTS) -c read_input.f

dimm.o: dimm.f
	$(FC) $(FCOPTS) -c dimm.f

deallocate_all.o: deallocate_all.f technical.o
	$(FC) $(FCOPTS) -c deallocate_all.f

#interaction.o: interaction.f technical.o geom.o Tcm_2b.o v3body.o v3b_no2b_bin.o
#	$(FC) $(FCOPTS) -c interaction.f

interaction_binary.o: interaction_binary.f technical.o geom.o Tcm_2b.o v3body.o v3b_no2b_bin.o
	$(FC) $(FCOPTS) -c interaction_binary.f

v2b_bin.o: v2b_bin.f90 technical.o Tcm_2b.o
	$(FC) $(FCOPTS) -c v2b_bin.f90

v3b_no2b_bin.o: v3b_no2b_bin.f90 technical.o
	$(FC) $(FCOPTS) -c v3b_no2b_bin.f90

transition.o: transition.f technical.o geom.o
	$(FC) $(FCOPTS) -c transition.f

v3body.o: v3body.f90 technical.o
	$(FC) $(FCOPTS) -c v3body.f90

hfb_iteration.o: hfb_iteration.f
	$(FC) $(FCOPTS) -c hfb_iteration.f

hfb_interaction.o: hfb_interaction.f technical.o geom.o
	$(FC) $(FCOPTS) -c hfb_interaction.f

initialize_basis.o: initialize_basis.f technical.o
	$(FC) $(FCOPTS) -c initialize_basis.f

occupations.o: occupations.f technical.o
	$(FC) $(FCOPTS) -c occupations.f

kinetic.o: kinetic.f technical.o
	$(FC) $(FCOPTS) -c kinetic.f

transf_interaction.o: transf_interaction.f technical.o
	$(FC) $(FCOPTS) -c transf_interaction.f

HFB_energy.o: HFB_energy.f technical.o
	$(FC) $(FCOPTS) -c HFB_energy.f

hp_field.o: hp_field.f technical.o
	$(FC) $(FCOPTS) -c hp_field.f

hn_field.o: hn_field.f technical.o
	$(FC) $(FCOPTS) -c hn_field.f

Dp_field.o: Dp_field.f technical.o
	$(FC) $(FCOPTS) -c Dp_field.f

Dn_field.o: Dn_field.f technical.o
	$(FC) $(FCOPTS) -c Dn_field.f

make_sp_levels.o: make_sp_levels.f technical.o
	$(FC) $(FCOPTS) -c make_sp_levels.f

make_densities.o: make_densities.f technical.o
	$(FC) $(FCOPTS) -c make_densities.f

check_densities.o: check_densities.f technical.o
	$(FC) $(FCOPTS) -c check_densities.f

R_val.o: R_val.f technical.o
	$(FC) $(FCOPTS) -c R_val.f

gauss_int.o: gauss_int.f technical.o
	$(FC) $(FCOPTS) -c gauss_int.f

diagonalization.o: diagonalization.f
	$(FC) $(FCOPTS) -c diagonalization.f

TDA.o: TDA.f math.o technical.o geom.o
	$(FC) $(FCOPTS) -c TDA.f

#pnTDApp.o: pnTDApp.f math.o technical.o geom.o
#	$(FC) $(FCOPTS) -c pnTDApp.f

#pnTDAhh.o: pnTDAhh.f math.o technical.o geom.o
#	$(FC) $(FCOPTS) -c pnTDAhh.f

#pnTDAph.o: pnTDAph.f math.o technical.o geom.o
#	$(FC) $(FCOPTS) -c pnTDAph.f

#pnTDAhp.o: pnTDAhp.f math.o technical.o geom.o
#	$(FC) $(FCOPTS) -c pnTDAhp.f

#ppTDApp.o: ppTDApp.f math.o technical.o geom.o
#	$(FC) $(FCOPTS) -c ppTDApp.f

#nnTDApp.o: nnTDApp.f math.o technical.o geom.o
#	$(FC) $(FCOPTS) -c nnTDApp.f

FNNtz.o: FNNtz.f technical.o
	$(FC) $(FCOPTS) -c FNNtz.f

VNNtz.o: VNNtz.f technical.o
	$(FC) $(FCOPTS) -c VNNtz.f

TDAtz.o: TDAtz.f FNNtz.o VNNtz.o math.o technical.o geom.o
	$(FC) $(FCOPTS) -c TDAtz.f

make_qsp_levels.o: make_qsp_levels.f technical.o
	$(FC) $(FCOPTS) -c make_qsp_levels.f

MBPT_energy.o: MBPT_energy.f technical.o geom.o
	$(FC) $(FCOPTS) -c MBPT_energy.f

transit_calc.o: transit_calc.f technical.o
	$(FC) $(FCOPTS) -c transit_calc.f

transit_calc_spur.o: transit_calc_spur.f technical.o
	$(FC) $(FCOPTS) -c transit_calc_spur.f

phonon_density_calc.o: phonon_density_calc.f math.o technical.o geom.o
	$(FC) $(FCOPTS) -c phonon_density_calc.f

diag_subspace.o: diag_subspace.f math.o technical.o geom.o
	$(FC) $(FCOPTS) -c diag_subspace.f

spur_vec.o: spur_vec.f math.o technical.o geom.o
	$(FC) $(FCOPTS) -c spur_vec.f

spur_vectz.o: spur_vectz.f math.o technical.o geom.o
	$(FC) $(FCOPTS) -c spur_vectz.f

technical.o: technical.f
	$(FC) $(FCOPTS) -c technical.f

math.o: math.f
	$(FC) $(FCOPTS) -c math.f

geom.o: geom.f
	$(FC) $(FCOPTS) -c geom.f

#geometric.o: geometric.f90
#	$(FC) $(FCOPTS) -c geometric.f90

Tcm_2b.o: geom.o Tcm_2b.f90 
	$(FC) $(FCOPTS) -c Tcm_2b.f90

#.f90.f.o:
#	$(FC) $(FCOPTS) -c $<

clean:
	rm -f *.o *.mod $(EXENAME)

#clobber:
#	rm -f $(EXENAME)
