n_cores:=4
#---------------------------------------------------------------
# Choose the MPI Fortran compiler with the corresponding options
#FF = mpiifort
# or
FF = mpifort
#FF = mpif90
FFLAGS = -fallow-argument-mismatch -O
#FFLAGS = -O
#---------------------------------------------------------------
# If FF = mpiifort, then uncomment
#       include 'mpif.h'
# and comment
#       use mpi
# in PSGen_test_mpi.f and vice versa if FF = mpif90 or FF = mpifort.
#---------------------------------------------------------------
OBJECTS =  inprms_ps.o engymom.o drivekk.o \
           fourmom.o distribs.o kincuts.o psnrmtest.o kinparams.o ntdisc.o \
           define_cuts.o prtinds.o qp_masses.o qrkinds.o \
           PSGen_test_mpi.o cutinfo.o set_indices.o \
           carlos.o impcuts.o \
           kinin.o kinin_qp.o boost.o boost_qp.o kinff.o kinff_qp.o \
           kinffg.o kinffg_qp.o mapakk.o mapbkk.o mapckk.o mapa_qp.o mapb_qp.o \
           mapc_qp.o \
	   kincls.o kintt3.o kintt_rest.o \
           kintt3_qp.o kintt_rest_qp.o kintt3a.o kintt3a_qp.o \
           kinschnl.o kingchnl.o \
           ini_photon.o \
           king2isr.o king2fsr.o king4isr.o king4fsr.o \
           king6isr.o king6fsr.o king6_33.o king8isr.o king8fsr.o king8_332.o \
           kin2_2.o kin3_2.o \
           kin4_22.o kin5_22.o kin6_33.o kin6_222.o kin6_32.o kin6_23.o kin7_222.o \
           kin7_33.o kin8_332.o kin8_323.o kin8_322.o kin8_2222.o kin8_233.o \
           kin8_232.o kin8_223.o kin9_332.o kin9_2222.o psncls.o psncls_qp.o \
           param_trans.o kintchnl.o kint4_2.o kint5_2.o kint6_22.o kint6_3.o \
           kint7_22.o kint7_3.o kint8_32.o kint8_22.o kint8_222.o kint9_222.o \
           kintmap.o kintmap_qp.o \
           psgen.o cs_sect.o ranlux.o logo_ps.o print_masses.o

.f.o: ; $(FF) $(FFLAGS) -c $*.f 

a.out: $(OBJECTS) makefile
	$(FF) -o a.out $(FFLAGS) $(OBJECTS) 
kinparams.o:  kinparams.f inprms_ps.f 
	$(FF) $(FFLAGS) -c kinparams.f 
prtinds.o:  prtinds.f kinparams.f 
	$(FF) $(FFLAGS) -c prtinds.f 
qp_masses.o:  qp_masses.f inprms_ps.f 
	$(FF) $(FFLAGS) -c qp_masses.f 
engymom.o:  engymom.f inprms_ps.f 
	$(FF) $(FFLAGS) -c engymom.f 
kincuts.o:  kincuts.f inprms_ps.f define_cuts.f 
	$(FF) $(FFLAGS) -c kincuts.f 
param_trans.o:  param_trans.f inprms_ps.f define_cuts.f set_indices.f
	$(FF) $(FFLAGS) -c param_trans.f 
PSGen_test_mpi.o: PSGen_test_mpi.f inprms_ps.f param_trans.f kinparams.f \
        drivekk.f kincuts.f cutinfo.f psnrmtest.f
	$(FF) $(FFLAGS) -c PSGen_test_mpi.f 
cutinfo.o: cutinfo.f kincuts.f prtinds.f define_cuts.f 
	$(FF) $(FFLAGS) -c cutinfo.f 
set_indices.o: set_indices.f prtinds.f
	$(FF) $(FFLAGS) -c set_indices.f 
define_cuts.o: define_cuts.f kincuts.f cutinfo.f
	$(FF) $(FFLAGS) -c define_cuts.f 
impcuts.o: impcuts.f inprms_ps.f kincuts.f prtinds.f define_cuts.f 
	$(FF) $(FFLAGS) -c impcuts.f 
kintchnl.o: kintchnl.f inprms_ps.f kinparams.f fourmom.f kincuts.f tchcalls.f
	$(FF) $(FFLAGS) -c kintchnl.f 
kinschnl.o: kinschnl.f inprms_ps.f kinparams.f fourmom.f
	$(FF) $(FFLAGS) -c kinschnl.f 
kingchnl.o: kingchnl.f inprms_ps.f kinparams.f fourmom.f
	$(FF) $(FFLAGS) -c kingchnl.f 
king2isr.o: king2isr.f fourmom.f engymom.f kincuts.f kinparams.f
	$(FF) $(FFLAGS) -c king2isr.f 
king2fsr.o: king2fsr.f fourmom.f engymom.f kincuts.f kinparams.f
	$(FF) $(FFLAGS) -c king2fsr.f 
king4isr.o: king4isr.f fourmom.f engymom.f kincuts.f kinparams.f
	$(FF) $(FFLAGS) -c king4isr.f 
king4fsr.o: king4fsr.f fourmom.f engymom.f kincuts.f kinparams.f
	$(FF) $(FFLAGS) -c king4fsr.f 
king6isr.o: king6isr.f fourmom.f kinparams.f
	$(FF) $(FFLAGS) -c king6isr.f 
king6fsr.o: king6fsr.f fourmom.f engymom.f kincuts.f kinparams.f
	$(FF) $(FFLAGS) -c king6fsr.f 
king6_33.o: king6_33.f fourmom.f kinparams.f
	$(FF) $(FFLAGS) -c king6_33.f 
king8isr.o: king8isr.f fourmom.f engymom.f kincuts.f kinparams.f
	$(FF) $(FFLAGS) -c king8isr.f 
king8fsr.o: king8fsr.f fourmom.f engymom.f kincuts.f kinparams.f
	$(FF) $(FFLAGS) -c king8fsr.f 
king8_332.o: king8_332.f fourmom.f kinparams.f
	$(FF) $(FFLAGS) -c king8_332.f 
ini_photon.o: ini_photon.f fourmom.f kincuts.f kinparams.f
	$(FF) $(FFLAGS) -c ini_photon.f 
kintt3.o: kintt3.f ntdisc.f
	$(FF) $(FFLAGS) -c kintt3.f 
kintt3_qp.o: kintt3_qp.f ntdisc.f 
	$(FF) $(FFLAGS) -c kintt3_qp.f 
kintt_rest.o: kintt_rest.f kinparams.f
	$(FF) $(FFLAGS) -c kintt_rest.f 
kintt_rest_qp.o: kintt_rest_qp.f kinparams.f
	$(FF) $(FFLAGS) -c kintt_rest_qp.f 
carlos.o: carlos.f drivekk.f distribs.f 
	$(FF) $(FFLAGS) -c carlos.f 
kincls.o: kincls.f kinparams.f fourmom.f kinschnl.f kintchnl.f \
          kingchnl.f kin_tests.f 
	$(FF) $(FFLAGS) -c kincls.f 
psncls.o: psncls.f inprms_ps.f kinparams.f psnrmtest.f 
	$(FF) $(FFLAGS) -c psncls.f 
psncls_qp.o: psncls_qp.f inprms_ps.f kinparams.f psnrmtest.f 
	$(FF) $(FFLAGS) -c psncls_qp.f 
kin2_2.o: kin2_2.f inprms_ps.f kinparams.f 
	$(FF) $(FFLAGS) -c kin2_2.f 
kin3_2.o: kin3_2.f inprms_ps.f kinparams.f 
	$(FF) $(FFLAGS) -c kin3_2.f 
kin4_22.o: kin4_22.f inprms_ps.f kinparams.f 
	$(FF) $(FFLAGS) -c kin4_22.f 
kin5_22.o: kin5_22.f inprms_ps.f kinparams.f 
	$(FF) $(FFLAGS) -c kin5_22.f 
kin6_33.o: kin6_33.f inprms_ps.f kinparams.f 
	$(FF) $(FFLAGS) -c kin6_33.f 
kin6_222.o: kin6_222.f inprms_ps.f kinparams.f 
	$(FF) $(FFLAGS) -c kin6_222.f 
kin6_32.o: kin6_32.f inprms_ps.f kinparams.f 
	$(FF) $(FFLAGS) -c kin6_32.f 
kin6_23.o: kin6_23.f inprms_ps.f kinparams.f 
	$(FF) $(FFLAGS) -c kin6_23.f 
kin7_222.o: kin7_222.f inprms_ps.f kinparams.f 
	$(FF) $(FFLAGS) -c kin7_222.f 
kin7_33.o: kin7_33.f inprms_ps.f kinparams.f 
	$(FF) $(FFLAGS) -c kin7_33.f 
kin8_332.o: kin8_332.f inprms_ps.f kinparams.f 
	$(FF) $(FFLAGS) -c kin8_332.f 
kin8_323.o: kin8_323.f inprms_ps.f kinparams.f 
	$(FF) $(FFLAGS) -c kin8_323.f 
kin8_322.o: kin8_322.f inprms_ps.f kinparams.f 
	$(FF) $(FFLAGS) -c kin8_322.f 
kin8_2222.o: kin8_2222.f inprms_ps.f kinparams.f 
	$(FF) $(FFLAGS) -c kin8_2222.f 
kin8_233.o: kin8_233.f inprms_ps.f kinparams.f 
	$(FF) $(FFLAGS) -c kin8_233.f 
kin8_232.o: kin8_232.f inprms_ps.f kinparams.f 
	$(FF) $(FFLAGS) -c kin8_232.f 
kin8_223.o: kin8_223.f inprms_ps.f kinparams.f 
	$(FF) $(FFLAGS) -c kin8_223.f 
kin9_332.o: kin9_332.f inprms_ps.f kinparams.f 
	$(FF) $(FFLAGS) -c kin9_332.f 
kin9_2222.o: kin9_2222.f inprms_ps.f kinparams.f 
	$(FF) $(FFLAGS) -c kin9_2222.f 
kint4_2.o: kint4_2.f inprms_ps.f kinparams.f 
	$(FF) $(FFLAGS) -c kint4_2.f 
kint5_2.o: kint5_2.f inprms_ps.f kinparams.f 
	$(FF) $(FFLAGS) -c kint5_2.f 
kint6_22.o: kint6_22.f inprms_ps.f kinparams.f 
	$(FF) $(FFLAGS) -c kint6_22.f 
kint6_3.o: kint6_3.f inprms_ps.f kinparams.f 
	$(FF) $(FFLAGS) -c kint6_3.f 
kint7_22.o: kint7_22.f inprms_ps.f kinparams.f 
	$(FF) $(FFLAGS) -c kint7_22.f 
kint7_3.o: kint7_3.f inprms_ps.f kinparams.f 
	$(FF) $(FFLAGS) -c kint7_3.f 
kint8_32.o: kint8_32.f inprms_ps.f kinparams.f 
	$(FF) $(FFLAGS) -c kint8_32.f 
kint8_22.o: kint8_22.f inprms_ps.f kinparams.f 
	$(FF) $(FFLAGS) -c kint8_22.f 
kint8_222.o: kint8_222.f inprms_ps.f kinparams.f 
	$(FF) $(FFLAGS) -c kint8_222.f 
kint9_222.o: kint9_222.f inprms_ps.f kinparams.f 
	$(FF) $(FFLAGS) -c kint9_222.f 
psgen.o: psgen.f kinparams.f kincuts.f fourmom.f fourmom.f
	$(FF) $(FFLAGS) -c psgen.f 
cs_sect.o: cs_sect.f inprms_ps.f kinparams.f drivekk.f \
           psnrmtest.f ntdisc.f fourmom.f
	$(FF) $(FFLAGS) -c cs_sect.f 
logo_ps.o: logo_ps.f process_nm.dat kinparams.f
	$(FF) $(FFLAGS) -c logo_ps.f 
print_masses.o: print_masses.f kinparams.f
	$(FF) $(FFLAGS) -c print_masses.f 
.PHONY: mc

mc: a.out
	mv tot_* ../test_output
	mv memos_* ../test_output
	time mpirun -np $(n_cores) ./a.out
	less tot_0_*
.PHONY: plot
#plot: 
#	gnuplot pion_ff.plo
#	ps2pdf test.ps
#	latex fig.tex
#	dvips -o fig.ps fig.dvi

