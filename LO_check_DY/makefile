#######################################################################################
###  runLO -- leading-order  p p -> e+ e-  (Drell-Yan)                              ###
###                                                                                 ###
###  Phase space: PSGen (Kolodziej), the generated 'u u~ -> e+ e-' kinematics in    ###
###               PSGen/.  PSGen supplies the MAPPING only -- it never draws a      ###
###               random number, so VEGAS (vsup.f) stays in charge of the loop and  ###
###               sees every integration variable, phase-space or not.              ###
###  Amplitudes:  OpenLoops.                                                        ###
#######################################################################################

FC = gfortran

OBJDIR      = objects
MODDIR      = modules
SRCDIR      = ..
LOCALSRCDIR = .
PSGEN_DIR   = PSGen

## Dependencies, resolved through pkgtool -- no hardcoded prefixes.
PKGTOOL = pkgtool

#========== L H A P D F ==========#
LHAPDF_LIBDIR := $(shell $(PKGTOOL) lhapdf --libdir)
LHFLAGS       := $(shell $(PKGTOOL) lhapdf --libs) -Wl,-rpath,$(LHAPDF_LIBDIR)

#========== O P E N L O O P S ==========#
OL_LIBDIR := $(shell $(PKGTOOL) openloops --libdir)
OL_PREFIX := $(abspath $(OL_LIBDIR)/..)
# pkgtool --incdir gives the C headers; `use openloops` needs the Fortran .mod tree.
# Override on the command line if your installation puts them elsewhere.
OL_MODDIR ?= $(OL_PREFIX)/lib_src/openloops/mod
OL_INC    := -I$(OL_MODDIR)
OL_LIB    := $(shell $(PKGTOOL) openloops --libs) -Wl,-rpath,$(OL_LIBDIR)
# Compiled in so the binary finds proclib/ without $OpenLoopsPath set in the shell.
OL_DEF    := -DOL_PREFIX=\"$(OL_PREFIX)\"

LDFLAGS = -lc

# Cuba is deliberately NOT linked: nothing in this tree calls it. The integrator
# is vsup/VEGAS.

#========== F L A G S ==========#
# -w on the F77 rule: the PSGen sources are fixed-form legacy and will not
# compile clean. The .f90 rule stays strict.
FFLAGS   = -O -w -std=legacy -J$(MODDIR) -I$(MODDIR)
F95FLAGS = -O -Wall -std=f2008 -cpp $(OL_DEF) -J$(MODDIR) -I$(MODDIR)

PROG = runLO

#========== P S G e n ==========#
# Trimmed to the set reachable at this multiplicity (nkin = 1):
#   psgen -> kincls -> kinschnl -> kin2_2 -> kinff
# plus psncls, kinin, param_trans and the cut/index bookkeeping. The other ~74
# generated files cover higher multiplicities and t-/g-channels; regenerate them
# from PSGen_1.1/ when a second process is added. carlos.f is gone on purpose --
# that is PSGen's own multichannel integrator, which VEGAS replaces here.
#
# Modules first: this makefile has no dependency tracking, so the order below IS
# the build order.
PSGEN_OBJS = $(OBJDIR)/inprms_ps.o   \
             $(OBJDIR)/qrkinds.o     \
             $(OBJDIR)/fourmom.o     \
             $(OBJDIR)/kincuts.o     \
             $(OBJDIR)/psnrmtest.o   \
             $(OBJDIR)/ntdisc.o      \
             $(OBJDIR)/drivekk.o     \
             $(OBJDIR)/kinparams.o   \
             $(OBJDIR)/prtinds.o     \
             $(OBJDIR)/qp_masses.o   \
             $(OBJDIR)/define_cuts.o \
             $(OBJDIR)/set_indices.o \
             $(OBJDIR)/cutinfo.o     \
             $(OBJDIR)/impcuts.o     \
             $(OBJDIR)/logo_ps.o     \
             $(OBJDIR)/print_masses.o\
             $(OBJDIR)/param_trans.o \
             $(OBJDIR)/kinin.o       \
             $(OBJDIR)/kinin_qp.o    \
             $(OBJDIR)/kinff.o       \
             $(OBJDIR)/kinff_qp.o    \
             $(OBJDIR)/kin2_2.o      \
             $(OBJDIR)/kinschnl.o    \
             $(OBJDIR)/kincls.o      \
             $(OBJDIR)/psncls.o      \
             $(OBJDIR)/psgen.o

#========== P R O J E C T ==========#
OBJECTS = $(PSGEN_OBJS)              \
          $(OBJDIR)/constants_mod.o  \
          $(OBJDIR)/global_mod.o     \
          $(OBJDIR)/input_card_mod.o \
          $(OBJDIR)/misc_utils_mod.o \
          $(OBJDIR)/observables_mod.o\
          $(OBJDIR)/histogram_mod.o  \
          $(OBJDIR)/analysis_mod.o   \
          $(OBJDIR)/printframe.o     \
          $(OBJDIR)/phasespace.o     \
          $(OBJDIR)/mat.amp.o        \
          $(OBJDIR)/ol_interface.o   \
          $(OBJDIR)/pdf_lha.o        \
          $(OBJDIR)/vsup.o           \
          $(OBJDIR)/psgen_interface.o\
          $(OBJDIR)/integrand.o

$(PROG): $(OBJDIR) $(OBJECTS) $(LOCALSRCDIR)/main.f90
	$(FC) $(LOCALSRCDIR)/main.f90 $(F95FLAGS) $(OL_INC) -o $(PROG) \
	      $(OBJECTS) $(LHFLAGS) $(LDFLAGS) $(OL_LIB)

$(OBJDIR):
	mkdir -p $(OBJDIR) $(MODDIR)

# PSGen fixed-form sources
$(OBJDIR)/%.o: $(PSGEN_DIR)/%.f | $(OBJDIR)
	$(FC) $(FFLAGS) -I$(PSGEN_DIR) -c $< -o $@

# Shared fixed-form sources one level up
$(OBJDIR)/%.o: $(SRCDIR)/%.f | $(OBJDIR)
	$(FC) $(FFLAGS) -c $< -o $@

# Shared free-form sources one level up
$(OBJDIR)/%.o: $(SRCDIR)/%.f90 | $(OBJDIR)
	$(FC) $(F95FLAGS) $(OL_INC) -c $< -o $@

# Local free-form sources
$(OBJDIR)/%.o: %.f90 | $(OBJDIR)
	$(FC) $(F95FLAGS) $(OL_INC) -c $< -o $@

.PHONY: clean
clean:
	rm -rf $(OBJDIR) $(MODDIR) $(PROG)
