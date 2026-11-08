#######################################################################################
###  Shared build rules for every study in this repository.                         ###
###                                                                                 ###
###  A study directory (decay/, onshell/) holds only what depends on the process:   ###
###  its generated PSGen/ kinematics, olp_mod.f90, integrand.f90, main.f90 and its  ###
###  tests. Everything that does not know which process is running lives in         ###
###  common/ and is compiled -- not archived -- into each study's own objects/,     ###
###  so the two studies never share a .mod file or a build directory.               ###
###                                                                                 ###
###  psgen_interface.f90 is the reason this is compiled per study rather than       ###
###  built once into a library: it takes nfsp/nexp/nkin from whichever generated    ###
###  PSGen/ it is compiled against, so its object is process-specific even though   ###
###  its source is not.                                                             ###
###                                                                                 ###
###  A study makefile must define OBJDIR, MODDIR and PSGEN_DIR before including     ###
###  this file, and must place $(COMMON_OBJS) after its $(PSGEN_OBJS) and before    ###
###  its own objects. THERE IS NO DEPENDENCY TRACKING -- the order IS the module    ###
###  dependency graph, here as in the study makefiles, and `make -j` is not safe.   ###
#######################################################################################

# A study makefile's `include` lands the rules below AHEAD of its own `all:`
# target, and make takes the first target it sees as the default goal. Without
# this line a bare `make` in a study directory would run the $(OBJDIR) rule and
# stop after mkdir, reporting success.
.DEFAULT_GOAL := all

# `?=` and not `=`: a study one level below the repository root gets ../common
# for free, while a program nested deeper -- onshell/kp_sa/ -- sets its
# own ../../common BEFORE this include and keeps it. Every pattern rule below
# resolves through this variable, and they are expanded when the rule is READ,
# so a study cannot fix it up afterwards.
COMMON ?= ../common

#========== S H A R E D  O B J E C T S ==========#
# Order is the dependency graph:
#   card_parser <- run_card <- {cuts, olp_mod}
#   scales_mod <- pdfvar_mod    (both consumed by histogram_mod and any
#                                integrand that replays a variation grid)
#   {histogram, observables} <- analysis <- checkpoint
# psgen_interface goes last of these because it `use`s the PSGen modules, which
# the study's $(PSGEN_OBJS) must therefore have built already.
COMMON_OBJS = $(OBJDIR)/card_parser_mod.o  \
              $(OBJDIR)/run_card_mod.o     \
              $(OBJDIR)/run_state_mod.o    \
              $(OBJDIR)/scales_mod.o       \
              $(OBJDIR)/pdfvar_mod.o       \
              $(OBJDIR)/observables_mod.o  \
              $(OBJDIR)/photon_isolation_mod.o \
              $(OBJDIR)/cuts_mod.o         \
              $(OBJDIR)/histogram_mod.o    \
              $(OBJDIR)/analysis_mod.o     \
              $(OBJDIR)/checkpoint_mod.o   \
              $(OBJDIR)/psgen_interface.o  \
              $(OBJDIR)/pdf_lha.o          \
              $(OBJDIR)/vsup.o

#========== R U L E S ==========#
# Study-local sources win over common ones: make tries the pattern rules in the
# order they are defined, and only one of the two files ever exists. A study that
# needs its own copy of a common module drops it in its own directory and it is
# picked up here without touching this file.

# Study-local free-form
$(OBJDIR)/%.o: %.f90 | $(OBJDIR)
	$(FC) $(F95FLAGS) $(INC_ORDER) -c $< -o $@

# Study-local fixed-form
$(OBJDIR)/%.o: %.f | $(OBJDIR)
	$(FC) $(FFLAGS) -c $< -o $@

# Generated PSGen kinematics, fixed-form. -I is needed: kincls.f includes
# kin_tests.f and logo_ps.f includes process_nm.dat.
$(OBJDIR)/%.o: $(PSGEN_DIR)/%.f | $(OBJDIR)
	$(FC) $(FFLAGS) -I$(PSGEN_DIR) -c $< -o $@

# Shared free-form
$(OBJDIR)/%.o: $(COMMON)/%.f90 | $(OBJDIR)
	$(FC) $(F95FLAGS) $(INC_ORDER) -c $< -o $@

# Shared fixed-form (vsup.f, pdf_lha.f)
$(OBJDIR)/%.o: $(COMMON)/%.f | $(OBJDIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJDIR):
	mkdir -p $(OBJDIR) $(MODDIR)
