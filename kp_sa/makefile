#######################################################################################
###  runKP -- the KP term.       No dependency tracking: the object order is the   ###
###  module dependency graph, `make -j` is not safe, `make clean` after a change.   ###
#######################################################################################

PROG      = runKP
OBJDIR    = objects
MODDIR    = modules
PSGEN_DIR = ../PSGen
COMMON    = ../common

include $(COMMON)/config.mk

override F95FLAGS += -ffree-line-length-none

# Sources shared with the other programs live at the repository root.
$(OBJDIR)/cuts_mod.o: ../cuts_mod.f90 | $(OBJDIR)
	$(FC) $(F95FLAGS) $(INC_ORDER) -c $< -o $@

$(OBJDIR)/psgen_interface.o: ../psgen_interface.f90 | $(OBJDIR)
	$(FC) $(F95FLAGS) $(INC_ORDER) -c $< -o $@

$(OBJDIR)/olp_mod.o: ../olp_mod.f90 | $(OBJDIR)
	$(FC) $(F95FLAGS) $(INC_ORDER) -c $< -o $@

include $(COMMON)/common.mk

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

OBJECTS = $(PSGEN_OBJS) $(COMMON_OBJS) $(OBJDIR)/olp_mod.o $(OBJDIR)/dy_kp_kernels_mod.o $(OBJDIR)/integrand.o

all: $(PROG)

$(PROG): $(OBJDIR) $(OBJECTS) main.f90 $(MAKEFILE_LIST)
	$(FC) main.f90 $(F95FLAGS) $(INC_ORDER) -o $(PROG) $(OBJECTS) $(LHFLAGS) $(LDFLAGS) $(LIB_ORDER)

clean:
	rm -rf $(OBJDIR) $(MODDIR) $(PROG) output_cll .collier

.PHONY: all clean
