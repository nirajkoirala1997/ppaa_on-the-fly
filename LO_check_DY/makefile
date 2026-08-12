# Compiler and flags
FC = gfortran
#FFLAGS = -O -w -ffree-form -Wall -J$(MODDIR) -I$(MODDIR) -I$(PSGEN_DIR)
FFLAGS = -O -w -J$(MODDIR) -I$(MODDIR) # -I$(PSGEN_DIR)
F95FLAGS = -O -Wall -std=f2008  -J$(MODDIR) -I$(MODDIR)

LDFLAGS = -lc -L/home/niraj/1TB-Disc/Workspace-IITG/Packages/Install/Cuba-4.2.2 -lcuba
export LHFLAGS = $(shell lhapdf-config --ldflags)

# OpenLoops paths
OPENLOOPS_DIR = /home/niraj/1TB-Disc/Workspace-IITG/Packages/Install/OpenLoops
OPENLOOPS_INCLUDE = $(OPENLOOPS_DIR)/lib_src/openloops/mod
OPENLOOPS_LIB = $(OPENLOOPS_DIR)/lib

# Directories
OBJDIR = objects
MODDIR = modules 
SRCDIR = ..
LOCALSRCDIR = .

# Path to PSGen files
PSGEN_DIR = PSGen

# PSGen Fortran sources and corresponding objects
#PSGEN_SRCS = $(wildcard $(PSGEN_DIR)/*.f)
#PSGEN_OBJS = $(patsubst $(PSGEN_DIR)/%.f,$(OBJDIR)/%.o,$(PSGEN_SRCS))

# List of object files located inside the PSGen directory
PSGEN_OBJS = $(OBJDIR)/inprms_ps.o \
             $(OBJDIR)/kinparams.o \
             $(OBJDIR)/engymom.o \
             $(OBJDIR)/fourmom.o \
             $(OBJDIR)/kincuts.o \
             $(OBJDIR)/psnrmtest.o \
             $(OBJDIR)/drivekk.o \
             $(OBJDIR)/distribs.o \
             $(OBJDIR)/ntdisc.o \
             $(OBJDIR)/define_cuts.o \
             $(OBJDIR)/prtinds.o \
             $(OBJDIR)/qp_masses.o \
             $(OBJDIR)/qrkinds.o \
             $(OBJDIR)/cutinfo.o \
             $(OBJDIR)/set_indices.o \
             $(OBJDIR)/carlos.o \
             $(OBJDIR)/impcuts.o \
             $(OBJDIR)/kinin.o \
             $(OBJDIR)/kinin_qp.o \
             $(OBJDIR)/boost.o \
             $(OBJDIR)/boost_qp.o \
             $(OBJDIR)/kinff.o \
             $(OBJDIR)/kinff_qp.o \
             $(OBJDIR)/kinffg.o \
             $(OBJDIR)/kinffg_qp.o \
             $(OBJDIR)/mapakk.o \
             $(OBJDIR)/mapbkk.o \
             $(OBJDIR)/mapckk.o \
             $(OBJDIR)/mapa_qp.o \
             $(OBJDIR)/mapb_qp.o \
             $(OBJDIR)/mapc_qp.o \
             $(OBJDIR)/kincls.o \
             $(OBJDIR)/kintt3.o \
             $(OBJDIR)/kintt_rest.o \
             $(OBJDIR)/kintt3_qp.o \
             $(OBJDIR)/kintt_rest_qp.o \
             $(OBJDIR)/kintt3a.o \
             $(OBJDIR)/kintt3a_qp.o \
             $(OBJDIR)/kinschnl.o \
             $(OBJDIR)/kingchnl.o \
             $(OBJDIR)/ini_photon.o \
             $(OBJDIR)/king2isr.o \
             $(OBJDIR)/king2fsr.o \
             $(OBJDIR)/king4isr.o \
             $(OBJDIR)/king4fsr.o \
             $(OBJDIR)/king6isr.o \
             $(OBJDIR)/king6fsr.o \
             $(OBJDIR)/king6_33.o \
             $(OBJDIR)/king8isr.o \
             $(OBJDIR)/king8fsr.o \
             $(OBJDIR)/king8_332.o \
             $(OBJDIR)/kin2_2.o \
             $(OBJDIR)/kin3_2.o \
             $(OBJDIR)/kin4_22.o \
             $(OBJDIR)/kin5_22.o \
             $(OBJDIR)/kin6_33.o \
             $(OBJDIR)/kin6_222.o \
             $(OBJDIR)/kin6_32.o \
             $(OBJDIR)/kin6_23.o \
             $(OBJDIR)/kin7_222.o \
             $(OBJDIR)/kin7_33.o \
             $(OBJDIR)/kin8_332.o \
             $(OBJDIR)/kin8_323.o \
             $(OBJDIR)/kin8_322.o \
             $(OBJDIR)/kin8_2222.o \
             $(OBJDIR)/kin8_233.o \
             $(OBJDIR)/kin8_232.o \
             $(OBJDIR)/kin8_223.o \
             $(OBJDIR)/kin9_332.o \
             $(OBJDIR)/kin9_2222.o \
             $(OBJDIR)/psncls.o \
             $(OBJDIR)/psncls_qp.o \
             $(OBJDIR)/param_trans.o \
             $(OBJDIR)/kintchnl.o \
             $(OBJDIR)/kint4_2.o \
             $(OBJDIR)/kint5_2.o \
             $(OBJDIR)/kint6_22.o \
             $(OBJDIR)/kint6_3.o \
             $(OBJDIR)/kint7_22.o \
             $(OBJDIR)/kint7_3.o \
             $(OBJDIR)/kint8_32.o \
             $(OBJDIR)/kint8_22.o \
             $(OBJDIR)/kint8_222.o \
             $(OBJDIR)/kint9_222.o \
             $(OBJDIR)/kintmap.o \
             $(OBJDIR)/kintmap_qp.o \
             $(OBJDIR)/psgen.o \
             $(OBJDIR)/cs_sect.o \
             $(OBJDIR)/ranlux.o \
             $(OBJDIR)/logo_ps.o \
             $(OBJDIR)/print_masses.o

#             $(PSGEN_DIR)/PSGen_test_mpi.o \

# Program and object files
PROG = runLO

# Object files (you can extend this as you convert more files)
OBJECTS = $(OBJDIR)/printframe.o     \
          $(OBJDIR)/misc_utils_mod.o \
          $(OBJDIR)/constants_mod.o  \
          $(OBJDIR)/global_mod.o     \
          $(OBJDIR)/mat.amp.o        \
          $(OBJDIR)/phasespace.o     \
          $(OBJDIR)/integrand.o      \
          $(OBJDIR)/ol_interface.o   \
          $(OBJDIR)/pdf_lha.o        \
          $(OBJDIR)/filehandle.o     \
	  $(OBJDIR)/vsup.o     \
          $(PSGEN_OBJS)

# OpenLoops linking
OPENLOOPS_FLAGS = -Wl,-rpath=$(OPENLOOPS_LIB) -I$(OPENLOOPS_INCLUDE) -L$(OPENLOOPS_LIB) -lopenloops

# Main target
$(PROG): $(OBJDIR) $(OBJECTS) $(LOCALSRCDIR)/main.f90
	$(FC) $(LOCALSRCDIR)/main.f90 $(FFLAGS) -o $(PROG) $(OBJECTS) $(LIB) $(LHFLAGS) $(LDFLAGS) $(OPENLOOPS_FLAGS)

# Create object directory if not present
$(OBJDIR):
	mkdir -p $(OBJDIR)
	mkdir -p $(MODDIR)

# Compile .f files from PSGen
$(OBJDIR)/%.o: $(PSGEN_DIR)/%.f | $(OBJDIR)
	$(FC) $(FFLAGS) -I$(OPENLOOPS_INCLUDE) -c $< -o $@

# Compile .f files (F77)
$(OBJDIR)/%.o: $(SRCDIR)/%.f | $(OBJDIR)
	$(FC) $(FFLAGS) -I$(OPENLOOPS_INCLUDE) -c $< -o $@

$(OBJDIR)/%.o: %.f | $(OBJDIR)
	$(FC) $(FFLAGS) -I$(OPENLOOPS_INCLUDE) -c $< -o $@

# Compile .f90 or .f95 files (F95)
$(OBJDIR)/%.o: $(SRCDIR)/%.f90 | $(OBJDIR)
	$(FC) $(F95FLAGS) -I$(OPENLOOPS_INCLUDE) -c $< -o $@

$(OBJDIR)/%.o: %.f90 | $(OBJDIR)
	$(FC) $(F95FLAGS) -I$(OPENLOOPS_INCLUDE) -c $< -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.f95 | $(OBJDIR)
	$(FC) $(F95FLAGS) -I$(OPENLOOPS_INCLUDE) -c $< -o $@

$(OBJDIR)/%.o: %.f95 | $(OBJDIR)
	$(FC) $(F95FLAGS) -I$(OPENLOOPS_INCLUDE) -c $< -o $@

# Clean
clean:
	rm -f $(OBJDIR)/*.o core *~ $(PROG)

