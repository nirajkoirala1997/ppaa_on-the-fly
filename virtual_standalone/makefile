# Compiler and flags
FC = gfortran
FFLAGS = -O -w -J$(MODDIR) -I$(MODDIR)
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

# Program and object files
PROG = runVirtual

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
          $(OBJDIR)/vsup.o         

# OpenLoops linking
OPENLOOPS_FLAGS = -Wl,-rpath=$(OPENLOOPS_LIB) -I$(OPENLOOPS_INCLUDE) -L$(OPENLOOPS_LIB) -lopenloops

# Main target
$(PROG): $(OBJDIR) $(OBJECTS) $(LOCALSRCDIR)/main.f90
	$(FC) $(LOCALSRCDIR)/main.f90 $(FFLAGS) -o $(PROG) $(OBJECTS) $(LIB) $(LHFLAGS) $(LDFLAGS) $(OPENLOOPS_FLAGS)

# Create object directory if not present
$(OBJDIR):
	mkdir -p $(OBJDIR)
	mkdir -p $(MODDIR)

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









## Compiler and flags
#FC = gfortran
#FFLAGS = -O -w
#LDFLAGS = -lc -L/home/niraj/1TB-Disc/Workspace-IITG/Packages/Install/Cuba-4.2.2 -lcuba 
#export LHFLAGS = $(shell lhapdf-config --ldflags)
#
## OpenLoops paths
#OPENLOOPS_DIR = /home/niraj/1TB-Disc/Workspace-IITG/Packages/Install/OpenLoops
#OPENLOOPS_INCLUDE = $(OPENLOOPS_DIR)/lib_src/openloops/mod
#OPENLOOPS_LIB = $(OPENLOOPS_DIR)/lib
#
## Directories
#OBJDIR = objects
#SRCDIR = ..
#LOCALSRCDIR = .
#
## Program and object files
#PROG = runLO
#
#OBJECTS = $(OBJDIR)/mat.amp.o    \
#          $(OBJDIR)/printframe.o \
#          $(OBJDIR)/filehandle.o \
#          $(OBJDIR)/misc.o       \
#          $(OBJDIR)/lum.o        \
#          $(OBJDIR)/pdf_lha.o    \
#          $(OBJDIR)/vsup.o       \
#          $(OBJDIR)/integrand.o  \
#          $(OBJDIR)/ol_interface.o  \
#          $(OBJDIR)/phasespace.o
#
## Add OpenLoops library to the linking process
#OPENLOOPS_FLAGS = -Wl,-rpath=$(OPENLOOPS_LIB) -I$(OPENLOOPS_INCLUDE) -L$(OPENLOOPS_LIB) -lopenloops
#
## Main target: Build executable from object files
#$(PROG): $(OBJDIR) $(OBJECTS) $(LOCALSRCDIR)/main.f
#	$(FC) $(LOCALSRCDIR)/main.f $(FFLAGS) -o $(PROG) $(OBJECTS) $(LIB) $(LHFLAGS) $(LDFLAGS) $(OPENLOOPS_FLAGS)
#
## Create the object directory if not present.
#$(OBJDIR):
#	mkdir -p $(OBJDIR)
#
## Ensure the object directory exists and compile .f files into .o files
#$(OBJDIR)/%.o: $(SRCDIR)/%.f | $(OBJDIR)
#	$(FC) $(FFLAGS) -I$(OPENLOOPS_INCLUDE) -c $< -o $@
#
## For files in the current directory
#$(OBJDIR)/%.o: %.f | $(OBJDIR)
#	$(FC) $(FFLAGS) -I$(OPENLOOPS_INCLUDE) -c $< -o $@
#
## Clean target
#clean:
#	rm -f $(OBJDIR)/*.o core *~ $(PROG)
