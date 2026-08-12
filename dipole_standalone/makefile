# Compiler and flags
FC = gfortran
FFLAGS = -O -J$(MODDIR) -I$(MODDIR)
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
PROG = runDipole

# Object files (you can extend this as you convert more files)
OBJECTS = $(OBJDIR)/printframe.o       \
          $(OBJDIR)/misc_utils_mod.o   \
          $(OBJDIR)/constants_mod.o    \
          $(OBJDIR)/global_mod.o       \
          $(OBJDIR)/mat.amp.o          \
          $(OBJDIR)/phasespace.o       \
          $(OBJDIR)/dipole.o           \
          $(OBJDIR)/integrand.o        \
          $(OBJDIR)/ol_interface.o     \
          $(OBJDIR)/pdf_lha.o          \
          $(OBJDIR)/filehandle.o       \
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
#
#
#
#
#
#
#
#
## Compile .f files (F77)
#
#
#
## Compiler and flags
#FC = gfortran
#FFLAGS = -O -w 
#LDFLAGS = -lc 
#export LHFLAGS = $(shell lhapdf-config --ldflags)
#
## Directories
#OBJDIR = objects
#SRCDIR = ..
#LOCALSRCDIR = .
#
## Program and object files
#PROG = runDipole
#
#OBJECTS = $(OBJDIR)/filehandle.o \
#          $(OBJDIR)/mat.amp.o    \
#          $(OBJDIR)/printframe.o \
#          $(OBJDIR)/misc.o       \
#          $(OBJDIR)/lum.o        \
#          $(OBJDIR)/pdf_lha.o    \
#          $(OBJDIR)/vsup.o       \
#          $(OBJDIR)/integrand.o  \
#          $(OBJDIR)/phasespace.o \
#          $(OBJDIR)/dipole.o	 
#
## Main target: Build executable from object files
#$(PROG): $(OBJDIR) $(OBJECTS) $(LOCALSRCDIR)/main.f
#	$(FC) $(LOCALSRCDIR)/main.f $(FFLAGS) -o $(PROG) $(OBJECTS) $(LIB) $(LHFLAGS) $(LDFLAGS)
#
## Create the object directory if not present.
#$(OBJDIR):
#	mkdir -p $(OBJDIR)
#
## Ensure the object directory exists and compile .f files into .o files
#$(OBJDIR)/%.o: $(SRCDIR)/%.f $(LOCALSRCDIR)/%.f| $(OBJDIR)
#	$(FC) $(FFLAGS) -c $< -o $@
#
## Pattern rule for all Fortran files in $(SRCDIR)
#$(OBJDIR)/%.o: $(SRCDIR)/%.f | $(OBJDIR)
#	$(FC) $(FFLAGS) -c $< -o $@
#
## For files in the current directory 
#$(OBJDIR)/%.o: %.f | $(OBJDIR)
#	$(FC) $(FFLAGS) -c $< -o $@
#
## Clean target
#clean:
#	rm -f $(OBJDIR)/*.o core *~ $(PROG)
#
#**********[ OLD MAKEFILE SCRIPT ]***********#
#FC = gfortran
#FFLAGS   = -O
#LDFLAGS = -L./../ -lfile
#export LHFLAGS=$(shell lhapdf-config --ldflags)
#LIB = -lm 
#LIB = 
##OBJECTS = Iterm.o matrixLO1.o couplings.o cmatrix.o PK.o dipole.o lum.o pdf_lha.o vsup.o mat.amp.o integrand.o misc.o phasespace.o
#OBJECTS = printframe.o dipole.o lum.o pdf_lha.o vsup.o mat.amp.o integrand.o misc.o phasespace.o
#PROG = runDipole
#$(PROG): $(OBJECTS) main.f  
#	$(FC) main.f $(FFLAGS) -o $(PROG) $(OBJECTS) $(LIB) $(LIBS) $(LHFLAGS) $(LDFLAGS) 
#lum.o: ../lum.f
#	$(FC) $(FFLAGS) -c ../lum.f -o lum.o
#vsup.o: ../vsup.f
#	$(FC) $(FFLAGS) -c ../vsup.f -o vsup.o
#pdf_lha.o: ../pdf_lha.f
#	$(FC) $(FFLAGS) -c ../pdf_lha.f -o pdf_lha.o
#mat.amp.o: ../mat.amp.f
#	$(FC) $(FFLAGS) -c ../mat.amp.f -o mat.amp.o
#misc.o: ../misc.f
#	$(FC) $(FFLAGS) -c ../misc.f -o misc.o
#phasespace.o: ../phasespace.f
#	$(FC) $(FFLAGS) -c ../phasespace.f -o phasespace.o
#printframe.o: ../printframe.f
#	$(FC) $(FFLAGS) -c ../printframe.f -o printframe.o
#clean: 
#	rm -f *.o core *~ *.x $(PROG)
#$(PROG): 
#$(PROG).o: $(PROG).f param_card.dat
