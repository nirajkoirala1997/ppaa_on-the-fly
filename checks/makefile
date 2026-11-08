#######################################################################################
###  check_recola_points -- Recola amplitudes at fixed phase-space points.          ###
###  There is no dependency tracking: `make clean` after changing a module.         ###
#######################################################################################

PROG   = check_recola_points
OBJDIR = objects
MODDIR = modules
COMMON = ../common

include $(COMMON)/config.mk

override F95FLAGS += -ffree-line-length-none

OBJECTS = $(OBJDIR)/card_parser_mod.o $(OBJDIR)/run_card_mod.o $(OBJDIR)/olp_mod.o

all: $(PROG)

$(OBJDIR):
	mkdir -p $(OBJDIR) $(MODDIR)

$(OBJDIR)/olp_mod.o: ../olp_mod.f90 | $(OBJDIR)
	$(FC) $(F95FLAGS) $(INC_ORDER) -c $< -o $@

$(OBJDIR)/%.o: $(COMMON)/%.f90 | $(OBJDIR)
	$(FC) $(F95FLAGS) $(INC_ORDER) -c $< -o $@

$(PROG): $(OBJECTS) $(PROG).f90
	$(FC) $(PROG).f90 $(F95FLAGS) $(INC_ORDER) -o $(PROG) $(OBJECTS) $(LDFLAGS) $(LIB_ORDER)

clean:
	rm -rf $(OBJDIR) $(MODDIR) $(PROG) output_cll .collier

.PHONY: all clean
