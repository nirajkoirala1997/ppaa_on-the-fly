#######################################################################################
###   Toolchain and external dependencies, shared by every study in this            ###
###   repository.  Nothing here knows which process is being computed.              ###
###                                                                                 ###
###   Phase space: PSGen (Kolodziej).  PSGen supplies the MAPPING only -- it        ###
###                never draws a random number, so VEGAS (vsup.f) stays in charge   ###
###                of the loop and sees every integration variable, the channel     ###
###                index included.  The generated kinematics are per study, in      ###
###                that study's PSGen/; the generator itself is common/PSGen_1.1.   ###
###   Integrator: VEGAS (vsup.f).                                                   ###
###   Amplitudes: Recola 1 -- ggxy's patched fork, vendored in                     ###
###               external/recola-ggxy/ -- is configured here for every study.  A    ###
###               study that also links another provider adds it in its own         ###
###               config fragment and appends to LIB_ORDER -- see                   ###
###               decay/config_mg5.mk, which carries the MadGraph library that      ###
###               used to sit in this file.                                         ###
#######################################################################################

FC = gfortran

## Dependencies. LHAPDF and handyG are resolved through pkgtool; Recola is NOT --
## it is built from this repository's own sources (see the RECOLA block below).

PKGTOOL = pkgtool

#========== L H A P D F ==========#
LHAPDF_LIBDIR := $(shell $(PKGTOOL) lhapdf --libdir)
LHFLAGS       := $(shell $(PKGTOOL) lhapdf --libs) -Wl,-rpath,$(LHAPDF_LIBDIR)

#========== R E C O L A  (required) ==========#
# RECOLA 1 -- NOT recola2, which pkgtool also carries
# and which this build used until now. The switch is deliberate: Recola 1
# exposes the model directly through its own API, where Recola 2 hides it
# behind a compiled model file. What that buys, and what olp_mod.f90 is now
# written against:
#
#   * INTERMEDIATE-STATE SELECTION IN THE PROCESS STRING. Round brackets
#     after a particle name restrict the amplitude to the diagrams in which
#     it appears as an internal line -- 'g g -> W+(e+ nu_e) W-(e- nu_e~)' is
#     the W+W- resonant subset, the same selection MadGraph's
#     `g g > w+ w- > e+ e- ve ve~` decay-chain syntax makes. Nothing else in
#     the build has to change to use it: it is just a different string.
#   * SINGLE VERTICES CAN BE SCALED OR SWITCHED OFF, at run time, through
#     switchoff_coupling3_rcl / switchoff_coupling4_rcl / scale_coupling3_rcl
#     / scale_coupling4_rcl.
#   * gs POWERS ONLY. Recola 1 selects coupling orders by power of gs alone
#     and infers the electroweak order, so the [QCD,QED] tuple -- and with it
#     the [2,4]-vs-[4,8] amplitude/squared-amplitude trap that cost this
#     project a debugging round -- does not exist in this API.
#
# The two Recolas export a Fortran module with the SAME name (`recola`) and
# the same symbol names, so they can never be linked together. There is no
# recola2 entry anywhere in this file any more, and reintroducing one would
# not work alongside this.
#
# WHICH Recola 1: ggxy's patched fork (iopt_loopCoef, force_zeroLO,
# dynamic_params), vendored with the COLLIER 1.2.8 it links in
# external/recola-ggxy/ -- see PROVENANCE.md there. One Recola for every
# program: lo_sa/kp_sa use it at its stock defaults, rs_sa switches the
# patches on through onshell/recola_ggxy_interface/recola_qp_mod.f90. It is
# built once with external/recola-ggxy/build.sh into a gitignored install/,
# and located from THIS file's own path, since config.mk is included from
# different depths. --disable-new-dtags: librecola.so needs libcollier.so
# transitively, which a DT_RUNPATH rpath does not reach.
CONFIG_MK_DIR  := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
RECOLA_PREFIX  := $(abspath $(CONFIG_MK_DIR)/../../ggZA/external/recola-ggxy/install)
RECOLA_LIBDIR  := $(RECOLA_PREFIX)/lib
RECOLA_INCFLAG := -I$(RECOLA_PREFIX)/include/recola
RECOLA_LIBS    := -L$(RECOLA_LIBDIR) -lrecola -lcollier

ifeq ($(wildcard $(RECOLA_LIBDIR)/librecola.so),)
$(warning Recola not built: run external/recola-ggxy/build.sh first ($(RECOLA_LIBDIR)/librecola.so missing))
endif

RECOLA_FLAGS = -Wl,--disable-new-dtags -Wl,-rpath,$(RECOLA_LIBDIR) $(RECOLA_INCFLAG) $(RECOLA_LIBS)

#========== H A N D Y G  (only for studies linking the ggza_analytic archive) ==========#
# Not appended to LIB_ORDER here -- only a study that links libggza_analytic.a
# needs handyG (ggza_polylog_safe_mod's polylogs), so that study's own
# makefile adds $(HANDYG_LIBS)/$(HANDYG_INCFLAG) itself, plus -no-pie at link
# time (libhandyg.a is not built -fPIE).
HANDYG_LIBDIR  := $(shell $(PKGTOOL) handyg --libdir)
HANDYG_INCFLAG := $(shell $(PKGTOOL) handyg --incflag)
HANDYG_LIBS    := $(shell $(PKGTOOL) handyg --libs)

LDFLAGS = -lc

#========== I N C L U D E  /  L I N K  O R D E R ==========#
INC_ORDER = $(RECOLA_INCFLAG)
LIB_ORDER = $(RECOLA_FLAGS)

#========== F L A G S ==========#
# -w on the F77 rule: the PSGen sources and vsup.f are fixed-form legacy and
# will not compile clean. The .f90 rule stays strict.
#
# There is no avhlib here, so none of the -fcheck workarounds the sibling
# build needs apply; -O is plain.
FFLAGS   = -O -w -std=legacy -J$(MODDIR) -I$(MODDIR)
F95FLAGS = -O -Wall -std=f2008 -J$(MODDIR) -I$(MODDIR)
