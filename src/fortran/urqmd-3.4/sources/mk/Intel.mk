# ===========================================================================
# Intel-specific variables for UrQMD-Makefile      Bjoern Baeuchle 21.08.2008
# ===========================================================================

FC          = ifort
LD          = $(FC)
FFLAGS      = -no-ipo -O3
LDFLAGS     = $(FFLAGS)
SYSTEMFILES = intranf.f erf.f

TYPE := $(TYPE).Intel
