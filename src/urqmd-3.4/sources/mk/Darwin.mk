# ===========================================================================
# Darwin/Mac specific variables for UrQMD-Makefile Bjoern Baeuchle 20.08.2008
# ===========================================================================

# copied from old Makefile. Not tested.
FC          = gfortran
LD          = gfortran
FFLAGS      = -O 
LDFLAGS     = -O

ifdef WARN
 FFLAGS    += -Wall -Wsurprising
endif

SYSTEMFILES = gnuranf.f
