# ===========================================================================
# Linux-specific variables for UrQMD-Makefile      Bjoern Baeuchle 21.08.2008
# ===========================================================================

# We need to define a fortran compiler $(FC) and a linker $(LD).
 FC = gfortran
 LD = $(FC)

# define compiler- and linker flags.
FFLAGS	= -O3 -mcmodel=medium
LDFLAGS = -O3 -mcmodel=medium


ifdef DEBUG
 FFLAGS        += -ggdb
 LDFLAGS       += -ggdb
endif
ifdef PROFILER
 FFLAGS        += -pg
 LDFLAGS       += -pg
endif
ifdef WARN
 FFLAGS        +=       -Wall -Wsurprising
endif

SYSTEMFILES = gnuranf.f
