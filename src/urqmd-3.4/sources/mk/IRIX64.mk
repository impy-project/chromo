# ===========================================================================
# IRIX64-specific variables for UrQMD-Makefile     Bjoern Baeuchle 21.08.2008
# ===========================================================================

# copied from old Makefile. Not tested.
FC          = f77
LD          = f77
FFLAGS      = -mips4 -64 -r10000  
LDFLAGS     = -mips4 -64 -r10000
SYSTEMFILES = genranf.f
