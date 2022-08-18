# ===========================================================================
# OSF1-specific variables for UrQMD-Makefile       Bjoern Baeuchle 21.08.2008
# ===========================================================================

# copied from old Makefile. Not tested.
FC          = f77
LD          = f77
FFLAGS      = -C -align dcommons  
LDFLAGS     = -C -align dcommons 
SYSTEMFILES = alpharanf.f
