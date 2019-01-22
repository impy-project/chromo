# ===========================================================================
# AIX-specific variables for UrQMD-Makefile        Bjoern Baeuchle 21.08.2008
# ===========================================================================

# copied from old Makefile. Not tested.
FC          = xlf
LD          = xlf
FFLAGS      = -O5 -qextname -qstrict -qipa=partition=large
LDFLAGS     = -O5 -qextname -qstrict -qipa=partition=large
SYSTEMFILES = genranf.f
