#
#      Makefile for Python-libs of Fortran programs
#
#                                        (A.F. 08/12)
#

CVendor = "GNU"
Config = "Release"

WORK_DIR = $(CURDIR)
RND = $(CURDIR)/src/rangen.o
# Shared library suffix
LEXT = so

# For f2py
LOGF = pho_openlogfile pho_closelogfile

#######################################################################
#
#   compiler
#
#######################################################################

ifeq ($(CVendor),"GNU")
	#  GNU
	FC = gfortran
	F2PY_C = gnu95
else
	#  Intel
	FC = ifort
	F2PY_C = intelem
endif

#######################################################################
#
#   compiler options for different platforms
#
#######################################################################
ifeq ($(CVendor),"GNU")
	ifeq ($(Config),"Debug")
		# GNU Debug
		OPT = -C -fPIC -Wall -fbounds-check -O0 -g \
			  -ffpe-trap=invalid,zero,overflow -Wuninitialized
		#OPT = -C -fPIC -Wall -Wno-uninitialized -Wno-unused-variable -O3 -g -ffpe-trap=invalid,zero,overflow
	else
		# GNU Release
		#OPT = -C -O0 -fPIC
		OPT = -C -ftree-vectorize -O3 -Wno-uninitialized -fPIC
	endif
else
	ifeq ($(Config),"Debug")
	# Intel Debug (-gen-interfaces -warn interfaces)
		OPT = -C -check bounds -O0 -g -check pointer -fpe0 -traceback
	else
		# Intel Release
		OPT = -C -fast -fpe0
	endif
endif

#######################################################################
#
#   F2PY
#
#######################################################################
#general version for signature file extraction and linking
F2PY = f2py --quiet
#additional flags for linker
F2PY_L = $(F2PY) -D_NPY_1_7_DEPRECATED_API_H
#######################################################################
#
#   Targets
#
#######################################################################

export

all: src

.PHONY: src
src:
	$(MAKE) -C $@

.PHONY: clean
clean:
	rm -f $(TARGET) *.o *.prj *.chk core  *.pyf
	$(MAKE) --directory=src clean
#*.$(LEXT)
.PHONY: distclean
distclean:
	rm -rf $(TARGET) *.o *.prj *.chk core *.$(LEXT) *.pyf *.dSYM
	$(MAKE) --directory=src distclean

.f.o:
	$(FC) -c $(OPT) -xf77-cpp-input $<
