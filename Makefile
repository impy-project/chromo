#
#      Makefile for Python-libs of Fortran programs
#
#                                        (A.F. 08/12)
#

CVendor = "GNU"
Config = "Release"

WORK_DIR = $(CURDIR)
LIB_DIR?=$(WORK_DIR)/lib
RND = $(CURDIR)/src/rangen.o
# Shared library suffix
LEXT?=$(shell python -c 'import sysconfig; print(sysconfig.get_config_var("EXT_SUFFIX"))')

# For f2py
LOGF = impy_openlogfile impy_closelogfile
export NPY_DISTUTILS_APPEND_FLAGS=0


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
		OPT = -fPIC -Wall -fbounds-check -O0 -g \
			  -ffpe-trap=invalid,zero,overflow -Wuninitialized
		OPTF90 = -fPIC -Wall -fbounds-check -O0 -g \
			  -ffpe-trap=invalid,zero,overflow -Wuninitialized \
			  -fno-second-underscore
		#OPT = -fPIC -Wall -Wno-uninitialized -Wno-unused-variable -O3 -g -ffpe-trap=invalid,zero,overflow
	else
		# GNU Release
		#OPT = -O0 -fPIC
		OPT = -O3 -Wno-uninitialized -fPIC -mtune=native
		OPTF90 = -O3 -Wno-uninitialized -fPIC -fno-second-underscore 
	endif
else
	ifeq ($(Config),"Debug")
	# Intel Debug (-gen-interfaces -warn interfaces)
		OPT = -check bounds -O0 -g -check pointer -fpe0 -traceback
		OPTF90 = -check bounds -O0 -g -check pointer -fpe0 -traceback \ 
			 -cpp -ffree-form -Wobsolescent -fno-second-underscore
	else
		# Intel Release
		OPTF90 = -fast -fpe0 \ 
		      -cpp -ffree-form -Wobsolescent -fno-second-underscore
		OPT = -fast -fpe0	
	endif
endif

#######################################################################
#
#   F2PY
#
#######################################################################
#general version for signature file extraction and linking
F2PY = python -m numpy.f2py --quiet
#additional flags for linker
F2PY_L = $(F2PY)

#######################################################################
#
#   Targets
#
#######################################################################

export

all: odir src

.PHONY: odir
odir:
	mkdir -p $(LIB_DIR)

.PHONY: src
src:
	$(MAKE) -C $@

.PHONY: clean
clean:
	rm -f $(TARGET) *.o *.prj *.chk core  *.pyf *$(LEXT)
	$(MAKE) --directory=src clean

.PHONY: distclean
distclean:
	rm -rf $(TARGET) *.o *.prj *.chk core *$(LEXT) *.pyf *.dSYM
	$(MAKE) --directory=src distclean

%.o: %.f
	$(FC) -c $(OPT) -cpp $<

%.o: %.f90
	$(FC) -c $(OPTF90) -cpp $<
