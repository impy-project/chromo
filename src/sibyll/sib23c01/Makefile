#
#      Makefile for SIBYLL   (R.Engel)
#
# options:
#
SYSTEM = $(shell uname -s)

# default settings
OPT = -c -O0 -g -pg
DEB = 

# Linux Intel settings
ifeq ("$(findstring Linux, $(SYSTEM))", "Linux")
# g77
#  OPT = -c -Wall -fno-automatic -fbounds-check -malign-double \
        -finit-local-zeros -march=pentium -O0 -g
# gfortran
#  OPT = -c -Wall -fno-automatic -fbounds-check -Wuninitialized -O0 -g -pg
#  OPT = -Wall -Wextra -std=legacy -O0 -ggdb3 -fno-inline -fbounds-check -fno-second-underscore -ffpe-trap=invalid,zero,overflow -finit-local-zero
  OPT = -Wall -Wextra -std=legacy -fstack-protector -O3 -ggdb3 -fno-inline -fbounds-check -fno-second-underscore -ffpe-trap=invalid,zero,overflow
  DEB = 
endif

# Mac OSX PPC settings
ifeq ("$(findstring Darwin, $(SYSTEM))", "Darwin")
  OPT = -c -Wall -fbounds-check -O0 -g
  DEB =
endif

#
# Libraries:
#
# (to link on HP libraries, use fort77 instead of f77 for LO)
# LIB = -L /usr/local/cern/96a/lib -lmathlib -lpacklib -lkernlib
#LIB =
#
# fortran compiler, linker/loader:
#
FC = gfortran
LO = gfortran
#
# C++ compiler
CCP = g++
CCPFLAGS = -g -I ./


OBJ = sib_cms_f.o rhisto_f.o

UTILDIR = ./utils

#SOURCEDIR = ./src
SOURCES = sibyll2.3c.f $(UTILDIR)/rndm_dbl.f
TESTDIR = ./test
#INCLUDES = $(wildcard $(SOURCEDIR)/*.inc)
BUILDDIR = .
#OBJECTS = $(patsubst $(SOURCEDIR)/%.f, $(BUILDDIR)/%_f.o, $(SOURCES))
OBJECTS = sibyll2.3c_f.o rndm_dbl_f.o #py_link_f.o pythia_decay_f.o

PYMODULE_FUNCTIONS = sibyll sibyll_ini ini_event decsib sib_sig string_frag_4flv sib_list sib_sigma_hp dec_ini prnt_prtn_stck sib_sigma_hair sib_sigma_hnuc fpni int_len_ini sig_rpp2014 sigma_pp sigma_kp sigma_pip sigma_ini fsighair sib_partpr sib_diff dec_ini rnd_ini signuc_ini nuc_nuc_ini


help:
	@echo ""
	@echo "  SIBYLL compile Makefile "
	@echo " "
	@echo "  Select components you want to compile "
	@echo "  Options are: 'sib_cms'     : FORTRAN example"
	@echo "               'sib_cms_cc'  : C++ example"
	@echo "               'sib23<X>crs' : compile Python module"
	@echo "                               where <X> stands for the current sub-version"
	@echo "                               crs is added to the name to indicate that no config. input is required (CORSIKA config)"
	@echo "                               the random number generator has to be initialized though"
	@echo "                               python example can then be run with 'python sib_cms.py'"
	@echo "  "
	@echo "  test programs:"	
	@echo "               'test_sibyll'      : print 10 events"
	@echo "               'test_partpr'      : print the particle ids in Sibyll"
	@echo "               'test_decpr'       : print all decay channels in Sibyll"
	@echo "                  ... many more. check ./test/ "
	@echo " "


all: help

clean:
	rm -vf *.o sib_cms sib_cms_cc test_* *.prj *.so *.pyf fort.7 *__genmod.f90 *__genmod.mod


test_%:  test_%_f.o $(OBJECTS) $(DEB) Makefile
	$(LO)  -o $@ $(OBJECTS) $@_f.o $(DEB) $(LIB)

sib_cms:  $(OBJECTS) $(OBJ) $(DEB) Makefile
	@echo "   linking program ..."
	$(LO)  -o sib_cms $(OBJECTS) $(OBJ) $(DEB) $(LIB)
	@echo "   run program with 'sib_cms' "

sib_cms_cc: $(OBJECTS) sib_cms_cc.o
	$(CCP) $(CCPFLAGS) -o sib_cms_cc $+ -lgfortran

sib23c01crs.pyf: $(SOURCES) $(INCLUDES)
	f2py -m sib23c01crs -h sib23c01crs.pyf --overwrite-signature only: $(PYMODULE_FUNCTIONS) : $(SOURCES)

sib23c01crs:  $(SOURCES) sib23c01crs.pyf
	@echo "wrapping SIBYLL 2.3c01crs ..."
	f2py -c --fcompiler=gnu95 --opt="$(OPT)"  $+

%_c.o: $(SOURCEDIR)/%.c
	gcc -c -o $@ $<

%_cc.o: %.cc
	$(CCP) $(CCPFLAGS) -c -o $@ $<

test_%_f.o:	$(TESTDIR)/test_%.f
	$(FC) -c $(OPT) -o $@ $<

%_f.o:	$(SOURCEDIR)/%.f
	$(FC) -c $(OPT) -o $@ $<

%_f.o:	$(UTILDIR)/%.f
	$(FC) -c $(OPT) -o $@ $<

%_f.o:	%.f
	$(FC) -c $(OPT) -o $@ $<

sib_cms.o: sib_cms.f
	$(FC) -c $(OPT) -o $@ $<


# tell make what a project file suffix is
.SUFFIXES: .prj

    # these options suppress global checks.
NOGLOBAL=-usage=no-ext-undefined

    # tell make how to create a .prj file from a .f file
.f.prj:
	ftnchek -project $(NOGLOBAL) -library $<

    # set up macro PRJS containing project filenames
PRJS= $(OBJ:.o=.prj)

    # "make check" will check everything that has been changed.
check: $(PRJS)
	ftnchek $(PRJS)
