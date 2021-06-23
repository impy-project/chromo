#
#      Makefile for Python-libs of Fortran programs
#
#                                        (A.F. 08/12)
#

DPM_SRCS :=$(wildcard sources/*.f)
# DPM_SRCS := $(filter-out ./user3.0-6.f, $(DPM_SRCS))
DPM_OBJS := $(DPM_SRCS:.f=.o)

# Options for native compilation
#Debug/Safe-mode
# FOPT_USR = -fPIC -fbounds-check -Og -g -Wno-unused \
# 			-Wno-conversion -Wno-compare-reals \
# 			-ffpe-trap=invalid,zero,overflow -std=legacy
# Optimized mode
FOPT_USR = -mtune=native -O3 -Wno-uninitialized -fPIC -std=legacy -Wno-argument-mismatch

OPT += -Wno-argument-mismatch

ifndef OPT
	OPT = $(FOPT_USR)
endif

FOPT = $(OPT) -std=legacy # -Wno-argument-mismatch this flag can suppress dozens of warnings for gfortran 7

DPMJET_FUNCS = pho_event dt_init dt_kkinc idt_icihad dt_xsglau \
pycomp dt_initjs dt_rndmst dt_rndm dt_inucas idt_ipdgha \
pho_init pho_setpar pho_pname pho_pmass pho_setmdl \
pho_setpdf pycomp pho_xsect pho_borncs pho_harmci pho_fitout \
pho_mcini pho_ptcut pytune pho_rregpar pho_sregpar pho_prevnt \
ipho_pdg2id ipho_id2pdg pho_harint pho_harxto pho_harxpt \
pho_setpcomb dt_phoxs dt_xshn dt_flahad $(LOGF)

all: install

.PHONY: clean
clean:
	$(DEL_COMMAND) sources$(PATHSEP)*.o $(CLEANEXT)

.PHONY: distclean
distclean: clean
	rm -f *.pyf

.PHONY: install
install: dpmjet306$(LEXT)
	$(COPY_COMMAND) *$(LEXT) $(LIB_DIR)

dpm3usr: sources$(PATHSEP)user3.0-6.o $(DPM_OBJS)
	$(FC) $(DPM_OBJS) -o dpm3usr

dpmjet306$(LEXT): dpmjet306.pyf $(DPM_OBJS) Makefile
	$(F2PY_L) -c $(F2PY_CCONF) --opt="$(OPT)" \
	dpmjet306.pyf $(DPM_OBJS) $(LOGO) >> dpmjet306.f2l 2>&1

dpmjet306.pyf: $(DPM_SRCS) Makefile
	$(F2PY) -m dpmjet306 -h dpmjet306.pyf --overwrite-signature \
	only: $(DPMJET_FUNCS) : $(DPM_SRCS) $(LOGS) > dpmjet306.f2l 2>&1

sources/%.o: sources/%.f
	$(FC) -c $(FOPT) -o $@ $<
