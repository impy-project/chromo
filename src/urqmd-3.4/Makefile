#
#      Makefile for Python-libs of Fortran programs
#
#                                        (A.F. 08/12)
#

URQMD_FILES = \
1fluid bessel delpart getmass hepcmp iso numrec pythia6409 \
siglookup upmerge addpart blockres coload detbal getspin \
hepnam ityp2pdg output string urqmd angdis boxprg dwidth \
gnuranf init jdecay2 paulibl saveinfo tabinit whichres anndec \
cascinit dectim error hepchg input make22 proppot scatter \
uhmerge urqinit

#$(filter-out ${NOT_NEEDED:=.f}, $(wildcard ./*.f))

URQMD_SRCS := ${URQMD_FILES:=.f}
URQMD_90SRCS := CFmax.f90 quadri.f90 cornelius.f90
URQMD_OBJS := $(URQMD_SRCS:.f=.o) $(URQMD_90SRCS:.f90=.o)



SFUNCS=urqmd init uinit set0 params uounit strini loginit loadwtab norm_init output cascinit urqini pdgid pdgid partname

all:  urqmd.$(LEXT) install

.PHONY: clean
clean:
	rm -rf $(TARGET) *.o *.prj *.chk core *.$(LEXT) *.dSYM

.PHONY: distclean
distclean: clean
	rm -f *.pyf

.PHONY: install
install: urqmd.$(LEXT)
	cp *.$(LEXT) $(WORK_DIR)

urqmd.$(LEXT): urqmd.pyf $(URQMD_OBJS) $(LOGO)
	$(F2PY_L) -c --fcompiler=$(F2PY_C) --opt="$(OPT)" urqmd.pyf \
	$(LOGO) $(URQMD_OBJS)

urqmd.pyf:
	$(F2PY) -m urqmd -h urqmd.pyf --overwrite-signature \
	only: $(SFUNCS) $(LOGF) : $(URQMD_SRCS) $(URQMD_90SRCS) $(LOGS)

%.o: %.f
	$(FC) -c $(OPT) $<

%.o: %.f90
	$(FC) -c $(OPTF90) $<
