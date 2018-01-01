#
#      Makefile for Python-libs of Fortran programs
#
#                                        (A.F. 08/12)
#

SRCS=  CFmax.f90 quadri.f90 cornelius.f90  bessel.f 1fluid.f     colltab.f  comstr.f   detbal.f     genranf.f  hepcmp.f intranf.f   make22.f   output.f      scatter.f    upmerge.f addpart.f    blockres.f  coload.f   comwid.f   dwidth.f     getmass.f  hepnam.f  iso.f       newpart.f  paulibl.f     siglookup.f  urqmd.f alpharanf.f  boxinc.f    comnorm.f  dectim.f   erf.f        getspin.f  init.f ityp2pdg.f  numrec.f   proppot.f     string.f     whichres.f angdis.f     boxprg.f    comres.f   defs.f     error.f      gnuranf.f  input.f   itypdata.f  options.f  pythia6409.f  tabinit.f anndec.f     cascinit.f  coms.f     delpart.f  freezeout.f  hepchg.f   inputs.f  jdecay2.f  outcom.f   saveinfo.f    uhmerge.f urqinit.f


OBJS= CFmax.o  quadri.o cornelius.o bessel.o  1fluid.o        delpart.o  getmass.o  hepcmp.o  iso.o       numrec.o   pythia6409.o  siglookup.o  upmerge.o addpart.o  blockres.o  coload.o     detbal.o  getspin.o  hepnam.o  ityp2pdg.o  output.o   string.o     urqmd.o angdis.o boxprg.o    dwidth.o   gnuranf.o  init.o    jdecay2.o   paulibl.o  saveinfo.o    tabinit.o    whichres.o anndec.o   cascinit.o  dectim.o     error.o    hepchg.o   input.o make22.o    proppot.o  scatter.o     uhmerge.o urqinit.o
      
SFUNCS=urqmd init uinit set0 params uounit strini loginit loadwtab norm_init output cascinit urqini

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

urqmd.$(LEXT): urqmd.pyf $(OBJS) $(LOGO)
	$(F2PY_L) -c --fcompiler=$(F2PY_C) --opt="$(OPT)" urqmd.pyf \
	$(LOGO) $(OBJS)

urqmd.pyf:
	$(F2PY) -m urqmd -h urqmd.pyf --overwrite-signature \
	only: $(SFUNCS) $(LOGF) : $(SRCS) $(LOGS)

%.o: %.f
	$(FC) -c $(OPT) $<

%.o: %.f90
	$(FC) -c $(OPTF90) $<
