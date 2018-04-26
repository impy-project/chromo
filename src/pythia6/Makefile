#
#      Makefile for Python-libs of Fortran programs
#
#                                        (A.F. 08/12)
#

PYTHIA6 = pythia-6.4.28
# LHAPDF_CONFIG = /Users/afedynitch/dev/prog/lhapdf/bin/lhapdf-config
# LHAPDF_LIBS = $(shell $(LHAPDF_CONFIG) --ldflags)

all:  pythia6.$(LEXT) install

.PHONY: clean
clean:
	rm -rf $(TARGET) *.o *.prj *.chk core *.$(LEXT) *.dSYM *.f2l

.PHONY: distclean
distclean: clean
	rm -f *.pyf

.PHONY: install
install: pythia6.$(LEXT)
	cp *.$(LEXT) $(WORK_DIR)

#PYTHIA
pythia6.$(LEXT): pythia6.pyf $(PYTHIA6).o $(LOGO) 
	$(F2PY_L) -c --fcompiler=$(F2PY_C) --opt="$(OPT)" pythia6.pyf \
	$(PYTHIA6).o $(LOGO) >> pythia6.f2l 2>&1

pythia6.pyf: $(PYTHIA6).f $(LOGS)
	$(F2PY) -m pythia6 -h pythia6.pyf --overwrite-signature \
	only: pyinit pyexec pytune pylist pyevnt pyevnw pystat pyedit \
	pychge pycomp pyk $(LOGF) : $(PYTHIA6).f $(LOGS) > pythia6.f2l 2>&1


.f.o:
	$(FC) -c $(OPT) $<
