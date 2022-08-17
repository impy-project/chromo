#
#      Makefile for Python-libs of Fortran programs
#
#                                        (A.F. 08/12)
#

PYTHIA6 = pythia-6.4.28
# LHAPDF_CONFIG = /Users/afedynitch/dev/prog/lhapdf/bin/lhapdf-config
# LHAPDF_LIBS = $(shell $(LHAPDF_CONFIG) --ldflags)

all:  install

.PHONY: clean
clean:
	$(DEL_COMMAND) $(TARGET) $(CLEANEXT)

.PHONY: distclean
distclean: clean
	rm -f *.pyf

.PHONY: install
install: pythia6$(LEXT)
	$(COPY_COMMAND) *$(LEXT) $(LIB_DIR)

#PYTHIA
pythia6$(LEXT): pythia6.pyf $(PYTHIA6).o $(LOGO) $(RND)
	$(F2PY_L) -c $(F2PY_CCONF) --opt="$(OPT)" pythia6.pyf \
	$(PYTHIA6).o $(RND) $(LOGO) >> pythia6.f2l 2>&1

pythia6.pyf: $(PYTHIA6).f $(LOGS) $(RNDS)
	$(F2PY) -m pythia6 -h pythia6.pyf --overwrite-signature \
	only: pyinit pyexec pytune pylist pyevnt pyevnw pystat pyedit pyhepc \
	pychge pycomp pyk init_rmmard $(LOGF) : $(PYTHIA6).f $(LOGS) $(RNDS) > pythia6.f2l 2>&1


.f.o:
	$(FC) -c $(OPT) $<
