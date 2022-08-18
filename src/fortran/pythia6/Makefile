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


pythia6$(LEXT): _pythia6.pyf $(PYTHIA6).o $(LOGO) $(RND)
	$(F2PY_L) -c $(F2PY_CCONF) --opt="$(OPT)" _pythia6.pyf \
	$(PYTHIA6).o $(RND) $(LOGO) >> pythia6.f2l 2>&1

_pythia6.pyf: $(PYTHIA6).f $(LOGS) $(RNDS)
	$(F2PY) -m _pythia6 -h _pythia6.pyf --overwrite-signature \
	only: pyinit pyexec pytune pylist pyevnt pyevnw pystat pyedit pyhepc \
	pychge pycomp pyk init_rmmard $(LOGF) : $(PYTHIA6).f $(LOGS) $(RNDS) > pythia6.f2l 2>&1
	$(F2PY) _pythia6.pyf


.f.o:
	$(FC) -c $(OPT) $<
