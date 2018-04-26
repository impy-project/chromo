#
#      Makefile for Python-libs of Fortran programs
#
#                                        (A.F. 08/12)
#

EPOS_SRCS :=$(wildcard ./*.f)
EPOS_SRCS := $(filter-out ./epos_example.f, $(EPOS_SRCS))
EPOS_OBJS := $(EPOS_SRCS:.f=.o)
EX_OBJS := $(filter-out ./epos_interface.o, $(EPOS_OBJS))
EPOS_FUNC = aaset ainit aepos afinal hepmcstore \
			getcharge idtrafo initializeepos initeposevt \
			setstable xsection $(LOGF)

all: install

.PHONY: clean
clean:
	rm -rf *.o *.prj *.chk core *.dSYM *.$(LEXT) *.f2l

.PHONY: distclean
distclean: clean
	rm -f *.pyf

.PHONY: install
install: eposlhc.$(LEXT)
	cp *.$(LEXT) ../../

example: ./epos_example.f $(EX_OBJS)
	gfortran ./epos_example.f $(EX_OBJS) -o epos_example.exe

eposlhc.$(LEXT): eposlhc.pyf $(EPOS_OBJS)
	$(F2PY_L) -c --fcompiler=$(F2PY_C) --opt="$(OPT)" \
	eposlhc.pyf $(EPOS_OBJS) $(LOGO) >> epos-lhc.f2l 2>&1

eposlhc.pyf: epos-bas-lhc.f epos_interface.f Makefile
	$(F2PY) -m eposlhc -h eposlhc.pyf --overwrite-signature \
	only: $(EPOS_FUNC) : epos-bas-lhc.f epos-ids-lhc.f \
	epos_interface.f $(LOGS) > epos-lhc.f2l 2>&1

.f.o:
	$(FC) -c $(OPT) $<
