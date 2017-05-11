#
#      Makefile for Python-libs of Fortran programs
#
#                                        (A.F. 08/12)
#

SIB_SRCS :=$(wildcard ./*.f)
SIB_OBJS :=$(SIB_SRCS:.f=.o)

all: install

.PHONY: clean
clean:
	rm -rf *.o *.prj *.chk core *.dSYM *.$(LEXT)

.PHONY: distclean
distclean: clean
	rm -f *.pyf

.PHONY: install
install: sib23_rc4a.$(LEXT)
	cp *.$(LEXT) ../../../

sib23_rc4a.$(LEXT): sib23_rc4a.pyf $(SIB_OBJS) $(PIM_FUNCS_O)
	$(F2PY_L) -c --fcompiler=$(F2PY_C) --opt="$(OPT)" \
	sib23_rc4a.pyf $(SIB_OBJS) $(PIM_FUNCS_O)

sib23_rc4a.pyf: $(SIB_SRCS) $(PIM_FUNCS_F)
	$(F2PY) -m sib23_rc4a -h sib23_rc4a.pyf --overwrite-signature \
	only: $(SIGFUNC) : $(SIB_SRCS) $(PIM_FUNCS_F)

.f.o:
	$(FC) -c $(OPT) $<
