#
#      Makefile for Python-libs of Fortran programs
#
#                                        (A.F. 08/12)
#
# If f2py fails, try with --verbose output. It is very likely
# that inline comments starting with ! after an expression cause problems


EPOS_SRCS :=$(wildcard sources/*.f)
EPOS_SRCS := $(filter-out sources/epos_example.f, $(EPOS_SRCS))
EPOS_OBJS := $(EPOS_SRCS:.f=.o)
EX_OBJS := $(filter-out sources/epos_interface.o, $(EPOS_OBJS))
EPOS_FUNC = aaset ainit aepos afinal hepmcstore \
			getcharge idtrafo initializeepos initeposevt \
			setstable setunstable xsection $(LOGF)

all: install

.PHONY: clean
clean:
	$(DEL_COMMAND) sources$(PATHSEP)*.o $(CLEANEXT)

.PHONY: distclean
distclean: clean
	$(DEL_COMMAND) *.pyf $(COPY_DUMP)

.PHONY: install
install: eposlhc$(LEXT)
	$(COPY_COMMAND) *$(LEXT) $(LIB_DIR)

example: sources/epos_example.f $(EX_OBJS)
	$(FC) sources$(PATHSEP)epos_example.f $(EX_OBJS) -o epos_example.exe

eposlhc$(LEXT): _eposlhc.pyf $(EPOS_OBJS)
	$(F2PY_L) -c $(F2PY_CCONF) --opt="$(OPT)" _eposlhc.pyf $(EPOS_OBJS) $(LOGO) >> epos-lhc.f2l 2>&1

_eposlhc.pyf: sources/epos-bas-lhc.f sources/epos_interface.f Makefile
	$(F2PY) -m _eposlhc -h _eposlhc.pyf --overwrite-signature only: $(EPOS_FUNC) : sources$(PATHSEP)epos-bas-lhc.f sources$(PATHSEP)epos-ids-lhc.f sources$(PATHSEP)epos_interface.f $(LOGS) > epos-lhc.f2l 2>&1
	$(F2PY) _eposlhc.pyf

sources/%.o: sources/%.f
	$(FC) -c $(OPT) -o $@ $<
