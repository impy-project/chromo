#
#      Makefile for Python-libs of Fortran programs
#
#                                        (A.F. 08/12)
#

SRCS= SOPHIA20.f eventgen.f sampling.f inpoutput.f\
      jetset74dp.f
OBJS= SOPHIA20.o eventgen.o sampling.o inpoutput.o\
      jetset74dp.o
      
SFUNCS=eventgen print_event crossection

all:  install

.PHONY: clean
clean:
	rm -rf $(TARGET) *.o *.prj *.chk core *$(LEXT) *.dSYM *.f2l

.PHONY: distclean
distclean: clean
	rm -f *.pyf

.PHONY: install
install: sophia$(LEXT)
	cp *$(LEXT) $(WORK_DIR)/lib

sophia$(LEXT): sophia.pyf $(RND) $(OBJS) $(LOGO)
	$(F2PY_L) -c --fcompiler=$(F2PY_C) --opt="$(OPT)" sophia.pyf \
	$(RND) $(LOGO) $(OBJS) >> sophia.f2l 2>&1

sophia.pyf:
	$(F2PY) -m sophia -h sophia.pyf --overwrite-signature \
	only: $(SFUNCS) $(LOGF) : $(SRCS) $(LOGS) > sophia.f2l 2>&1

.f.o:
	$(FC) -c $(OPT) $<
