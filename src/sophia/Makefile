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

all:  sophia.$(LEXT) install

.PHONY: clean
clean:
	rm -rf $(TARGET) *.o *.prj *.chk core *.$(LEXT) *.dSYM

.PHONY: distclean
distclean: clean
	rm -f *.pyf

.PHONY: install
install: sophia.$(LEXT)
	cp *.$(LEXT) $(WORK_DIR)

sophia.$(LEXT): sophia.pyf $(RND) $(OBJS) $(LOGO)
	$(F2PY_L) -c --fcompiler=$(F2PY_C) --opt="$(OPT)" sophia.pyf \
	$(RND) $(LOGO) $(OBJS)

sophia.pyf:
	$(F2PY) -m sophia -h sophia.pyf --overwrite-signature \
	only: $(SFUNCS) $(LOGF) : $(SRCS) $(LOGS)

.f.o:
	$(FC) -c $(OPT) $<
