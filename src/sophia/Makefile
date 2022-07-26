#
#      Makefile for Python-libs of Fortran programs
#
#                                        (A.F. 08/12)
#

SRCS= SOPHIA20.f eventgen.f sampling.f inpoutput.f\
      jetset74dp.f utility.f
OBJS= $(subst .f,.o,$(SRCS))
LIB_NAME = sophia
      
SFUNCS=eventgen print_event crossection initial\
      icon_pdg_sib init_rmmard toevt

all:  install

.PHONY: clean
clean:
	$(DEL_COMMAND) $(TARGET) $(CLEANEXT)

.PHONY: distclean
distclean: clean
	rm -f *.pyf

.PHONY: remove_generated_files
remove_generated_files:
	rm -f *.pyf *.f2l *.so

install: remove_generated_files
$(OBJS): remove_generated_files

$(LIB_DIR)/%$(LEXT): %$(LEXT)
	$(COPY_COMMAND) $< $@	


.PHONY: install
install: $(LIB_DIR)/$(LIB_NAME)$(LEXT)

sophia$(LEXT): sophia.pyf $(RND) $(OBJS) $(LOGO) $(SRCS)
	$(F2PY_L) -c $(F2PY_CCONF) --opt="$(OPT)" sophia.pyf \
	$(RND) $(LOGO) $(OBJS) >> sophia.f2l 2>&1

sophia.pyf: $(RNDS) $(SRCS) $(LOGS) $(OBJS) Makefile
	$(F2PY) -m sophia -h sophia.pyf --overwrite-signature \
	only: $(SFUNCS) $(LOGF) : $(SRCS) $(LOGS) $(RNDS) > sophia.f2l 2>&1

.f.o:
	$(FC) -c $(OPT) $<
	