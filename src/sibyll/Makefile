#
#      Makefile for Python-libs of Fortran programs
#
#                                        (A.F. 08/12)
#

export SIGFUNC = sibyll sibyll_ini sib_sigma_hp sib_sigma_hair int_nuc \
	decsib decpar sibini sib_list isib_pid2pdg isib_pdg2pid pdg_ini \
	$(LOGF)

# Fortran functions for PyInteractionModels
export PIM_FUNCS_F = $(CURDIR)/sibyll_init.f $(LOGS)
export PIM_FUNCS_O = $(CURDIR)/sibyll_init.o $(LOGO) $(RND)

export SIBMISC_21_F = sib21aux.f 
export SIBMISC_21_O = sib21aux.o $(CURDIR)/sibyll_init.o $(LOGO) sib21rnd.o

#all: install $(SIB_SUBDIRS)

all: install

.PHONY: clean
clean:
	rm -rf $(TARGET) *.o *.prj *.chk core *.$(LEXT) *.dSYM *.f2l

.PHONY: distclean
distclean: clean
	rm -f *.pyf

.PHONY: install
install: sib21.$(LEXT) sib23.$(LEXT) sib23c.$(LEXT)
	cp *.$(LEXT) $(WORK_DIR)

sib21.$(LEXT): sib21.pyf $(SIBMISC_21_O) $(PIM_FUNCS_O) sibyll_21.o
	$(F2PY_L) -c --fcompiler=$(F2PY_C) --opt="$(OPT)" sib21.pyf \
	sibyll_21.o $(SIBMISC_21_O) >> sib21.f2l 2>&1

sib21.pyf: sibyll_21.f $(SIBMISC_21_F) $(PIM_FUNCS_F) Makefile
	$(F2PY) -m sib21 -h sib21.pyf --overwrite-signature \
	only: $(SIGFUNC) : sibyll_21.f $(SIBMISC_21_F) \
	$(PIM_FUNCS_F) > sib21.f2l 2>&1

sib21rnd.o:
	$(FC) -c $(OPT) -cpp -DSIBYLL_SP ../rangen.f -o sib21rnd.o

sib23.$(LEXT): sib23.pyf $(PIM_FUNCS_O) sibyll2.3.o 
	$(F2PY_L) -c --fcompiler=$(F2PY_C) --opt="$(OPT)" \
	sib23.pyf \
	sibyll2.3.o $(PIM_FUNCS_O) >> sib23.f2l 2>&1

sib23.pyf: sibyll2.3.f $(PIM_FUNCS_F) Makefile
	$(F2PY) -m sib23 -h sib23.pyf --overwrite-signature \
	only: $(SIGFUNC) : sibyll2.3.f $(PIM_FUNCS_F) > sib23.f2l 2>&1

sib23c.$(LEXT): sib23c.pyf $(PIM_FUNCS_O) sibyll2.3c.o 
	$(F2PY_L) -c --fcompiler=$(F2PY_C) --opt="$(OPT)" \
	sib23c.pyf sibyll2.3c.o $(PIM_FUNCS_O) >> sib23c.f2l 2>&1

sib23c.pyf: sibyll2.3c.f $(PIM_FUNCS_F) Makefile
	$(F2PY) -m sib23c -h sib23c.pyf --overwrite-signature \
	only: $(SIGFUNC) : sibyll2.3c.f $(PIM_FUNCS_F) > sib23c.f2l 2>&1

.f.o:
	$(FC) -c -cpp $(OPT) $<
