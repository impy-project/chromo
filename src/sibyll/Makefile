#
#      Makefile for Python-libs of Fortran programs
#
#                                        (A.F. 08/12)
#

export SIGFUNC = sibyll sibyll_ini sib_sigma_hp sib_sigma_hair int_nuc \
	decsib decpar sibini sib_list isib_pid2pdg isib_pdg2pid pdg_ini \
	$(LOGF)

# Fortran functions for PyInteractionModels
export PIM_FUNCS_F = $(CURDIR)$(PATHSEP)sibyll_init.f $(LOGS)
export PIM_FUNCS_O = $(CURDIR)$(PATHSEP)sibyll_init.o $(LOGO) $(RND)

export SIBMISC_21_F = sib21aux.f 
export SIBMISC_21_O = sib21aux.o $(CURDIR)$(PATHSEP)sibyll_init.o $(LOGO) sib21rnd.o

all: install

.PHONY: clean
clean:
	$(DEL_COMMAND) $(TARGET) *.o *.prj *.chk core *$(LEXT) *.dSYM *.f2l

.PHONY: distclean
distclean: clean
	rm -f *.pyf

.PHONY: install
install: sib21$(LEXT) sib23$(LEXT) sib23c01$(LEXT) sib23c03$(LEXT)
	$(COPY_COMMAND) *$(LEXT) $(LIB_DIR)

sib21$(LEXT): sib21.pyf $(SIBMISC_21_O) $(PIM_FUNCS_O) sibyll_21.o
	$(F2PY_L) -c $(F2PY_CCONF) --opt="$(OPT)" sib21.pyf sibyll_21.o $(SIBMISC_21_O) >> sib21.f2l 2>&1

sib21.pyf: sibyll_21.f $(SIBMISC_21_F) $(PIM_FUNCS_F) Makefile
	$(F2PY) -m sib21 -h sib21.pyf --overwrite-signature only: $(SIGFUNC) sibhep1 : sibyll_21.f $(SIBMISC_21_F) $(PIM_FUNCS_F) > sib21.f2l 2>&1

sib21rnd.o:
	$(FC) -c $(OPT) -cpp -DSIBYLL_SP ..$(PATHSEP)rangen.f -o sib21rnd.o

sib23$(LEXT): sib23.pyf $(PIM_FUNCS_O) sibyll2.3.o 
	$(F2PY_L) -c $(F2PY_CCONF) --opt="$(OPT)" \
	sib23.pyf \
	sibyll2.3.o $(PIM_FUNCS_O) >> sib23.f2l 2>&1

sib23.pyf: sibyll2.3.f $(PIM_FUNCS_F) Makefile
	$(F2PY) -m sib23 -h sib23.pyf --overwrite-signature \
	only: $(SIGFUNC) sibhep3 : sibyll2.3.f $(PIM_FUNCS_F) > sib23.f2l 2>&1

sib23c01$(LEXT): sib23c01.pyf $(PIM_FUNCS_O) sibyll2.3c01.o 
	$(F2PY_L) -c $(F2PY_CCONF) --opt="$(OPT)" \
	sib23c01.pyf sibyll2.3c01.o $(PIM_FUNCS_O) >> sib23c01.f2l 2>&1

sib23c01$(LEXT): sib23c01.pyf $(PIM_FUNCS_O) sibyll2.3c01.o 
	$(F2PY_L) -c $(F2PY_CCONF) --opt="$(OPT)" \
	sib23c01.pyf sibyll2.3c01.o $(PIM_FUNCS_O) >> sib23c01.f2l 2>&1

sib23c01.pyf: sibyll2.3c01.f $(PIM_FUNCS_F) Makefile
	$(F2PY) -m sib23c01 -h sib23c01.pyf --overwrite-signature \
	only: $(SIGFUNC) sibhep3 : sibyll2.3c01.f $(PIM_FUNCS_F) > sib23c01.f2l 2>&1

sib23c03$(LEXT): sib23c03.pyf $(PIM_FUNCS_O) sibyll2.3c03.o 
	$(F2PY_L) -c $(F2PY_CCONF) --opt="$(OPT)" \
	sib23c03.pyf sibyll2.3c03.o $(PIM_FUNCS_O) >> sib23c03.f2l 2>&1

sib23c03$(LEXT): sib23c03.pyf $(PIM_FUNCS_O) sibyll2.3c03.o 
	$(F2PY_L) -c $(F2PY_CCONF) --opt="$(OPT)" \
	sib23c03.pyf sibyll2.3c03.o $(PIM_FUNCS_O) >> sib23c03.f2l 2>&1

sib23c03.pyf: sibyll2.3c03.f $(PIM_FUNCS_F) Makefile
	$(F2PY) -m sib23c03 -h sib23c03.pyf --overwrite-signature \
	only: $(SIGFUNC) sibhep3 : sibyll2.3c03.f $(PIM_FUNCS_F) > sib23c03.f2l 2>&1

.f.o:
	$(FC) -c -cpp $(OPT) $<
