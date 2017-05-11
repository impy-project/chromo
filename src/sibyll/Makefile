#
#      Makefile for Python-libs of Fortran programs
#
#                                        (A.F. 08/12)
#

export SIGFUNC = sibyll sibyll_ini sib_sigma_hp sib_sigma_hair \
	decsib decpar sibini sib_list isib_pid2pdg isib_pdg2pid pdg_ini \
	$(LOGF)

# Fortran functions for PyInteractionModels
export PIM_FUNCS_F = $(CURDIR)/sibyll_init.f $(LOGS)
export PIM_FUNCS_O = $(CURDIR)/sibyll_init.o $(LOGO) $(RND)

export SIBMISC_21_F = sib21aux.f 
export SIBMISC_21_O = sib21aux.o 

SIB_SUBDIRS = sibyll23rc3a  sibyll23rc4a 

#all: install $(SIB_SUBDIRS)

all: install

.PHONY: subdirs $(SIB_SUBDIRS)
     
subdirs: $(SIB_SUBDIRS)
     
$(SIB_SUBDIRS): $(PIM_FUNCS_O)
	$(MAKE) -C $@

.PHONY: clean
clean:
	rm -rf $(TARGET) *.o *.prj *.chk core *.$(LEXT) *.dSYM
	for d in $(SIB_SUBDIRS); \
	do \
	$(MAKE) --directory=$$d clean; \
	done

.PHONY: distclean
distclean: clean
	rm -f *.pyf
	for d in $(SIB_SUBDIRS); \
	do \
	$(MAKE) --directory=$$d distclean; \
	done

.PHONY: install
install: sib21.$(LEXT) sib23_rc1.$(LEXT) sib23.$(LEXT) sib232.$(LEXT) \
	sib233.$(LEXT) sib234.$(LEXT) sib235.$(LEXT)
	cp *.$(LEXT) $(WORK_DIR)

sib21.$(LEXT): sib21.pyf $(SIBMISC_21_O) $(PIM_FUNCS_O) sibyll_21.o
	$(F2PY_L) -c --fcompiler=$(F2PY_C) --opt="$(OPT)" sib21.pyf \
	sibyll_21.o $(SIBMISC_21_O) $(PIM_FUNCS_O)

sib21.pyf: sibyll_21.f $(SIBMISC_21_F) $(PIM_FUNCS_F)
	$(F2PY) -m sib21 -h sib21.pyf --overwrite-signature \
	only: $(SIGFUNC) : sibyll_21.f $(SIBMISC_21_F) $(PIM_FUNCS_F)

sib21rnd.o:
	$(FC) -c $(OPT) -xf77-cpp-input -DSIBYLL_SP ../../rangen.f -o sib21rnd.o

sib23_rc1.$(LEXT): sib23_rc1.pyf $(SIBMISC_23_O) $(PIM_FUNCS_O) sibyll_23_rc1.o 
	$(F2PY_L) -c --fcompiler=$(F2PY_C) --opt="$(OPT)" \
	sib23_rc1.pyf sibyll_23_rc1.o $(SIBMISC_23_O) $(PIM_FUNCS_O)

sib23_rc1.pyf: sibyll_23_rc1.f $(SIBMISC_23_F) $(PIM_FUNCS_F) 
	$(F2PY) -m sib23_rc1 -h sib23_rc1.pyf --overwrite-signature \
	only: $(SIGFUNC) : sibyll_23_rc1.f $(SIBMISC_23_F) $(PIM_FUNCS_F)

sib23.$(LEXT): sib23.pyf $(PIM_FUNCS_O) sibyll2.3.o 
	$(F2PY_L) -c --fcompiler=$(F2PY_C) --opt="$(OPT)" \
	sib23.pyf \
	sibyll2.3.o $(PIM_FUNCS_O)

sib23.pyf: sibyll2.3.f $(PIM_FUNCS_F)
	$(F2PY) -m sib23 -h sib23.pyf --overwrite-signature \
	only: $(SIGFUNC) : sibyll2.3.f $(PIM_FUNCS_F)

sib232.$(LEXT): sib232.pyf $(PIM_FUNCS_O) sibyll2.3.2.o 
	$(F2PY_L) -c --fcompiler=$(F2PY_C) --opt="$(OPT)" \
	sib232.pyf \
	sibyll2.3.2.o $(PIM_FUNCS_O)

sib232.pyf: sibyll2.3.2.f $(PIM_FUNCS_F)
	$(F2PY) -m sib232 -h sib232.pyf --overwrite-signature \
	only: $(SIGFUNC) : sibyll2.3.2.f $(PIM_FUNCS_F)

sib233.$(LEXT): sib233.pyf $(PIM_FUNCS_O) sibyll2.3.3.o 
	$(F2PY_L) -c --fcompiler=$(F2PY_C) --opt="$(OPT)" \
	sib233.pyf \
	sibyll2.3.3.o $(PIM_FUNCS_O)

sib233.pyf: sibyll2.3.3.f $(PIM_FUNCS_F)
	$(F2PY) -m sib233 -h sib233.pyf --overwrite-signature \
	only: $(SIGFUNC) : sibyll2.3.3.f $(PIM_FUNCS_F)

sib234.$(LEXT): sib234.pyf $(PIM_FUNCS_O) sibyll2.3.4.o 
	$(F2PY_L) -c --fcompiler=$(F2PY_C) --opt="$(OPT)" \
	sib234.pyf \
	sibyll2.3.4.o $(PIM_FUNCS_O)

sib234.pyf: sibyll2.3.4.f $(PIM_FUNCS_F)
	$(F2PY) -m sib234 -h sib234.pyf --overwrite-signature \
	only: $(SIGFUNC) : sibyll2.3.4.f $(PIM_FUNCS_F)

sib235.$(LEXT): sib235.pyf $(PIM_FUNCS_O) sibyll2.3.5.o 
	$(F2PY_L) -c --fcompiler=$(F2PY_C) --opt="$(OPT)" \
	sib235.pyf \
	sibyll2.3.5.o $(PIM_FUNCS_O)

sib235.pyf: sibyll2.3.5.f $(PIM_FUNCS_F)
	$(F2PY) -m sib235 -h sib235.pyf --overwrite-signature \
	only: $(SIGFUNC) : sibyll2.3.5.f $(PIM_FUNCS_F)

.f.o:
	$(FC) -c $(OPT) $<
