#
#      Makefile for Python-libs of Fortran programs
#
#                                        (A.F. 08/12)
#

DPMJET_O = dpmjet253c.o dpmjet254.o dpmjet255.o dpmjet256.o testdpmjetII.o dummy.o
DPM2FUNCS = dpjsig dpmjin dpmjlk dpmjet2_event mcihad mpdgha \
	dpmjet2_lt2lab pycomp pylist dechkk xshn siinel

all: dpmjet2$(LEXT) install

.PHONY: clean
clean:
	rm -rf $(TARGET) *.o *.prj *.chk core *$(LEXT) *.dSYM *.f2l

.PHONY: distclean
distclean: clean
	rm -f *.pyf

.PHONY: install
install: dpmjet2$(LEXT)
	cp *$(LEXT) $(WORK_DIR)

dpmjet2$(LEXT): dpmjet2.pyf $(LOGO) $(RND) $(DPMJET_O) dpmjet_iface.o
	@echo "RND is " $(RND) $(WORK_DIR)
	$(F2PY_L) -c --fcompiler=$(F2PY_C) --opt="$(OPT)" dpmjet2.pyf \
	$(LOGS) $(DPMJET_O) $(RND) dpmjet_iface.o >> dpmjet2.f2l 2>&1

dpmjet2.pyf: dpmjet_iface.f90 dpmjet253c.f dpmjet255.f Makefile
	$(F2PY) -m dpmjet2 -h dpmjet2.pyf --overwrite-signature \
	only: $(DPM2FUNCS) : \
	dpmjet_iface.f90 dpmjet253c.f dpmjet255.f > dpmjet2.f2l 2>&1
#dpmjet255.f dpmjet256.f

dpmjettest: $(DPMJET_O) $(RND) $(LOGO) dpmjet_iface.o
	gfortran -O0 -g -o dpmjettest $(DPMJET_O) $(RND) $(LOGO) dpmjet_iface.o

.f.o:
	$(FC) -c $(OPT) $<

%.o: %.f90
	$(FC) $(OPT) -c -o $@ $<   
