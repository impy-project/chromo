#
#      Makefile for Python-libs of Fortran programs
#
#                                        (A.F. 08/12)
#

all:  qgsII03.$(LEXT) qgsII04.$(LEXT) qgs01.$(LEXT) install

.PHONY: clean
clean:
	rm -rf $(TARGET) *.o *.prj *.chk core *.$(LEXT) *.dSYM

.PHONY: distclean
distclean: clean
	rm -f *.pyf

.PHONY: install
install: qgsII03.$(LEXT) qgsII04.$(LEXT) qgs01.$(LEXT)
	cp *.$(LEXT) $(WORK_DIR)

#QGSJET
qgsII03.$(LEXT): qgsII-03.pyf $(RND) qgsjet-II-03.o qgsjet-II_init.o
	$(F2PY_L) -c --fcompiler=$(F2PY_C) --opt="$(OPT)" qgsII-03.pyf \
	qgsjet-II-03.o $(RND) $(LOGO) qgsjet-II_init.o

qgsII-03.pyf: qgsjet-II-03.f qgsjet-II_init.f
	$(F2PY) -m qgsII03 -h qgsII-03.pyf --overwrite-signature \
	only: cqgsini qgsect qgini qgconf qgreg $(LOGF) : qgsjet-II-03.f \
	qgsjet-II_init.f $(LOGS)

qgsII04.$(LEXT): qgsII-04.pyf $(RND) qgsjet-II_init.o qgsjet-II-04.o 
	$(F2PY_L) -c --fcompiler=$(F2PY_C) --opt="$(OPT)" qgsII-04.pyf \
	qgsjet-II-04.o $(RND) $(LOGO) qgsjet-II_init.o

qgsII-04.pyf: qgsjet-II-04.f qgsjet-II_init.f
	$(F2PY) -m qgsII04 -h qgsII-04.pyf --overwrite-signature \
	only: cqgsini qgsect qgini qgconf qgreg $(LOGF) : qgsjet-II-04.f \
	qgsjet-II_init.f $(LOGS)

qgs01.$(LEXT): qgs01.pyf $(RND) qgsjet01_init.o qgsjet01d.o
	$(F2PY_L) -c --fcompiler=$(F2PY_C) --opt="$(OPT)" qgs01.pyf \
	qgsjet01d.o $(RND) $(LOGO) qgsjet01_init.o
	
qgs01.pyf: qgsjet01d.f qgsjet01_init.f
	$(F2PY) -m qgs01 -h qgs01.pyf --overwrite-signature \
	only: cqgsini sectnu xxaini psconf xxreg psaini $(LOGF) : qgsjet01d.f \
	qgsjet01_init.f $(LOGS)

.f.o:
	$(FC) -c $(OPT) $<
