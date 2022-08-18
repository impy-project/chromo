#
#      Makefile for Python-libs of Fortran programs
#
#                                        (A.F. 08/12)
#

all:  install

.PHONY: clean
clean:
	$(DEL_COMMAND) $(TARGET) $(CLEANEXT)

.PHONY: distclean
distclean: clean
	rm -f *.pyf

.PHONY: install
install: qgsII03$(LEXT) qgsII04$(LEXT) qgs01$(LEXT)
	$(COPY_COMMAND) *$(LEXT) $(LIB_DIR)

#QGSJET
qgsII03$(LEXT): _qgsII03.pyf $(RND) qgsjet-II-03.o impy_qgsII.o
	$(F2PY_L) -c $(F2PY_CCONF) --opt="$(OPT)" _qgsII03.pyf \
	qgsjet-II-03.o $(RND) $(LOGO) impy_qgsII.o >> qgsII03.f2l 2>&1

_qgsII03.pyf: qgsjet-II-03.f impy_qgsII.f
	$(F2PY) -m _qgsII03 -h _qgsII03.pyf --overwrite-signature \
	only: cqgsini qgsect qgini qgconf qgreg chepevt $(LOGF) : qgsjet-II-03.f \
	impy_qgsII.f $(LOGS) > qgsII03.f2l 2>&1
	$(F2PY) _qgsII03.pyf

qgsII04$(LEXT): _qgsII04.pyf $(RND) impy_qgsII.o qgsjet-II-04.o 
	$(F2PY_L) -c $(F2PY_CCONF) --opt="$(OPT)" _qgsII04.pyf \
	qgsjet-II-04.o $(RND) $(LOGO) impy_qgsII.o >> qgsII04.f2l 2>&1

_qgsII04.pyf: qgsjet-II-04.f impy_qgsII.f
	$(F2PY) -m _qgsII04 -h _qgsII04.pyf --overwrite-signature \
	only: cqgsini qgsect qgini qgconf qgreg chepevt $(LOGF) : qgsjet-II-04.f \
	impy_qgsII.f $(LOGS) > qgsII04.f2l 2>&1
	$(F2PY) _qgsII04.pyf

qgs01$(LEXT): _qgs01.pyf $(RND) impy_qgs1.o qgsjet01d.o
	$(F2PY_L) -c $(F2PY_CCONF) --opt="$(OPT)" _qgs01.pyf \
	qgsjet01d.o $(RND) $(LOGO) impy_qgs1.o >> qgs01.f2l 2>&1
	
_qgs01.pyf: qgsjet01d.f impy_qgs1.f
	$(F2PY) -m _qgs01 -h _qgs01.pyf --overwrite-signature \
	only: cqgsini sectnu xxaini psconf xxreg psaini chepevt $(LOGF) : qgsjet01d.f \
	impy_qgs1.f $(LOGS) > qgs01.f2l 2>&1
	$(F2PY) _qgs01.pyf

.f.o:
	$(FC) -c $(OPT) $<
