#
#      Makefile for Python-libs of Fortran programs
#
#                                        (A.F. 08/12)
#

URQMD_FILES = \
1fluid bessel delpart getmass hepcmp iso numrec pythia6409 \
siglookup upmerge addpart blockres coload detbal getspin \
hepnam ityp2pdg output string urqmd angdis boxprg dwidth \
init jdecay2 paulibl saveinfo tabinit whichres anndec \
cascinit dectim error hepchg input make22 proppot scatter \
uhmerge urqinit

#$(filter-out ${NOT_NEEDED:=.f}, $(wildcard ./*.f))

# URQMD_SRCS := ${URQMD_FILES:=.f}
URQMD_SRCS := $(addprefix sources/,${URQMD_FILES:=.f})
# URQMD_90SRCS := CFmax.f90 quadri.f90 cornelius.f90
URQMD_90SRCS := $(addprefix sources/, CFmax.f90 quadri.f90 cornelius.f90)
URQMD_OBJS := $(URQMD_SRCS:.f=.o) $(URQMD_90SRCS:.f90=.o)



SFUNCS=urqmd init uinit set0 params uounit strini loginit loadwtab \
norm_init output cascinit nucrad urqini partname chepevt ptsigtot init_rmmard

all: install

.PHONY: clean
clean:
	$(DEL_COMMAND) $(TARGET) *.o sources$(PATHSEP)*.o $(CLEANEXT)

.PHONY: distclean
distclean: clean
	$(DEL_COMMAND) *.pyf $(COPY_DUMP)

.PHONY: install
install: urqmd34$(LEXT)
	$(COPY_COMMAND) *$(LEXT) $(LIB_DIR)

urqmd34$(LEXT): _urqmd34.pyf $(URQMD_OBJS) impy_urqmd.o $(LOGO) $(RND)
	$(F2PY_L) -c $(F2PY_CCONF) --opt="$(OPT)" _urqmd34.pyf \
	$(LOGO) $(RND) $(URQMD_OBJS) impy_urqmd.o >> urqmd34.f2l 2>&1

_urqmd34.pyf: $(URQMD_SRCS) $(URQMD_90SRCS) $(RNDS) Makefile
	$(F2PY) -m _urqmd34 -h _urqmd34.pyf --overwrite-signature \
	only: $(SFUNCS) $(LOGF)  : $(URQMD_SRCS) $(URQMD_90SRCS) \
	$(LOGS) $(RNDS) impy_urqmd.f > urqmd34.f2l 2>&1
	$(F2PY) _urqmd34.pyf

sources/%.o: sources/%.f
	$(FC) -c $(OPT) -o $@ $<

sources/%.o: sources/%.f90
	$(FC) -c $(OPTF90) -o $@ $<

impy_urqmd.o: impy_urqmd.f
	$(FC) -c $(OPT) impy_urqmd.f