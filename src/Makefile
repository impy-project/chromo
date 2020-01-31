#
#      Makefile for Python-libs of Fortran programs
#
#                                        (A.F. 08/12)
#
#

# Global variables for all recipes
export RNDS = $(CURDIR)/rangen.f
export RND = $(CURDIR)/rangen.o
export LOGS = $(CURDIR)/logging.f
export LOGO = $(CURDIR)/logging.o

SUBDIRS = EPOS-LHC pythia6 qgsjet sophia sibyll dpmjet3.0-6 urqmd-3.4 

all: $(SUBDIRS)

.PHONY: subdirs $(SUBDIRS)
     
subdirs: $(SUBDIRS)
     
$(SUBDIRS): $(RND) $(LOGO)
	$(MAKE) -C $@

%.o: %.f
	$(FC) -c $(OPT) -cpp $<

%.o: %.f90
	$(FC) -c $(OPTF90) -cpp $<

.PHONY: clean
clean:
	rm -f $(TARGET) *.o *.prj *.chk core *$(LEXT) *.dSYM
	for d in $(SUBDIRS); \
	do \
	$(MAKE) --directory=$$d clean; \
	done

.PHONY: distclean
distclean:
	rm -rf $(TARGET) *.o *.prj *.chk core *$(LEXT) *.pyf *.dSYM
	for d in $(SUBDIRS); \
	do \
	$(MAKE) --directory=$$d distclean; \
	done
