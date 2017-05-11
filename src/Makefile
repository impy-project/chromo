#
#      Makefile for Python-libs of Fortran programs
#
#                                        (A.F. 08/12)
#
#

# Global variables for all recipes
export RND = $(CURDIR)/rangen.o
export LOGS = $(CURDIR)/logging.f
export LOGO = $(CURDIR)/logging.o

SUBDIRS = EPOS-LHC dpmjet2 pythia6 qgsjet sophia sibyll   

all: $(SUBDIRS)

.PHONY: subdirs $(SUBDIRS)
     
subdirs: $(SUBDIRS)
     
$(SUBDIRS): $(RND) $(LOGO)
	$(MAKE) -C $@

.f.o:
	$(FC) -c $(OPT) -xf77-cpp-input $<

.PHONY: clean
clean:
	rm -f $(TARGET) *.o *.prj *.chk core *.$(LEXT) *.dSYM
	for d in $(SUBDIRS); \
	do \
	$(MAKE) --directory=$$d clean; \
	done

.PHONY: distclean
distclean:
	rm -rf $(TARGET) *.o *.prj *.chk core *.$(LEXT) *.pyf *.dSYM
	for d in $(SUBDIRS); \
	do \
	$(MAKE) --directory=$$d distclean; \
	done