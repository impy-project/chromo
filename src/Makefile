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

define \n


endef

SUBDIRS = EPOS-LHC pythia6 qgsjet sophia sibyll dpmjet3.0-6 dpmjetIII-19.1 urqmd-3.4 

all: $(SUBDIRS)

.PHONY: subdirs $(SUBDIRS)
	 
subdirs: $(SUBDIRS)
	 
$(SUBDIRS): $(RND) $(LOGO)
	$(MAKE) -C $@ install

%.o: %.f
	$(FC) -c $(OPT) -cpp $<

%.o: %.f90
	$(FC) -c $(OPTF90) -cpp $<

.PHONY: clean
clean:
	$(DEL_COMMAND) $(TARGET) *.o *.prj *.chk core *$(LEXT) $(COPY_DUMP) *.cmdx *.cmod *.ilm *.stb *.pyd
	$(foreach d, $(SUBDIRS), $(MAKE) --directory=$(d) clean${\n})

.PHONY: distclean
distclean:
	$(DEL_COMMAND) *$(LEXT) *.pyf *.dSYM lib$(PATHSEP)* $(COPY_DUMP)
	$(foreach d, $(SUBDIRS), $(MAKE) --directory=$(d) distclean${\n})
