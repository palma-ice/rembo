.SUFFIXES: .f .F .F90 .f90 .o .mod
.SHELL: /bin/sh

# PATH options
srcdir = src
objdir = librembo/include
bindir = librembo/bin
libdir = libs
testdir = tests

# Command-line options at make call
debug ?= 0
openmp ?= 1

## COMPILER CONFIGURATION ##
# (should be loaded from config directory)

<COMPILER_CONFIGURATION>

# Determine whether to use normal flags or debugging flags
DFLAGS   = $(DFLAGS_NODEBUG)
ifeq ($(debug), 1)
	DFLAGS   = $(DFLAGS_DEBUG)
endif

# Debugging flags with profiling output enabled
ifeq ($(debug), 2)
	DFLAGS   = $(DFLAGS_PROFILE)
endif

###############################################
##							
## List of rembo rules and source files
##
###############################################

include config/Makefile_rembo.mk

###############################################
##
## Compilation of complete programs
##
###############################################

# Static libraries 
rembo-static: $(rembo_libs) $(rembo_base) $(rembo_tests)
	ar rc $(objdir)/librembo.a $(rembo_libs) $(rembo_base) $(rembo_tests)
	ranlib $(objdir)/librembo.a
	@echo " "
	@echo "    $(objdir)/librembo.a is ready."
	@echo " "

test-rembo : rembo-static
		$(FC) $(DLAGS) $(FFLAGS) -o $(bindir)/test_rembo.x test_rembo.f90 \
			-L${CURDIR}/librembo/include -lrembo $(LFLAGS)
		@echo " "
		@echo "    $(bindir)/test_rembo.x is ready."
		@echo " "

test-monthlydaily : $(objdir)/monthlydaily.o $(objdir)/precision.o 
		$(FC) $(DLAGS) $(FFLAGS) -o $(bindir)/test_monthlydaily.x tests/test_monthlydaily.f90 \
				$^	$(LFLAGS)
		@echo " "
		@echo "    $(bindir)/test_monthlydaily.x is ready."
		@echo " "

test-solvers : $(objdir)/solvers.o $(objdir)/solver_linear.o $(objdir)/precision.o $(objdir)/ncio.o
		$(FC) $(DLAGS) $(FFLAGS) -o $(bindir)/test_solvers.x tests/test_solvers.f90 \
				$^	$(LFLAGS)
		@echo " "
		@echo "    $(bindir)/test_solvers.x is ready."
		@echo " "

check: $(rembo_libs) $(rembo_tests)
	@echo $(rembo_libs)
	@echo $(rembo_tests)
	@echo ""

# NOTE: eventually ncio should be extracted from all external libraries (like coord?), and then
# all external libraries should be linked by adding them to the end of the program compilation 
# command, instead of just `$(objdir)/nml.o` as it is now. 

.PHONY : usage
usage:
	@echo ""
	@echo "    * USAGE * "
	@echo ""
	@echo " make rembo 	   : compiles test_rembo.x, rembo on a given domain."
	@echo " make clean     : cleans object files"
	@echo ""

clean:
	rm -f $(bindir)/*.x
	rm -f  *.x gmon.out $(objdir)/*.o $(objdir)/*.mod $(objdir)/*.a $(objdir)/*.so
	rm -rf *.x.dSYM
