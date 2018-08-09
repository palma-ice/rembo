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
rembo-static: $(rembo_libs) $(rembo_physics) $(rembo_base) $(rembo_tests)
	ar rc $(objdir)/librembo.a $(rembo_libs) $(rembo_physics) $(rembo_base) $(rembo_tests)
	ranlib $(objdir)/librembo.a
	@echo " "
	@echo "    $(objdir)/librembo.a is ready."
	@echo " "

rembo_test : rembo-static
		$(FC) $(DLAGS) $(FFLAGS) $(INC_COORD) $(INC_LIS) -o $(bindir)/rembo_test.x rembo_test.f90 \
			-L${CURDIR}/librembo/include -lrembo $(LFLAGS) $(objdir)/nml.o
		@echo " "
		@echo "    rembo_test.x is ready."
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
