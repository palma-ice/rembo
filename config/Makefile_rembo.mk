#############################################################
##							
## Rules for individual libraries or modules
##
#############################################################

## EXTERNAL LIBRARIES #######################################

$(objdir)/ncio.o: $(libdir)/ncio.f90
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_NC) -c -o $@ $<

$(objdir)/nml.o: $(libdir)/nml.f90
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/solvers.o: $(libdir)/solvers.f90
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/insolation.o: $(libdir)/insol/insolation.f90
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

## REMBO BASE ###############################################

$(objdir)/rembo_defs.o: $(srcdir)/rembo_defs.f90 $(objdir)/nml.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/rembo_physics.o : $(srcdir)/rembo_physics.f90 $(objdir)/rembo_defs.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/rembo_atm.o: $(srcdir)/rembo_atm.f90 $(objdir)/rembo_defs.o $(objdir)/solvers.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/rembo_api.o: $(srcdir)/rembo_api.f90 $(objdir)/rembo_atm.o $(objdir)/insolation.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/rembo.o: $(srcdir)/rembo.f90 $(objdir)/rembo_defs.o $(objdir)/rembo_physics.o $(objdir)/rembo_api.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

## YELMO TESTS ###############################################


#############################################################
##							
## List of rembo files
##
#############################################################

rembo_libs = 		   $(objdir)/nml.o \
			 		   $(objdir)/ncio.o \
			 		   $(objdir)/solvers.o \
			 		   $(objdir)/insolation.o
			 	   
rembo_base = 		   $(objdir)/rembo_defs.o \
					   $(objdir)/rembo_physics.o \
					   $(objdir)/rembo_atm.o \
					   $(objdir)/rembo_api.o \
	         		   $(objdir)/rembo.o

rembo_tests = 		   

