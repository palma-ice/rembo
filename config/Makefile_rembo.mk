#############################################################
##							
## Rules for individual libraries or modules
##
#############################################################

## EXTERNAL LIBRARIES #######################################

$(objdir)/precision.o: $(libdir)/precision.f90
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_NC) -c -o $@ $<

$(objdir)/ncio.o: $(libdir)/ncio.f90
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_NC) -c -o $@ $<

$(objdir)/nml.o: $(libdir)/nml.f90
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/solvers.o: $(libdir)/solvers.f90 $(objdir)/precision.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/monthlydaily.o: $(libdir)/monthlydaily.f90 $(objdir)/precision.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/insolation.o: $(libdir)/insol/insolation.f90 $(objdir)/interp1D.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

# coordinates-light

$(objdir)/gaussian_filter.o: $(libdir)/coordinates-light/gaussian_filter.f90
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/grid_to_cdo.o: $(libdir)/coordinates-light/grid_to_cdo.f90
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/index.o: $(libdir)/coordinates-light/index.f90
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/interp1D.o: $(libdir)/coordinates-light/interp1D.f90
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/interp2D.o: $(libdir)/coordinates-light/interp2D.f90
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/mapping_scrip.o: $(libdir)/coordinates-light/mapping_scrip.f90 $(objdir)/ncio.o $(objdir)/interp2D.o \
								$(objdir)/gaussian_filter.o $(objdir)/index.o $(objdir)/grid_to_cdo.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

## REMBO BASE ###############################################

$(objdir)/rembo_defs.o: $(srcdir)/rembo_defs.f90 $(objdir)/nml.o $(objdir)/precision.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/rembo_grid.o : $(srcdir)/rembo_grid.f90 $(objdir)/rembo_defs.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/rembo_physics.o : $(srcdir)/rembo_physics.f90 $(objdir)/rembo_defs.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/rembo_atm.o: $(srcdir)/rembo_atm.f90 $(objdir)/rembo_defs.o $(objdir)/rembo_grid.o \
					$(objdir)/solvers.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/rembo_api.o: $(srcdir)/rembo_api.f90 $(objdir)/rembo_defs.o $(objdir)/rembo_grid.o \
					$(objdir)/rembo_atm.o  $(objdir)/insolation.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/rembo.o: $(srcdir)/rembo.f90 $(objdir)/rembo_defs.o $(objdir)/rembo_physics.o $(objdir)/rembo_api.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

## YELMO TESTS ###############################################


#############################################################
##							
## List of rembo files
##
#############################################################

rembo_libs = 		   $(objdir)/precision.o \
					   $(objdir)/nml.o \
			 		   $(objdir)/ncio.o \
			 		   $(objdir)/solvers.o \
			 		   $(objdir)/insolation.o \
					   $(objdir)/gaussian_filter.o \
					   $(objdir)/grid_to_cdo.o \
					   $(objdir)/index.o \
					   $(objdir)/interp1D.o \
					   $(objdir)/interp2D.o \
					   $(objdir)/mapping_scrip.o \
					   $(objdir)/monthlydaily.o

			 	   
rembo_base = 		   $(objdir)/rembo_defs.o \
					   $(objdir)/rembo_grid.o \
					   $(objdir)/rembo_physics.o \
					   $(objdir)/rembo_atm.o \
					   $(objdir)/rembo_api.o \
	         		   $(objdir)/rembo.o

rembo_tests = 		   

