#############################################################
##							
## Rules for individual libraries or modules
##
#############################################################

## EXTERNAL LIBRARIES #######################################

$(objdir)/basal_hydrology.o: $(libdir)/basal_hydrology.f90 $(objdir)/nml.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/geothermal.o: $(libdir)/geothermal.f90 $(objdir)/nml.o $(objdir)/ncio.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/hyster.o: $(libdir)/hyster.f90 $(objdir)/nml.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/latinhypercube.o: $(libdir)/latinhypercube.f90
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/marine_shelf.o: $(libdir)/marine_shelf.f90 $(objdir)/nml.o $(objdir)/ncio.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/isostasy.o: $(libdir)/isos/isostasy.f90 $(objdir)/nml.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/ncio.o: $(libdir)/ncio.f90
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_NC) -c -o $@ $<

$(objdir)/nml.o: $(libdir)/nml.f90
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/sealevel.o: $(libdir)/sealevel.f90 $(objdir)/nml.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/sediments.o: $(libdir)/sediments.f90 $(objdir)/nml.o $(objdir)/ncio.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/snapclim.o: $(libdir)/snapclim.f90 $(objdir)/nml.o $(objdir)/ncio.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/stommel.o: $(libdir)/stommel.f90 $(objdir)/yelmo_defs.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

# insol library 
$(objdir)/insolation.o: $(libdir)/insol/insolation.f90
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<


# smbpal library
$(objdir)/smbpal_precision.o: $(libdir)/smbpal/smbpal_precision.f90
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/smb_itm.o: $(libdir)/smbpal/smb_itm.f90 $(objdir)/smbpal_precision.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/smb_pdd.o: $(libdir)/smbpal/smb_pdd.f90 $(objdir)/smbpal_precision.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/smbpal.o: $(libdir)/smbpal/smbpal.f90 $(objdir)/smbpal_precision.o $(objdir)/nml.o $(objdir)/insolation.o  \
					$(objdir)/ncio.o \
					$(objdir)/smb_pdd.o $(objdir)/smb_itm.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

## INTERNAL PHYSICS LIBRARIES ###############################

$(objdir)/basal_dragging.o: $(srcdir)/physics/basal_dragging.f90 $(objdir)/yelmo_defs.o \
							$(objdir)/yelmo_tools.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/calving.o: $(srcdir)/physics/calving.f90 $(objdir)/yelmo_defs.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/deformation.o: $(srcdir)/physics/deformation.f90 $(objdir)/yelmo_defs.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/enthalpy.o : $(srcdir)/physics/enthalpy.f90 $(objdir)/yelmo_defs.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/icetemp.o : $(srcdir)/physics/icetemp.f90 $(objdir)/yelmo_defs.o \
					  $(objdir)/solver_tridiagonal.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/icetemp_new.o : $(srcdir)/physics/icetemp_new.f90 $(objdir)/yelmo_defs.o \
					  $(objdir)/solver_tridiagonal.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/mass_conservation.o : $(srcdir)/physics/mass_conservation.f90 $(objdir)/yelmo_defs.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/solver_ssa_grisli.o: $(srcdir)/physics/solver_ssa_grisli.f90 $(objdir)/yelmo_defs.o \
							$(objdir)/yelmo_tools.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/solver_ssa_sico32.o: $(srcdir)/physics/solver_ssa_sico32.F90 $(objdir)/yelmo_defs.o \
							$(objdir)/yelmo_tools.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_LIS) -c -o $@ $<

$(objdir)/solver_ssa_sor.o: $(srcdir)/physics/solver_ssa_sor_ij.f90 $(objdir)/yelmo_defs.o \
							$(objdir)/yelmo_tools.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/solver_tridiagonal.o: $(srcdir)/physics/solver_tridiagonal.f90 $(objdir)/yelmo_defs.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/velocity_hybrid_g11.o: $(srcdir)/physics/velocity_hybrid_g11.f90 $(objdir)/nml.o \
						  	$(objdir)/yelmo_defs.o \
						  	$(objdir)/yelmo_tools.o \
						  	$(objdir)/solver_ssa_sor.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/velocity_hybrid_pd12.o: $(srcdir)/physics/velocity_hybrid_pd12.f90 $(objdir)/nml.o \
						  	$(objdir)/yelmo_defs.o \
						  	$(objdir)/yelmo_tools.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

## YELMO BASE ###############################################

$(objdir)/yelmo_defs.o: $(srcdir)/yelmo_defs.f90 $(objdir)/nml.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/yelmo_tools.o: $(srcdir)/yelmo_tools.f90 $(objdir)/yelmo_defs.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/yelmo_topography.o: $(srcdir)/yelmo_topography.f90 $(objdir)/yelmo_defs.o \
 							  $(objdir)/yelmo_tools.o $(objdir)/mass_conservation.o $(objdir)/calving.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/yelmo_dynamics.o: $(srcdir)/yelmo_dynamics.f90 $(objdir)/yelmo_defs.o \
							$(objdir)/velocity_hybrid_pd12.o \
							$(objdir)/solver_ssa_sico32.o \
							$(objdir)/solver_ssa_grisli.o \
							$(objdir)/basal_dragging.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/yelmo_material.o: $(srcdir)/yelmo_material.f90 $(objdir)/yelmo_defs.o $(objdir)/deformation.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/yelmo_thermodynamics.o: $(srcdir)/yelmo_thermodynamics.f90 $(objdir)/yelmo_defs.o \
								  $(objdir)/icetemp.o $(objdir)/icetemp_new.o \
								  $(objdir)/enthalpy.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/yelmo_boundaries.o: $(srcdir)/yelmo_boundaries.f90 $(objdir)/yelmo_defs.o $(objdir)/ncio.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/yelmo_data.o: $(srcdir)/yelmo_data.f90 $(objdir)/yelmo_defs.o $(objdir)/nml.o $(objdir)/ncio.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/yelmo_regions.o: $(srcdir)/yelmo_regions.f90 $(objdir)/ncio.o $(objdir)/yelmo_defs.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/yelmo_io.o: $(srcdir)/yelmo_io.f90 $(objdir)/ncio.o $(objdir)/yelmo_defs.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/control.o: control.f90 $(objdir)/nml.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/yelmo_ice.o: $(srcdir)/yelmo_ice.f90 $(objdir)/yelmo_defs.o  \
					   $(objdir)/icetemp.o \
					   $(objdir)/icetemp_new.o \
				   	   $(objdir)/yelmo_topography.o \
	                   $(objdir)/yelmo_dynamics.o \
	                   $(objdir)/velocity_hybrid_pd12.o \
	                   $(objdir)/yelmo_material.o \
	                   $(objdir)/yelmo_thermodynamics.o \
	                   $(objdir)/yelmo_boundaries.o \
	                   $(objdir)/yelmo_data.o \
	                   $(objdir)/yelmo_regions.o \
	                   $(objdir)/yelmo_io.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/yelmo.o: $(srcdir)/yelmo.f90 $(objdir)/yelmo_ice.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

## YELMO TESTS ###############################################

$(objdir)/eismint.o: $(testdir)/eismint.f90 $(objdir)/ncio.o $(objdir)/yelmo_defs.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/mismip3D.o: $(testdir)/mismip3D.f90 $(objdir)/ncio.o $(objdir)/yelmo_defs.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

#############################################################
##							
## List of yelmo files
##
#############################################################

yelmo_libs = 		   $(objdir)/basal_hydrology.o \
					   $(objdir)/geothermal.o \
					   $(objdir)/hyster.o \
					   $(objdir)/insolation.o \
					   $(objdir)/isostasy.o \
					   $(objdir)/marine_shelf.o \
					   $(objdir)/nml.o \
			 		   $(objdir)/ncio.o \
			 		   $(objdir)/sealevel.o \
			 		   $(objdir)/sediments.o \
			 		   $(objdir)/smbpal_precision.o \
					   $(objdir)/smb_itm.o \
					   $(objdir)/smb_pdd.o \
					   $(objdir)/smbpal.o \
					   $(objdir)/snapclim.o

#$(objdir)/control.o 

yelmo_physics =  	   $(objdir)/basal_dragging.o \
					   $(objdir)/calving.o \
					   $(objdir)/deformation.o \
					   $(objdir)/enthalpy.o \
					   $(objdir)/icetemp.o \
					   $(objdir)/icetemp_new.o \
					   $(objdir)/mass_conservation.o \
					   $(objdir)/solver_ssa_grisli.o \
					   $(objdir)/solver_ssa_sico32.o \
					   $(objdir)/solver_ssa_sor.o \
					   $(objdir)/solver_tridiagonal.o \
					   $(objdir)/velocity_hybrid_g11.o \
					   $(objdir)/velocity_hybrid_pd12.o
				 	   
yelmo_base = 		   $(objdir)/yelmo_defs.o \
	                   $(objdir)/yelmo_tools.o \
			 		   $(objdir)/yelmo_ice.o \
			 	       $(objdir)/yelmo_topography.o \
	         		   $(objdir)/yelmo_dynamics.o \
	         		   $(objdir)/yelmo_material.o \
	         		   $(objdir)/yelmo_thermodynamics.o \
	         		   $(objdir)/yelmo_boundaries.o \
	                   $(objdir)/yelmo_data.o \
	                   $(objdir)/yelmo_regions.o \
	                   $(objdir)/yelmo_io.o \
	         		   $(objdir)/yelmo.o

yelmo_tests = 		   $(objdir)/eismint.o \
					   $(objdir)/mismip3D.o

