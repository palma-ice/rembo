program rembo_test

    use coord 
    use rembo 

    implicit none 

    type(grid_class)  :: grid 
    type(rembo_class) :: rembo1 

    type rembo_forcing_class
        ! Climatology and forcing data for a whole year 

        real(wp), allocatable :: z_srf(:,:)      ! [m]     Surface elevation
        real(wp), allocatable :: f_ice(:,:)      ! [--]    Fraction of land-based ice coverage in cell
        real(wp), allocatable :: f_shlf(:,:)     ! [--]    Fraction of floating (shelf) ice coverage in cell
        
        real(wp), allocatable :: S(:,:,:)        ! [W m-2] Insolation top-of-atmosphere
        real(wp), allocatable :: t2m(:,:,:)      ! [K]     Near-surface temperature (used for boundary)
        real(wp), allocatable :: al_s(:,:,:)     ! [--]    Surface albedo 
        real(wp), allocatable :: co2_a           ! [ppm]   Atmospheric CO2 concentration
        real(wp), allocatable :: Z(:,:,:)        ! [m?]    Geopotential height of 750 Mb layer
    
    end type 

    type(rembo_forcing_class) :: forc 


    call grid_init(grid,name="GRL-20KM",mtype="stereographic",units="kilometers", &
                                   lon180=.TRUE.,dx=20.d0,nx=91,dy=20.d0,ny=151, &
                                   lambda=-40.d0,phi=72.d0,alpha=8.4d0)

    call rembo_init(rembo1,par_path="par/Greenland.nml",domain="Greenland",grid=grid)

    ! Load forcing
    call rembo_forc_alloc(forc,grid%G%nx,grid%G%ny)
    call load_clim_monthly_era(forc,path="ice_data/Greenland",grid_name=trim(grid%name))

    !call rembo_update(dom,z_srf,f_ice,f_shlf,t2m,Z,co2_a,year)

    write(*,*) "rembo_test.x finished succesfully."

contains 

    subroutine load_clim_monthly_era(forc,path,grid_name)
        ! Load the data into the era_daily_class object,
        ! Should output daily data for climatological mean
        ! *Actually loads hybrid boundary conditions: ECMWF (climate) + RCM (surface)
        ! RCM = MAR (Greenland), RCM = RACMO2 (Antarctica)

        implicit none 

        type(rembo_forcing_class), intent(INOUT) :: forc  
        character(len=*),          intent(IN)    :: path 
        character(len=*),          intent(IN)    :: grid_name 

        ! Local variables  
        real(wp), parameter :: T0 = 273.15d0 
        real(wp), allocatable :: var3D(:,:,:)
        real(dp), allocatable :: var2Ddp(:,:)
        
        character(len=512) :: filename, filename_pres 
        integer :: nx, ny, i0, i1, i2, m, nm  

        nx = size(forc%z_srf,1)
        ny = size(forc%z_srf,2)
        nm = 12 

        allocate(var3D(nx,ny,nm))
        allocate(var2Ddp(nx,ny))
        
        filename = trim(path)//"/"//trim(grid_name)//"/ERA-INT/"//trim(grid_name)//"_ERA-INT_1981-2010.nc"

        ! Static fields

        ! ## Surface elevation ##
        call nc_read(filename,"zs",forc%z_srf)
        
        ! Ice fractions
        forc%f_ice  = 0.0 
        where (forc%z_srf .gt. 0.0) forc%f_ice = 1.0 

        forc%f_shlf = 0.0 

        ! Monthly fields

        ! ## tas ## 
        call nc_read(filename,"t2m",forc%t2m)

        ! ## ttcorr ## 
        ! Get p_t[plev=650Mb]
        i0 = index(filename,"ERA-INT",back=.TRUE.)
        i1 = i0 + 7
        i2 = len_trim(filename)

        ! ## zg (750Mb) ## 
        filename_pres = filename(1:i0-1)//"ERA-INT-750Mb"//filename(i1:i2)

        call nc_read(filename_pres,"p_z",var3D)
        forc%Z  = calc_geo_height(var3D)
        do m = 1, nm 
            var2Ddp = real(forc%Z(:,:,m),dp)
            call diffuse(var2Ddp,iter=2,missing_value=-9999.d0)
            forc%Z(:,:,m) = real(var2Ddp,wp)
        end do 
        
!         ! ## dzgdx, dzgdy ##
!         ! Get horizontal gradient of geopotential height
!         do m = 1, nm 
!             call d_dx(var3Da(:,:,m),var3D(:,:,m),dx=dx)
!             call d_dy(var3Db(:,:,m),var3D(:,:,m),dx=dx)
!         end do 
!         call convert_monthly_daily_3D(var3Da,erad%dzgdx,days=erad%days)
!         call convert_monthly_daily_3D(var3Db,erad%dzgdy,days=erad%days)

!         ! ## prw ## 
!         call nc_read(filename,"tcw",var3D)
!         call convert_monthly_daily_3D(var3D,erad%prw,days=erad%days)

!         ! ## clwvi ## 
!         call nc_read(filename,"clw",var3Da)
!         call nc_read(filename,"ciw",var3Db)
!         var3D = var3Da + var3Db   ! clwvi (total = liquid + ice content)
!         call convert_monthly_daily_3D(var3D,erad%clwvi,days=erad%days)

!         ! ## hurs ## (calculated from temp, elevation)
!         ! Calculate relative humidity for each day 
!         ! Can go greater than one, where we mis-calculate tcw_sat 
!         rho_a = calc_airdens(erad%orog)
!         do k = 1, nd 
!             ! Get saturated specific humidity and total water content
!             q_sat   = calc_qsat(erad%tas(:,:,k)+erad%tascor(:,:,k),erad%orog) 
!             tcw_sat = q_sat * rho_a * H_e
!             erad%hurs(:,:,k) = erad%prw(:,:,k) / tcw_sat
!         end do 

!         ! ## u (750 Mb) ##
!         filename_pres = filename(1:i0-1)//"ERA-INT-750Mb"//filename(i1:i2)
!         call nc_read(filename_pres,"p_u",var3D)
!         call convert_monthly_daily_3D(var3D,erad%ua,days=erad%days)

!         ! ## v (750 Mb) ## 
!         call nc_read(filename_pres,"p_v",var3D)
!         call convert_monthly_daily_3D(var3D,erad%va,days=erad%days)

!         erad%ts = 273.15d0 

        return 

    end subroutine load_clim_monthly_era

    subroutine rembo_forc_alloc(forc,nx,ny)

        implicit none 

        type(rembo_forcing_class), intent(INOUT) :: forc 
        integer,                   intent(IN)  :: nx, ny 

        if (allocated(forc%z_srf))  deallocate(forc%z_srf)
        if (allocated(forc%f_ice))  deallocate(forc%f_ice)
        if (allocated(forc%f_shlf)) deallocate(forc%f_shlf)
        if (allocated(forc%S))      deallocate(forc%S)
        if (allocated(forc%t2m))    deallocate(forc%t2m)
        if (allocated(forc%al_s))   deallocate(forc%al_s)
        if (allocated(forc%Z))      deallocate(forc%Z)
        
        allocate(forc%z_srf(nx,ny))
        allocate(forc%f_ice(nx,ny))
        allocate(forc%f_shlf(nx,ny))
        
        allocate(forc%S(nx,ny,12))
        allocate(forc%t2m(nx,ny,12))
        allocate(forc%al_s(nx,ny,12))
        allocate(forc%Z(nx,ny,12))
        
        forc%z_srf  = 0.0 
        forc%f_ice  = 0.0 
        forc%f_shlf = 0.0 
        forc%S      = 0.0 
        forc%t2m    = 0.0
        forc%al_s   = 0.0
        forc%Z      = 0.0

        return 

    end subroutine rembo_forc_alloc 

    ! # Convert geopotential into geopotential height, (m2/s2) => (m)
    elemental function calc_geo_height(phi) result(Z)
        implicit none 
        real(wp), intent(IN)  :: phi
        real(wp) :: Z
        real(wp), parameter :: g0 = 9.80665d0

        Z = Phi/g0

        return
    end function calc_geo_height

end program rembo_test 

