program test_rembo

    use precision, only : sp, dp, wp

    use ncio
    use rembo 
    use insolation 

    use coordinates_mapping_scrip, only : map_scrip_class, map_scrip_init, &
                                            gen_map_filename, nc_read_interp

    implicit none 

    type(rembo_class)           :: rembo1 
    type(rembo_forcing_class)   :: forc 
    
    real(wp), allocatable :: z_srf(:,:)         ! [m]     Surface elevation
    real(wp), allocatable :: f_ice(:,:)         ! [--]    Fraction of land-based ice coverage in cell
    real(wp), allocatable :: f_shlf(:,:)        ! [--]    Fraction of floating (shelf) ice coverage in cell
    logical,  allocatable :: mask_domain(:,:)   ! [--]    Maximum model region of interest
    
    character(len=56)  :: domain 
    character(len=56)  :: grid_name 
    character(len=512) :: infldr
    character(len=512) :: path_par
    character(len=512) :: outfldr 
    character(len=512) :: file_out 
    real(wp)           :: time 
    integer            :: m, day  
    integer            :: nx, ny 

    real(8) :: cpu_start_time, cpu_end_time, cpu_dtime  
    
    ! Start timing 
    call rembo_cpu_time(cpu_start_time)
    
    ! Determine the parameter file from the command line 
    call load_command_line_args(path_par)

    ! Assume program is running from the output folder
    outfldr = "./"

    ! TESTING
    call test_physics()
    stop 

    ! Initialize rembo
    call rembo_init(rembo1,path_par=trim(path_par))
    
    domain    = trim(rembo1%par%domain)
    grid_name = trim(rembo1%grid%name) 

    nx = rembo1%grid%nx 
    ny = rembo1%grid%ny 

    ! Define output folder 
    file_out = trim(outfldr)//"rembo.nc"
    
    ! Allocate topo data and load it 
    allocate(z_srf(nx,ny))
    allocate(f_ice(nx,ny))
    allocate(f_shlf(nx,ny))
    allocate(mask_domain(nx,ny))

    infldr    = "ice_data/"//trim(domain)//"/"//trim(grid_name)
    call load_topo_rtopo2(z_srf,f_ice,f_shlf,mask_domain,path=trim(infldr), &
                domain=domain,grid_name=grid_name,set_ice_free=.FALSE.)

    ! Load forcing
    call rembo_forc_alloc(forc,nx,ny)

    ! Load ERA5 data already interpolated onto Greenland grid
    !infldr    = "ice_data/"//trim(domain)//"/"//trim(grid_name)
    !call load_clim_monthly_era(forc,path=trim(infldr),grid_name=grid_name,z_srf=z_srf)

    ! Load ERA5 data from global latlon grid (best)
    infldr    = "ice_data/"
    call load_clim_monthly_era_latlon(forc,path=trim(infldr),grid_name=grid_name, &
                                                           z_srf=z_srf,dx=rembo1%grid%dx)

    ! Loading ERA40 instead...
    ! infldr    = "ice_data/"
    ! call load_clim_monthly_era40_latlon(forc,path=trim(infldr),grid_name=grid_name, &
    !                                                         z_srf=z_srf,dx=rembo1%grid%dx)
    
    ! Define additional forcing values 
    forc%co2_a = 350.0    ! [ppm]

    ! Define current year and update rembo (including insolation)
    time       = 0.0      ! [kyr ago]   

if (.TRUE.) then
    ! REMBO1
    call rembo1_update(rembo1,z_srf,f_ice,f_shlf,mask_domain,forc%t2m,forc%Z,forc%tcwv,forc%co2_a,int(time),forc%pr)
else
    ! REMBO2
    call rembo_update(rembo1,z_srf,f_ice,f_shlf,mask_domain,forc%t2m,forc%Z,forc%tcwv,forc%co2_a,int(time))
end if 

    ! Write final result 
    call rembo_write_init(rembo1,file_out,time,units="kyr ago")
    call rembo_write_step(rembo1,forc,file_out,time)

    ! Write a restart file too
    call rembo_write_state(rembo1,"rembo_restart.nc",time,units="kyr ago",init=.TRUE.)


if (.TRUE.) then
    ! Write lots of reanalysis data for offline analysis...
    infldr    = "ice_data/"
    call write_clim_monthly_era_latlon(rembo1,path=trim(infldr),grid_name=grid_name, &
                                                            z_srf=z_srf,dx=rembo1%grid%dx)
end if
    
    call rembo_end(rembo1)

    ! Stop timing 
    call rembo_cpu_time(cpu_end_time,cpu_start_time,cpu_dtime)

    write(*,"(a,f12.3,a)") "Time  = ",cpu_dtime/60.0 ," min"
    !write(*,"(a,f12.1,a)") "Speed = ",(1e-3*(time_end-time_init))/(cpu_dtime/3600.0), " kiloyears / hr"

contains 

    subroutine load_topo_rtopo2(z_srf,f_ice,f_shlf,mask_domain,path,domain,grid_name,set_ice_free)
        ! Load the data into the rembo_class object

        implicit none 

        real(wp),          intent(INOUT) :: z_srf(:,:) 
        real(wp),          intent(INOUT) :: f_ice(:,:) 
        real(wp),          intent(INOUT) :: f_shlf(:,:) 
        logical,           intent(INOUT) :: mask_domain(:,:) 
        character(len=*),  intent(IN)    :: path 
        character(len=*),  intent(IN)    :: domain
        character(len=*),  intent(IN)    :: grid_name 
        logical,           intent(IN)    :: set_ice_free

        ! Local variables
        character(len=512) :: filename
        real(wp), allocatable :: H_ice(:,:)  
        real(wp), allocatable :: z_bed(:,:)  
        real(wp), allocatable :: z_sl(:,:) 
        real(wp), allocatable :: mask_reg(:,:) 
        integer :: nx, ny 
        real(wp) :: region_number 

        nx = size(z_srf,1)
        ny = size(z_srf,2)

        allocate(H_ice(nx,ny))
        allocate(z_bed(nx,ny))
        allocate(z_sl(nx,ny))
        allocate(mask_reg(nx,ny))

        filename = trim(path)//"/"//trim(grid_name)//"_TOPO-RTOPO-2.0.1.nc"

        ! Static fields

        ! ## Surface elevation ##
        call nc_read(filename,"z_srf",z_srf)
        
        ! ## Ice thickness ##
        call nc_read(filename,"H_ice",H_ice)
        
        ! ## Bedrock elevation ##
        call nc_read(filename,"z_bed",z_bed)
        
        ! ## Sea level ##
        z_sl = 0.0 

        if (set_ice_free) then
            ! Adjust the bedrock elevation to reflect new ice thickness
            z_bed = z_bed + 910.d0/3300.d0*(H_ice)
            z_srf = z_bed
            where (z_srf .lt. z_sl) z_srf = z_sl
            H_ice = 0.0
        end if

        ! Ice fractions
        f_ice  = 0.0 
        where ( (z_srf-z_sl) .gt. 0.0 .and. H_ice .gt. 10.0) f_ice = 1.0 
        where ( (z_srf-z_sl) .gt. 0.0 .and. H_ice .lt. 10.0) f_ice = H_ice / 10.0 

        f_shlf = 0.0

        ! Load regions to delete regions out of interest 
        filename = trim(path)//"/"//trim(grid_name)//"_REGIONS.nc"
        call nc_read(filename,"mask",mask_reg)
        
        if (trim(domain) .eq. "Greenland") then 
            region_number = 1.3 
        else if (trim(domain) .eq. "Antarctica") then  
            region_number = 2.11
        else 
            write(*,*) "Domain not recognized: "//trim(domain)
            stop 
        end if 
        
        where (mask_reg .ne. region_number) 
            mask_domain = .FALSE.
        elsewhere 
            mask_domain = .TRUE.
        end where
        
        return 

    end subroutine load_topo_rtopo2

    subroutine load_clim_monthly_era(forc,path,grid_name,z_srf)
        ! Load the data into the era_daily_class object,
        ! Should output daily data for climatological mean
        ! *Actually loads hybrid boundary conditions: ECMWF (climate) + RCM (surface)
        ! RCM = MAR (Greenland), RCM = RACMO2 (Antarctica)

        implicit none 

        type(rembo_forcing_class), intent(INOUT) :: forc  
        character(len=*),          intent(IN)    :: path 
        character(len=*),          intent(IN)    :: grid_name 
        real(wp),                  intent(IN)    :: z_srf(:,:)   ! Topography to be used in model 
        ! Local variables  
        real(wp), parameter :: T0 = 273.15d0 
        real(wp), allocatable :: var3D(:,:,:)
        real(dp), allocatable :: var2Ddp(:,:)
        
        character(len=512) :: filename, filename_pres 
        integer :: nx, ny, i0, i1, i2, m, nm  
        real(wp) :: als_max, als_min, afac, tmid 

        nx = size(forc%z_srf,1)
        ny = size(forc%z_srf,2)
        nm = 12 

        allocate(var3D(nx,ny,nm))
        allocate(var2Ddp(nx,ny))
        
        filename = trim(path)//"/ERA-INT/"//trim(grid_name)//"_ERA-INT_1981-2010.nc"

        ! Static fields

        ! ## Surface elevation ##
        call nc_read(filename,"zs",forc%z_srf)
        where (forc%z_srf .lt. 0.0) forc%z_srf = 0.0 

        ! Monthly fields

        ! ## tas ## 
        call nc_read(filename,"t2m",forc%t2m)

        ! ## tsl and then correct temperature for model topography (instead of ERA topography)
        do m = 1, nm 
            forc%tsl(:,:,m) = forc%t2m(:,:,m) + 0.0065*forc%z_srf
            forc%t2m(:,:,m) = forc%tsl(:,:,m) - 0.0065*z_srf  
        end do 

        ! ## al_s ## 
        als_max =   0.80
        als_min =   0.69
        afac     =  -0.18
        tmid     = 275.35
        forc%al_s = calc_albedo_t2m(forc%t2m,als_min,als_max,afac,tmid)

        ! Define pressure-level filenames
        i0 = index(filename,"ERA-INT",back=.TRUE.)
        i1 = i0 + 7
        i2 = len_trim(filename)

        ! ## zg (750Mb) ## 
        filename_pres = filename(1:i0-1)//"ERA-INT-750Mb"//filename(i1:i2)

        call nc_read(filename_pres,"p_z",var3D)
        forc%Z  = calc_geo_height(var3D,g=real(9.80665,wp))
        do m = 1, nm 
            var2Ddp = real(forc%Z(:,:,m),dp)
            !ajr: disabled below with update to coordinates-light... to do:
            !call diffuse(var2Ddp,iter=2,missing_value=-9999.d0)
            forc%Z(:,:,m) = real(var2Ddp,wp)
        end do 
        
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

    subroutine load_clim_monthly_era_latlon(forc,path,grid_name,z_srf,dx)
        ! Load the monthly data from ERA5 on global latlon grid,
        ! interpolate online to current grid using maps.

        implicit none 

        type(rembo_forcing_class), intent(INOUT) :: forc  
        character(len=*),          intent(IN)    :: path 
        character(len=*),          intent(IN)    :: grid_name 
        real(wp),                  intent(IN)    :: z_srf(:,:)   ! Topography to be used in model 
        real(wp),                  intent(IN)    :: dx 
        
        ! Local variables  
        type(map_scrip_class) :: mps 

        character(len=512) :: filename 
        integer  :: m, nm  
        real(wp) :: als_max, als_min, afac, tmid 

        nm = 12 

        ! Intialize the map (ERA5 to grid_name)
        call map_scrip_init(mps,"ERA5",grid_name,method="con",fldr="maps",load=.TRUE.)
        

        ! ## Surface elevation ##

        filename = trim(path)//"/ERA5/era5_orography.nc"
        call nc_read_interp(filename,"z",forc%z_srf,mps=mps,method="mean")
        forc%z_srf = forc%z_srf / 9.80665_wp 
        where (forc%z_srf .lt. 0.0) forc%z_srf = 0.0 

        ! ## Near-surface air temperature (monthly) ##

        filename = trim(path)//"/ERA5/clim/"&
                        //"era5_monthly-single-levels_2m_temperature_1961-1990.nc"
        call nc_read_interp(filename,"t2m",forc%t2m,mps=mps,method="mean") !, &
                                !filt_method="gaussian",filt_par=[32e3_wp,dx])
        
        ! Get tsl and then correct temperature for model topography (instead of ERA topography)
        do m = 1, nm 
            forc%tsl(:,:,m) = forc%t2m(:,:,m) + 0.0065*forc%z_srf
            forc%t2m(:,:,m) = forc%tsl(:,:,m) - 0.0065*z_srf  
        end do 

        ! # Calculate surface albedo too
        als_max =   0.80
        als_min =   0.69
        afac     =  -0.18
        tmid     = 275.35
        forc%al_s = calc_albedo_t2m(forc%t2m,als_min,als_max,afac,tmid)
        
        ! ## Geopotential height 750Mb (monthly) ##

        filename = trim(path)//"/ERA5/clim/"&
                        //"era5_monthly-single-levels_geopotential_750_1961-1990.nc"
        call nc_read_interp(filename,"z",forc%Z,mps=mps,method="mean", &
                        filt_method="gaussian",filt_par=[32e3_wp,dx])
        forc%Z = forc%Z / 9.80665_wp

        ! ## Total column water vapour (monthly) ##

        filename = trim(path)//"/ERA5/clim/"&
                        //"era5_monthly-single-levels_total_column_water_vapour_1961-1990.nc"
        call nc_read_interp(filename,"tcwv",forc%tcwv,mps=mps,method="mean") !, &
                                !filt_method="gaussian",filt_par=[32e3_wp,dx])        
        
        ! ## Precipitation rate (monthly) ##

        filename = trim(path)//"/ERA5/clim/"&
                        //"era5_monthly-single-levels_mean_total_precipitation_rate_1961-1990.nc"
        call nc_read_interp(filename,"mtpr",forc%pr,mps=mps,method="mean") !, &
                                !filt_method="gaussian",filt_par=[32e3_wp,dx])        
        
        write(*,*) "z_srf:  ", minval(forc%z_srf),       maxval(forc%z_srf)
        write(*,*) "t2m 1:  ", minval(forc%t2m(:,:,1)),  maxval(forc%t2m(:,:,1))
        write(*,*) "t2m 7:  ", minval(forc%t2m(:,:,7)),  maxval(forc%t2m(:,:,7))
        write(*,*) "al_s 1: ", minval(forc%al_s(:,:,1)), maxval(forc%al_s(:,:,1))
        write(*,*) "al_s 7: ", minval(forc%al_s(:,:,7)), maxval(forc%al_s(:,:,7))
        write(*,*) "Z 1:    ", minval(forc%Z(:,:,1)),    maxval(forc%Z(:,:,1))
        write(*,*) "Z 7:    ", minval(forc%Z(:,:,7)),    maxval(forc%Z(:,:,7))
        write(*,*) "tcwv 1: ", minval(forc%tcwv(:,:,1)), maxval(forc%tcwv(:,:,1))
        write(*,*) "tcwv 7: ", minval(forc%tcwv(:,:,7)), maxval(forc%tcwv(:,:,7))
        write(*,*) "pr 1:   ", minval(forc%pr(:,:,1)), maxval(forc%pr(:,:,1))
        write(*,*) "pr 7:   ", minval(forc%pr(:,:,7)), maxval(forc%pr(:,:,7))
        
        write(*,*) "Loaded ERA5 boundary climate dataset."
        
        return 

    end subroutine load_clim_monthly_era_latlon

    subroutine load_clim_monthly_era40_latlon(forc,path,grid_name,z_srf,dx)
        ! Load the monthly data from ERA-40 on global latlon grid,
        ! interpolate online to current grid using maps.

        implicit none 

        type(rembo_forcing_class), intent(INOUT) :: forc  
        character(len=*),          intent(IN)    :: path 
        character(len=*),          intent(IN)    :: grid_name 
        real(wp),                  intent(IN)    :: z_srf(:,:)   ! Topography to be used in model 
        real(wp),                  intent(IN)    :: dx 
        
        ! Local variables  
        type(map_scrip_class) :: mps 

        character(len=512) :: filename 
        integer  :: m, nm  
        real(wp) :: als_max, als_min, afac, tmid 

        nm = 12 

        ! Intialize the map (ERA5 to grid_name)
        call map_scrip_init(mps,"ERA40",grid_name,method="con",fldr="maps",load=.TRUE.)
        

        ! ## Surface elevation ##

        filename = trim(path)//"/ERA40/era40-invariant.nc"
        call nc_read_interp(filename,"z",forc%z_srf,mps=mps,method="mean")
        forc%z_srf = forc%z_srf / 9.80665_wp 
        where (forc%z_srf .lt. 0.0) forc%z_srf = 0.0 

        ! ## Near-surface air temperature (monthly) ##

        filename = trim(path)//"/ERA40/clim/"&
                        //"era40-monthly-surface_1958-2001.nc"
        call nc_read_interp(filename,"t2m",forc%t2m,mps=mps,method="mean", &
                                filt_method="gaussian",filt_par=[32e3_wp,dx])
        
        ! Get tsl and then correct temperature for model topography (instead of ERA topography)
        do m = 1, nm 
            forc%tsl(:,:,m) = forc%t2m(:,:,m) + 0.0065*forc%z_srf
            forc%t2m(:,:,m) = forc%tsl(:,:,m) - 0.0065*z_srf  
        end do 

        ! # Calculate surface albedo too
        als_max =   0.80
        als_min =   0.69
        afac     =  -0.18
        tmid     = 275.35
        forc%al_s = calc_albedo_t2m(forc%t2m,als_min,als_max,afac,tmid)
        
        ! ## Geopotential height 750Mb (monthly) ##
        
        forc%Z = 0.0_wp         ! Variable not available/needed for ERA40 boundary forcing (ie, REMBO1 configuration)

        write(*,*) "z_srf:  ", minval(forc%z_srf),       maxval(forc%z_srf)
        write(*,*) "t2m 1:  ", minval(forc%t2m(:,:,1)),  maxval(forc%t2m(:,:,1))
        write(*,*) "t2m 7:  ", minval(forc%t2m(:,:,7)),  maxval(forc%t2m(:,:,7))
        write(*,*) "al_s 1: ", minval(forc%al_s(:,:,1)), maxval(forc%al_s(:,:,1))
        write(*,*) "al_s 7: ", minval(forc%al_s(:,:,7)), maxval(forc%al_s(:,:,7))
        !write(*,*) "Z 1:    ", minval(forc%Z(:,:,1)),    maxval(forc%Z(:,:,1))
        !write(*,*) "Z 7:    ", minval(forc%Z(:,:,7)),    maxval(forc%Z(:,:,7))
        
        write(*,*) "Loaded ERA-40 boundary climate dataset."
        
        return 

    end subroutine load_clim_monthly_era40_latlon

    subroutine rembo_write_step(dom,forc,filename,time)

        implicit none 
        
        type(rembo_class),          intent(IN) :: dom
        type(rembo_forcing_class),  intent(IN) :: forc 
        character(len=*),           intent(IN) :: filename
        real(wp),                   intent(IN) :: time

        ! Local variables
        integer    :: ncid, n, nx, ny, m, nm 
        real(wp) :: time_prev 

        ! Open the file for writing
        call nc_open(filename,ncid,writable=.TRUE.)

        ! Determine current writing time step (all the same year)
        n  = 1
        nx = dom%grid%nx 
        ny = dom%grid%ny 
        nm = size(dom%mon)

        ! Update the time step
        call nc_write(filename,"time",time,dim1="time",start=[n],count=[1],ncid=ncid)

        ! == Forcing fields ==
        
        ! == REMBO boundary fields ==

        call nc_write(filename,"z_srf",dom%bnd%z_srf,units="m",long_name="Surface elevation", &
                      dim1="xc",dim2="yc",ncid=ncid)
        call nc_write(filename,"f_ice",dom%bnd%f_ice,units="1",long_name="Ice fraction (grounded)", &
                      dim1="xc",dim2="yc",ncid=ncid)
        call nc_write(filename,"f_shlf",dom%bnd%f_shlf,units="1",long_name="Ice fraction (floating)", &
                      dim1="xc",dim2="yc",ncid=ncid)
        call nc_write(filename,"mask",dom%bnd%mask,units="1",long_name="Mask (solve REMBO or boundary)", &
                      dim1="xc",dim2="yc",ncid=ncid)

        call nc_write(filename,"dzsdx",dom%bnd%dzsdx,units="m/m",long_name="Surface elevation gradient (x-dir)", &
                      dim1="xc",dim2="yc",ncid=ncid)
        call nc_write(filename,"dzsdy",dom%bnd%dzsdy,units="m/m",long_name="Surface elevation gradient (y-dir)", &
                      dim1="xc",dim2="yc",ncid=ncid)
        call nc_write(filename,"dzsdxy",dom%bnd%dzsdxy,units="m/m",long_name="Surface elevation gradient (magnitude)", &
                      dim1="xc",dim2="yc",ncid=ncid)

        call nc_write(filename,"f",dom%bnd%f,units="m/m",long_name="Coriolis parameter", &
                      dim1="xc",dim2="yc",ncid=ncid)

        ! == REMBO fields == 

        call nc_write(filename,"tcm",dom%now%tcm,units="kg m**-2",long_name="Total column mass", &
                      dim1="xc",dim2="yc",ncid=ncid)
                      
        do m = 1, nm

            ! Forcing fields
            call nc_write(filename,"S",dom%mon(m)%S,units="W m**-2",long_name="Insolation TOA (boundary)", &
                          dim1="xc",dim2="yc",dim3="month",start=[1,1,m],count=[nx,ny,1],ncid=ncid)
            call nc_write(filename,"t2m_bnd",dom%mon(m)%t2m_bnd,units="K",long_name="Near-surface temperature (boundary)", &
                          dim1="xc",dim2="yc",dim3="month",start=[1,1,m],count=[nx,ny,1],ncid=ncid)
            call nc_write(filename,"tsl_bnd",dom%mon(m)%tsl_bnd,units="K",long_name="Near-surface temperature at sea level (boundary)", &
                          dim1="xc",dim2="yc",dim3="month",start=[1,1,m],count=[nx,ny,1],ncid=ncid)
            call nc_write(filename,"al_s",dom%mon(m)%al_s,units="K",long_name="Surface albedo (boundary)", &
                          dim1="xc",dim2="yc",dim3="month",start=[1,1,m],count=[nx,ny,1],ncid=ncid)
            call nc_write(filename,"Z",dom%mon(m)%Z,units="m",long_name="Geopotential height 750 Mb layer (boundary)", &
                          dim1="xc",dim2="yc",dim3="month",start=[1,1,m],count=[nx,ny,1],ncid=ncid)
            call nc_write(filename,"co2_a",dom%mon(m)%co2_a,units="ppm",long_name="Atmospheric CO2 (boundary)", &
                      dim1="xc",dim2="yc",dim3="month",start=[1,1,m],count=[nx,ny,1],ncid=ncid)
            
            ! Intermediate fields
            call nc_write(filename,"dZdx",dom%mon(m)%dZdx,units="m",long_name="Gradient Geopotential height 750 Mb layer (boundary)", &
                          dim1="xc",dim2="yc",dim3="month",start=[1,1,m],count=[nx,ny,1],ncid=ncid)
            call nc_write(filename,"dZdy",dom%mon(m)%dZdy,units="m",long_name="Gradient Geopotential height 750 Mb layer (boundary)", &
                          dim1="xc",dim2="yc",dim3="month",start=[1,1,m],count=[nx,ny,1],ncid=ncid)
            call nc_write(filename,"gamma",dom%mon(m)%gamma,units="K m**-1",long_name="Atmospheric lapse rate", &
                      dim1="xc",dim2="yc",dim3="month",start=[1,1,m],count=[nx,ny,1],ncid=ncid)
            
            ! Climate fields
            call nc_write(filename,"t2m",dom%mon(m)%t2m,units="K",long_name="Near-surface air temperature", &
                          dim1="xc",dim2="yc",dim3="month",start=[1,1,m],count=[nx,ny,1],ncid=ncid)
            call nc_write(filename,"tsl",dom%mon(m)%tsl,units="K",long_name="Near-surface air temperature at sea level", &
                          dim1="xc",dim2="yc",dim3="month",start=[1,1,m],count=[nx,ny,1],ncid=ncid)
            call nc_write(filename,"tsurf",dom%mon(m)%tsurf,units="K",long_name="Surface temperature", &
                          dim1="xc",dim2="yc",dim3="month",start=[1,1,m],count=[nx,ny,1],ncid=ncid)
            call nc_write(filename,"uv_s",dom%mon(m)%uv_s,units="m s**-1",long_name="Near-surface wind (magnitude)", &
                          dim1="xc",dim2="yc",dim3="month",start=[1,1,m],count=[nx,ny,1],ncid=ncid)
            
            call nc_write(filename,"swd",dom%mon(m)%swd,units="W m**-2",long_name="Shortwave radiation down (toa)", &
                          dim1="xc",dim2="yc",dim3="month",start=[1,1,m],count=[nx,ny,1],ncid=ncid)
            call nc_write(filename,"al_p",dom%mon(m)%al_p,units="1",long_name="Planetary albedo", &
                          dim1="xc",dim2="yc",dim3="month",start=[1,1,m],count=[nx,ny,1],ncid=ncid)
            call nc_write(filename,"at",dom%mon(m)%at,units="1",long_name="Atmospheric transmissivity", &
                          dim1="xc",dim2="yc",dim3="month",start=[1,1,m],count=[nx,ny,1],ncid=ncid)
            call nc_write(filename,"lwu",dom%mon(m)%lwu,units="W m**-2",long_name="Longwave radiation up (toa)", &
                          dim1="xc",dim2="yc",dim3="month",start=[1,1,m],count=[nx,ny,1],ncid=ncid)
            call nc_write(filename,"rco2",dom%mon(m)%rco2_a,units="W m**-2",long_name="Radiative forcing, CO2", &
                          dim1="xc",dim2="yc",dim3="month",start=[1,1,m],count=[nx,ny,1],ncid=ncid)
            

            call nc_write(filename,"swd_s",dom%mon(m)%swd_s,units="W m**-2",long_name="Shortwave radiation down (surface)", &
                          dim1="xc",dim2="yc",dim3="month",start=[1,1,m],count=[nx,ny,1],ncid=ncid)
            call nc_write(filename,"lwd_s",dom%mon(m)%lwd_s,units="W m**-2",long_name="Longwave radiation down (surface)", &
                          dim1="xc",dim2="yc",dim3="month",start=[1,1,m],count=[nx,ny,1],ncid=ncid)
            call nc_write(filename,"shf_s",dom%mon(m)%shf_s,units="W m**-2",long_name="Sensible heat flux (surface)", &
                          dim1="xc",dim2="yc",dim3="month",start=[1,1,m],count=[nx,ny,1],ncid=ncid)
            call nc_write(filename,"lhf_s",dom%mon(m)%lhf_s,units="W m**-2",long_name="Latent heat flux (surface)", &
                          dim1="xc",dim2="yc",dim3="month",start=[1,1,m],count=[nx,ny,1],ncid=ncid)
            call nc_write(filename,"lwu_s",dom%mon(m)%lwu_s,units="W m**-2",long_name="Longwave radiation up (surface)", &
                          dim1="xc",dim2="yc",dim3="month",start=[1,1,m],count=[nx,ny,1],ncid=ncid)
            
            call nc_write(filename,"ug",dom%mon(m)%ug,units="m s**-1",long_name="Geostrophic horizontal velocity, x-direction", &
                          dim1="xc",dim2="yc",dim3="month",start=[1,1,m],count=[nx,ny,1],ncid=ncid)
            call nc_write(filename,"vg",dom%mon(m)%vg,units="m s**-1",long_name="Geostrophic horizontal velocity, y-direction", &
                          dim1="xc",dim2="yc",dim3="month",start=[1,1,m],count=[nx,ny,1],ncid=ncid)
            call nc_write(filename,"uvg",dom%mon(m)%uvg,units="m s**-1",long_name="Geostrophic horizontal velocity, magnitude", &
                          dim1="xc",dim2="yc",dim3="month",start=[1,1,m],count=[nx,ny,1],ncid=ncid)
            
            call nc_write(filename,"q_sat",dom%mon(m)%q_sat,units="kg/kg",long_name="Saturated specific humidity", &
                          dim1="xc",dim2="yc",dim3="month",start=[1,1,m],count=[nx,ny,1],ncid=ncid)
            call nc_write(filename,"tcw_sat",dom%mon(m)%tcw_sat,units="kg m^-2",long_name="Saturated total column water content", &
                          dim1="xc",dim2="yc",dim3="month",start=[1,1,m],count=[nx,ny,1],ncid=ncid)
            
            call nc_write(filename,"q_s",dom%mon(m)%q_s,units="kg/kg",long_name="Specific humidity", &
                          dim1="xc",dim2="yc",dim3="month",start=[1,1,m],count=[nx,ny,1],ncid=ncid)
            call nc_write(filename,"q_r",dom%mon(m)%q_r,units="1",long_name="Relative humidity", &
                          dim1="xc",dim2="yc",dim3="month",start=[1,1,m],count=[nx,ny,1],ncid=ncid)
            call nc_write(filename,"tcw",dom%mon(m)%tcw,units="kg m^-2",long_name="Total column water content", &
                          dim1="xc",dim2="yc",dim3="month",start=[1,1,m],count=[nx,ny,1],ncid=ncid)
            call nc_write(filename,"ccw",dom%mon(m)%ccw,units="kg m^-2",long_name="Total column cloud water content", &
                          dim1="xc",dim2="yc",dim3="month",start=[1,1,m],count=[nx,ny,1],ncid=ncid)
            call nc_write(filename,"ccw_bnd",dom%mon(m)%ccw_bnd,units="kg m^-2",long_name="Total column cloud water content (boundary)", &
                          dim1="xc",dim2="yc",dim3="month",start=[1,1,m],count=[nx,ny,1],ncid=ncid)
            

            call nc_write(filename,"pr",dom%mon(m)%pr*dom%par%c%sec_day,units="mm d**-1",long_name="Precipitation", &
                          dim1="xc",dim2="yc",dim3="month",start=[1,1,m],count=[nx,ny,1],ncid=ncid)
            call nc_write(filename,"sf",dom%mon(m)%sf*dom%par%c%sec_day,units="mm d**-1",long_name="Snowfall", &
                          dim1="xc",dim2="yc",dim3="month",start=[1,1,m],count=[nx,ny,1],ncid=ncid)
            
            ! Error compared to forcing, assuming boundary field is the target
            call nc_write(filename,"t2m_err",dom%mon(m)%t2m-dom%mon(m)%t2m_bnd,units="K",long_name="Near-surface air temperature error", &
                          dim1="xc",dim2="yc",dim3="month",start=[1,1,m],count=[nx,ny,1],ncid=ncid)
            call nc_write(filename,"tsl_err",dom%mon(m)%tsl-dom%mon(m)%tsl_bnd,units="K",long_name="Near-surface air temperature at sea level error", &
                          dim1="xc",dim2="yc",dim3="month",start=[1,1,m],count=[nx,ny,1],ncid=ncid)
            call nc_write(filename,"ccw_err",dom%mon(m)%ccw-dom%mon(m)%ccw_bnd,units="kg m^-2",long_name="Total column cloud water content error", &
                          dim1="xc",dim2="yc",dim3="month",start=[1,1,m],count=[nx,ny,1],ncid=ncid)
        
        end do 
        
!         call nc_write(filename,"Ta_ann",dom%ann%t2m,units="K",long_name="Near-surface air temperature (ann)", &
!                       dim1="xc",dim2="yc",ncid=ncid)
!         call nc_write(filename,"Ta_sum",dom%ann%t2m,units="K",long_name="Near-surface air temperature (sum)", &
!                       dim1="xc",dim2="yc",ncid=ncid)
!         call nc_write(filename,"pr_ann",dom%ann%pr,units="m/a",long_name="Precipitation (ann)", &
!                       dim1="xc",dim2="yc",ncid=ncid)
!         call nc_write(filename,"sf_ann",dom%ann%sf,units="m/a",long_name="Snowfall (ann)", &
!                       dim1="xc",dim2="yc",ncid=ncid)
        
        ! Close the netcdf file
        call nc_close(ncid)

        return 

    end subroutine rembo_write_step

    subroutine load_command_line_args(path_par)

        implicit none 

        character(len=*), intent(OUT) :: path_par 

        ! Local variables 
        integer :: narg 

        narg = command_argument_count()

        if (narg .ne. 1) then 
            write(*,*) "load_command_line_args:: Error: The following &
            &argument must be provided: path_par"
            stop 
        end if 

        call get_command_argument(1,path_par)

        return 

    end subroutine load_command_line_args




    subroutine write_clim_monthly_era_latlon(dom,path,grid_name,z_srf,dx)
        ! Load the monthly data from ERA5 on global latlon grid,
        ! interpolate online to current grid using maps.

        implicit none 

        type(rembo_class),          intent(IN)  :: dom
        character(len=*),           intent(IN)  :: path 
        character(len=*),           intent(IN)  :: grid_name 
        real(wp),                   intent(IN)  :: z_srf(:,:)   ! Topography to be used in model 
        real(wp),                   intent(IN)  :: dx 
        
        ! Local variables  
        type(map_scrip_class) :: mps 
        type(rembo_forcing_class) :: forc 

        character(len=512) :: filename   
        real(wp) :: als_max, als_min, afac, tmid 

        integer    :: ncid, n, nx, ny, m, nm

        real(wp), allocatable :: z_srf_orig(:,:)
        real(wp), allocatable :: tmp2D(:,:)
        real(wp), allocatable :: tmp3D(:,:,:)

        nx = size(z_srf,1)
        ny = size(z_srf,2)
        nm = 12 

        allocate(tmp2D(nx,ny))
        allocate(tmp3D(nx,ny,nm))

        call rembo_forc_alloc(forc,nx,ny)

        ! Intialize the map (ERA5 to grid_name)
        call map_scrip_init(mps,"ERA5",grid_name,method="con",fldr="maps",load=.TRUE.)
        

        ! ## Surface elevation ##

        filename = trim(path)//"/ERA5/era5_orography.nc"
        call nc_read_interp(filename,"z",forc%z_srf,mps=mps,method="mean")
        forc%z_srf = forc%z_srf / 9.80665_wp 
        where (forc%z_srf .lt. 0.0) forc%z_srf = 0.0 

        ! ## Near-surface air temperature (monthly) ##

        filename = trim(path)//"/ERA5/clim/"&
                        //"era5_monthly-single-levels_2m_temperature_1961-1990.nc"
        call nc_read_interp(filename,"t2m",forc%t2m,mps=mps,method="mean") !, &
                                !filt_method="gaussian",filt_par=[32e3_wp,dx])
        
        ! Get tsl and then correct temperature for model topography (instead of ERA topography)
        do m = 1, nm 
            forc%tsl(:,:,m) = forc%t2m(:,:,m) + 0.0065*forc%z_srf
            forc%t2m(:,:,m) = forc%tsl(:,:,m) - 0.0065*z_srf  
        end do 

        ! # Calculate surface albedo too
        als_max =   0.80
        als_min =   0.69
        afac     =  -0.18
        tmid     = 275.35
        forc%al_s = calc_albedo_t2m(forc%t2m,als_min,als_max,afac,tmid)
        
        ! ## Geopotential height 750Mb (monthly) ##

        filename = trim(path)//"/ERA5/clim/"&
                        //"era5_monthly-single-levels_geopotential_750_1961-1990.nc"
        call nc_read_interp(filename,"z",forc%Z,mps=mps,method="mean") !, &
                        !filt_method="gaussian",filt_par=[32e3_wp,dx])
        forc%Z = forc%Z / 9.80665_wp

        ! ## Total column water vapour (monthly) ##

        filename = trim(path)//"/ERA5/clim/"&
                        //"era5_monthly-single-levels_total_column_water_vapour_1961-1990.nc"
        call nc_read_interp(filename,"tcwv",forc%tcwv,mps=mps,method="mean") !, &
                                !filt_method="gaussian",filt_par=[32e3_wp,dx])        
        
        write(*,*) "z_srf:  ", minval(forc%z_srf),       maxval(forc%z_srf)
        write(*,*) "t2m 1:  ", minval(forc%t2m(:,:,1)),  maxval(forc%t2m(:,:,1))
        write(*,*) "t2m 7:  ", minval(forc%t2m(:,:,7)),  maxval(forc%t2m(:,:,7))
        write(*,*) "al_s 1: ", minval(forc%al_s(:,:,1)), maxval(forc%al_s(:,:,1))
        write(*,*) "al_s 7: ", minval(forc%al_s(:,:,7)), maxval(forc%al_s(:,:,7))
        write(*,*) "Z 1:    ", minval(forc%Z(:,:,1)),    maxval(forc%Z(:,:,1))
        write(*,*) "Z 7:    ", minval(forc%Z(:,:,7)),    maxval(forc%Z(:,:,7))
        write(*,*) "tcwv 1: ", minval(forc%tcwv(:,:,1)), maxval(forc%tcwv(:,:,1))
        write(*,*) "tcwv 7: ", minval(forc%tcwv(:,:,7)), maxval(forc%tcwv(:,:,7))
        
        write(*,*) "Loaded ERA5 boundary climate dataset."
        

        ! === WRITE DATA TO NETCDF OUTPUT ===

        filename = trim(path)//"/ERA5/"//trim(grid_name)//"_era5.nc"

        ! Initialize netcdf file and dimensions
        call nc_create(filename)
        call nc_write_dim(filename,"xc",    x=dom%grid%xc,  units="kilometers")
        call nc_write_dim(filename,"yc",    x=dom%grid%yc,  units="kilometers")
        call nc_write_dim(filename,"month", x=1,dx=1,nx=12, units="month")

        ! Write grid information
        call rembo_grid_write(dom%grid,filename,dom%par%domain,dom%par%grid_name,create=.FALSE.)


        ! Open the file for writing actual data
        call nc_open(filename,ncid,writable=.TRUE.)

        
        ! == 2D fields ==

        call nc_write(filename,"z_srf",dom%bnd%z_srf,units="m",long_name="Surface elevation", &
                      dim1="xc",dim2="yc",ncid=ncid)
        call nc_write(filename,"f_ice",dom%bnd%f_ice,units="1",long_name="Ice fraction (grounded)", &
                      dim1="xc",dim2="yc",ncid=ncid)
        call nc_write(filename,"f_shlf",dom%bnd%f_shlf,units="1",long_name="Ice fraction (floating)", &
                      dim1="xc",dim2="yc",ncid=ncid)
        call nc_write(filename,"mask",dom%bnd%mask,units="1",long_name="Mask (solve REMBO or boundary)", &
                      dim1="xc",dim2="yc",ncid=ncid)

        call nc_write(filename,"dzsdx",dom%bnd%dzsdx,units="m/m",long_name="Surface elevation gradient (x-dir)", &
                      dim1="xc",dim2="yc",ncid=ncid)
        call nc_write(filename,"dzsdy",dom%bnd%dzsdy,units="m/m",long_name="Surface elevation gradient (y-dir)", &
                      dim1="xc",dim2="yc",ncid=ncid)
        call nc_write(filename,"dzsdxy",dom%bnd%dzsdxy,units="m/m",long_name="Surface elevation gradient (magnitude)", &
                      dim1="xc",dim2="yc",ncid=ncid)

        call nc_write(filename,"f",dom%bnd%f,units="m/m",long_name="Coriolis parameter", &
                      dim1="xc",dim2="yc",ncid=ncid)


        call nc_write(filename,"z_srf_orig",forc%z_srf,units="m",long_name="Surface elevation (input dataset)", &
                      dim1="xc",dim2="yc",ncid=ncid)
        
        ! Monthly fields
        do m = 1, nm 
            call nc_write(filename,"S",dom%mon(m)%S,units="W m**-2",long_name="Insolation TOA (boundary)", &
                        dim1="xc",dim2="yc",dim3="month",start=[1,1,m],count=[nx,ny,1],ncid=ncid)
        end do

        call nc_write(filename,"tsl",forc%tsl,units="K",long_name="Sea-level temperature", &
                        dim1="xc",dim2="yc",dim3="month",start=[1,1,1],count=[nx,ny,nm],ncid=ncid)
        call nc_write(filename,"t2m",forc%t2m,units="K",long_name="2m temperature", &
                        dim1="xc",dim2="yc",dim3="month",start=[1,1,1],count=[nx,ny,nm],ncid=ncid)
        call nc_write(filename,"al_s",forc%al_s,units="K",long_name="Surface albedo", &
                        dim1="xc",dim2="yc",dim3="month",start=[1,1,1],count=[nx,ny,nm],ncid=ncid)

        call nc_write(filename,"z750",forc%Z,units="m",long_name="Geopotential height (750 Mb)", &
                        dim1="xc",dim2="yc",dim3="month",start=[1,1,1],count=[nx,ny,nm],ncid=ncid)

        call nc_write(filename,"tcwv",forc%tcwv,units="kg m**-2",long_name="Total column water vapour", &
                        dim1="xc",dim2="yc",dim3="month",start=[1,1,1],count=[nx,ny,nm],ncid=ncid)


        ! Load and write remaining data simultaneously...

        ! ## TOA net long-wave radiation (monthly) ##

        filename = trim(path)//"/ERA5/clim/"&
                        //"era5_monthly-single-levels_mean_top_net_long_wave_radiation_flux_1961-1990.nc"
        call nc_read_interp(filename,"mtnlwrf",tmp3D,mps=mps,method="mean") !, &
                                !filt_method="gaussian",filt_par=[32e3_wp,dx])        
        call nc_write(filename,"mtnlwrf",tmp3D,units="W m**-2",long_name="Mean top net long-wave radiation flux", &
                        dim1="xc",dim2="yc",dim3="month",start=[1,1,1],count=[nx,ny,nm],ncid=ncid)

        ! ## TOA net short-wave radiation (monthly) ##

        filename = trim(path)//"/ERA5/clim/"&
                        //"era5_monthly-single-levels_mean_top_net_short_wave_radiation_flux_1961-1990.nc"
        call nc_read_interp(filename,"mtnswrf",tmp3D,mps=mps,method="mean") !, &
                                !filt_method="gaussian",filt_par=[32e3_wp,dx])        
        call nc_write(filename,"mtnswrf",tmp3D,units="W m**-2",long_name="Mean top net short-wave radiation flux", &
                        dim1="xc",dim2="yc",dim3="month",start=[1,1,1],count=[nx,ny,nm],ncid=ncid)

        ! ## TOA incident short-wave radiation (monthly) ##

        filename = trim(path)//"/ERA5/clim/"&
                        //"era5_monthly-single-levels_toa_incident_solar_radiation_1961-1990.nc"
        call nc_read_interp(filename,"tisr",tmp3D,mps=mps,method="mean") !, &
                                !filt_method="gaussian",filt_par=[32e3_wp,dx])    
        tmp3D = tmp3D / 86400.0_wp    
        call nc_write(filename,"tisr",tmp3D,units="W m**-2",long_name="Incident top net short-wave radiation flux", &
                        dim1="xc",dim2="yc",dim3="month",start=[1,1,1],count=[nx,ny,nm],ncid=ncid)
        
        ! ## Surface sensible heat flux (monthly) ##

        filename = trim(path)//"/ERA5/clim/"&
                        //"era5_monthly-single-levels_mean_surface_sensible_heat_flux_1961-1990.nc"
        call nc_read_interp(filename,"msshf",tmp3D,mps=mps,method="mean") !, &
                                !filt_method="gaussian",filt_par=[32e3_wp,dx])        
        call nc_write(filename,"msshf",tmp3D,units="W m**-2",long_name="Mean surface sensible heat flux", &
                        dim1="xc",dim2="yc",dim3="month",start=[1,1,1],count=[nx,ny,nm],ncid=ncid)

        ! ## Surface latent heat flux (monthly) ##

        filename = trim(path)//"/ERA5/clim/"&
                        //"era5_monthly-single-levels_mean_surface_latent_heat_flux_1961-1990.nc"
        call nc_read_interp(filename,"mslhf",tmp3D,mps=mps,method="mean") !, &
                                !filt_method="gaussian",filt_par=[32e3_wp,dx])        
        call nc_write(filename,"mslhf",tmp3D,units="W m**-2",long_name="Mean surface latent heat flux", &
                        dim1="xc",dim2="yc",dim3="month",start=[1,1,1],count=[nx,ny,nm],ncid=ncid)

        ! ## Surface net long-wave radiation flux (monthly) ##

        filename = trim(path)//"/ERA5/clim/"&
                        //"era5_monthly-single-levels_mean_surface_net_long_wave_radiation_flux_1961-1990.nc"
        call nc_read_interp(filename,"msnlwrf",tmp3D,mps=mps,method="mean") !, &
                                !filt_method="gaussian",filt_par=[32e3_wp,dx])        
        call nc_write(filename,"msnlwrf",tmp3D,units="W m**-2",long_name="Mean surface net long-wave radiation flux", &
                        dim1="xc",dim2="yc",dim3="month",start=[1,1,1],count=[nx,ny,nm],ncid=ncid)

        ! ## Surface net short-wave radiation flux (monthly) ##

        filename = trim(path)//"/ERA5/clim/"&
                        //"era5_monthly-single-levels_mean_surface_net_short_wave_radiation_flux_1961-1990.nc"
        call nc_read_interp(filename,"msnswrf",tmp3D,mps=mps,method="mean") !, &
                                !filt_method="gaussian",filt_par=[32e3_wp,dx])        
        call nc_write(filename,"msnswrf",tmp3D,units="W m**-2",long_name="Mean surface net short-wave radiation flux", &
                        dim1="xc",dim2="yc",dim3="month",start=[1,1,1],count=[nx,ny,nm],ncid=ncid)

        ! ## Surface pressure (monthly) ##

        filename = trim(path)//"/ERA5/clim/"&
                        //"era5_monthly-single-levels_surface_pressure_1961-1990.nc"
        call nc_read_interp(filename,"sp",tmp3D,mps=mps,method="mean") !, &
                                !filt_method="gaussian",filt_par=[32e3_wp,dx])        
        call nc_write(filename,"sp",tmp3D,units="Pa",long_name="Surface pressure", &
                        dim1="xc",dim2="yc",dim3="month",start=[1,1,1],count=[nx,ny,nm],ncid=ncid)

        ! ## Total cloud cover (monthly) ##

        filename = trim(path)//"/ERA5/clim/"&
                        //"era5_monthly-single-levels_total_cloud_cover_1961-1990.nc"
        call nc_read_interp(filename,"tcc",tmp3D,mps=mps,method="mean") !, &
                                !filt_method="gaussian",filt_par=[32e3_wp,dx])        
        call nc_write(filename,"tcc",tmp3D,units="(0 - 1)",long_name="Total cloud cover", &
                        dim1="xc",dim2="yc",dim3="month",start=[1,1,1],count=[nx,ny,nm],ncid=ncid)
        
        ! ## Total column water content (monthly) ##

        filename = trim(path)//"/ERA5/clim/"&
                        //"era5_monthly-single-levels_total_cloud_cover_1961-1990.nc"
        call nc_read_interp(filename,"tcc",tmp3D,mps=mps,method="mean") !, &
                                !filt_method="gaussian",filt_par=[32e3_wp,dx])        
        call nc_write(filename,"tcc",tmp3D,units="(0 - 1)",long_name="Total cloud cover", &
                        dim1="xc",dim2="yc",dim3="month",start=[1,1,1],count=[nx,ny,nm],ncid=ncid)

        ! ## Mean total precipitation rate (monthly) ##

        filename = trim(path)//"/ERA5/clim/"&
                        //"era5_monthly-single-levels_mean_total_precipitation_rate_1961-1990.nc"
        call nc_read_interp(filename,"mtpr",tmp3D,mps=mps,method="mean") !, &
                                !filt_method="gaussian",filt_par=[32e3_wp,dx])        
        call nc_write(filename,"mtpr",tmp3D,units="(0 - 1)",long_name="Mean total precipitation rate", &
                        dim1="xc",dim2="yc",dim3="month",start=[1,1,1],count=[nx,ny,nm],ncid=ncid)
        

        ! ## Sea surface temperature (monthly) ##

        filename = trim(path)//"/ERA5/clim/"&
                        //"era5_monthly-single-levels_sea_surface_temperature_1961-1990.nc"
        call nc_read_interp(filename,"sst",tmp3D,mps=mps,method="mean") !, &
                                !filt_method="gaussian",filt_par=[32e3_wp,dx])        
        call nc_write(filename,"sst",tmp3D,units="K",long_name="Sea surface temperature", &
                        dim1="xc",dim2="yc",dim3="month",start=[1,1,1],count=[nx,ny,nm],ncid=ncid)
        
        ! ## Sea ice cover (monthly) ##

        filename = trim(path)//"/ERA5/clim/"&
                        //"era5_monthly-single-levels_sea_ice_cover_1961-1990.nc"
        call nc_read_interp(filename,"siconc",tmp3D,mps=mps,method="mean") !, &
                                !filt_method="gaussian",filt_par=[32e3_wp,dx])        
        call nc_write(filename,"siconc",tmp3D,units="(0 - 1)",long_name="Sea ice area fraction", &
                        dim1="xc",dim2="yc",dim3="month",start=[1,1,1],count=[nx,ny,nm],ncid=ncid)
        
        ! Close the netcdf file
        call nc_close(ncid)

        write(*,*) "write_clim_monthly_era_latlon:: done."

        return 

    end subroutine write_clim_monthly_era_latlon


    subroutine test_physics()

        implicit none

        ! Local variables 
        real(wp) :: tcm

        ! Testing 
        call calc_total_column_mass(tcm,z_srf=0.0_wp,rho_0=1.3_wp,H_a=8e3_wp,H_toa=20e3_wp)
        write(*,*) "tcm = ", tcm 
        call calc_total_column_mass(tcm,z_srf=3000.0_wp,rho_0=1.3_wp,H_a=8e3_wp,H_toa=20e3_wp)
        write(*,*) "tcm = ", tcm 

        return

    end subroutine test_physics

end program test_rembo

