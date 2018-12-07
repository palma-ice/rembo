program rembo_test

    use ncio
    use coord

    use rembo_defs  
    use rembo 

    use insolation 

    implicit none 

    type rembo_forcing_class
        ! Climatology and forcing data for a whole year 

        real(wp), allocatable :: z_srf(:,:)      ! [m]     Surface elevation
        real(wp), allocatable :: S(:,:,:)        ! [W m-2] Insolation top-of-atmosphere
        real(wp), allocatable :: t2m(:,:,:)      ! [K]     Near-surface temperature (used for boundary)
        real(wp), allocatable :: tsl(:,:,:)      ! [K]     Sea-level temperature (used for boundary)
        real(wp), allocatable :: al_s(:,:,:)     ! [--]    Surface albedo 
        real(wp), allocatable :: co2_a           ! [ppm]   Atmospheric CO2 concentration
        real(wp), allocatable :: Z(:,:,:)        ! [m?]    Geopotential height of 750 Mb layer
    
    end type 

    type(grid_class)            :: grid 
    type(rembo_class)           :: rem1 
    type(rembo_forcing_class)   :: forc 
    
    real(wp), allocatable :: z_srf(:,:)      ! [m]     Surface elevation
    real(wp), allocatable :: f_ice(:,:)      ! [--]    Fraction of land-based ice coverage in cell
    real(wp), allocatable :: f_shlf(:,:)     ! [--]    Fraction of floating (shelf) ice coverage in cell
    real(wp), allocatable :: reg_mask(:,:)   ! [--]    Maximum region of interest
    
    character(len=56)  :: domain 
    character(len=56)  :: grid_name 
    character(len=512) :: infldr
    character(len=512) :: path_par
    character(len=512) :: path_const
    character(len=512) :: outfldr 
    character(len=512) :: file_out 
    real(wp)           :: time 
    integer            :: m, day  

    real(4) :: start, finish

    ! Start timing 
    call cpu_time(start)

    ! Determine the parameter file from the command line 
    call load_command_line_args(path_par)

    ! Assume program is running from the output folder
    outfldr = "./"

    ! Define input and output locations 
    path_const = trim(outfldr)//"rembo_const_Earth.nml"
    
    ! Initialize rembo
    call rembo_global_init(trim(path_const))
    call rembo_init(rem1,path_par=trim(path_par))
    
    domain    = trim(rem1%par%domain)
    grid_name = trim(rem1%grid%name)
    grid      = rem1%grid 

    infldr    = "ice_data/"//trim(domain)//"/"//trim(grid_name)

    ! Define output folder 
    file_out = trim(outfldr)//"rembo.nc"
    
    ! Allocate topo data and load it 
    call grid_allocate(grid,z_srf)
    call grid_allocate(grid,f_ice)
    call grid_allocate(grid,f_shlf)
    call grid_allocate(grid,reg_mask)
    call load_topo_rtopo2(z_srf,f_ice,f_shlf,reg_mask,path=trim(infldr),domain=domain,grid_name=grid%name)

    ! Load forcing
    call rembo_forc_alloc(forc,grid%G%nx,grid%G%ny)
    call load_clim_monthly_era(forc,path=trim(infldr),grid_name=trim(grid%name),z_srf=z_srf)
    
    ! Define additional forcing values 
    forc%co2_a = 350.0    ! [ppm]

    ! Define current year and update rembo (including insolation)
    time       = 0.0      ! [kyr ago]   

    call rembo_update(rem1,z_srf,f_ice,f_shlf,reg_mask,forc%t2m,forc%Z,forc%co2_a,int(time))

    ! Write final result 
    call rembo_write_init(rem1,file_out,time,units="kyr ago")
    call rembo_write_step(rem1,forc,file_out,time)

    call rembo_end(rem1)

    ! Stop timing 
    call cpu_time(finish)

    !print '("Time = ",f12.3," min.")', (finish-start)/60.0 
    print '("Time = ",f12.3," sec.")', (finish-start) 
    
contains 

    subroutine load_topo_rtopo2(z_srf,f_ice,f_shlf,reg_mask,path,domain,grid_name)
        ! Load the data into the rembo_class object

        implicit none 

        real(wp),          intent(INOUT) :: z_srf(:,:) 
        real(wp),          intent(INOUT) :: f_ice(:,:) 
        real(wp),          intent(INOUT) :: f_shlf(:,:) 
        real(wp),          intent(INOUT) :: reg_mask(:,:) 
        character(len=*),  intent(IN)    :: path 
        character(len=*),  intent(IN)    :: domain
        character(len=*),  intent(IN)    :: grid_name 

        ! Local variables
        character(len=512) :: filename
        real(wp), allocatable :: H_ice(:,:)  
        integer :: nx, ny 
        real(wp) :: region_number 

        nx = size(z_srf,1)
        ny = size(z_srf,2)

        allocate(H_ice(nx,ny))

        filename = trim(path)//"/"//trim(grid_name)//"_TOPO-RTOPO-2.0.1.nc"

        ! Static fields

        ! ## Surface elevation ##
        call nc_read(filename,"z_srf",z_srf)
        
        ! ## Ice thickness ##
        call nc_read(filename,"H_ice",H_ice)
        
        ! Ice fractions
        f_ice  = 0.0 
        where (z_srf .gt. 0.0 .and. H_ice .gt. 10.0) f_ice = 1.0 
        where (z_srf .gt. 0.0 .and. H_ice .lt. 10.0) f_ice = H_ice / 10.0 

        f_shlf = 0.0

        ! Load regions to delete regions out of interest 
        filename = trim(path)//"/"//trim(grid_name)//"_REGIONS.nc"
        call nc_read(filename,"mask",reg_mask)
        
        if (trim(domain) .eq. "Greenland") then 
            region_number = 1.3 
        else if (trim(domain) .eq. "Antarctica") then  
            region_number = 2.2 
        else 
            write(*,*) "Domain not recognized: "//trim(domain)
            stop 
        end if 
         
        where (reg_mask .ne. region_number) 
            reg_mask = 0.0 
        elsewhere 
            reg_mask = 1.0 
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
        forc%Z  = calc_geo_height(var3D)
        do m = 1, nm 
            var2Ddp = real(forc%Z(:,:,m),dp)
            call diffuse(var2Ddp,iter=2,missing_value=-9999.d0)
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

    subroutine rembo_forc_alloc(forc,nx,ny)

        implicit none 

        type(rembo_forcing_class), intent(INOUT) :: forc 
        integer,                   intent(IN)  :: nx, ny 

        if (allocated(forc%z_srf))  deallocate(forc%z_srf)
        if (allocated(forc%S))      deallocate(forc%S)
        if (allocated(forc%t2m))    deallocate(forc%t2m)
        if (allocated(forc%tsl))    deallocate(forc%tsl)
        if (allocated(forc%al_s))   deallocate(forc%al_s)
        if (allocated(forc%Z))      deallocate(forc%Z)
        
        allocate(forc%z_srf(nx,ny))

        allocate(forc%S(nx,ny,12))
        allocate(forc%t2m(nx,ny,12))
        allocate(forc%tsl(nx,ny,12))
        allocate(forc%al_s(nx,ny,12))
        allocate(forc%Z(nx,ny,12))
        
        forc%z_srf  = 0.0 
        forc%S      = 0.0 
        forc%t2m    = 0.0
        forc%tsl    = 0.0
        forc%al_s   = 0.0
        forc%Z      = 0.0

        return 

    end subroutine rembo_forc_alloc 

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
        nx = dom%grid%G%nx 
        ny = dom%grid%G%ny 
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

        do m = 1, nm

            ! Forcing fields
            call nc_write(filename,"S",dom%mon(m)%S,units="W m**-2",long_name="Insolation TOA (boundary)", &
                          dim1="xc",dim2="yc",dim3="month",start=[1,1,m],count=[nx,ny,1],ncid=ncid)
            call nc_write(filename,"t2m_bnd",dom%mon(m)%t2m_bnd,units="K",long_name="Near-surface temperature (boundary)", &
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
            call nc_write(filename,"tsurf",dom%mon(m)%tsurf,units="K",long_name="Surface temperature", &
                          dim1="xc",dim2="yc",dim3="month",start=[1,1,m],count=[nx,ny,1],ncid=ncid)
            call nc_write(filename,"uv_s",dom%mon(m)%uv_s,units="m s**-1",long_name="Near-surface wind (magnitude)", &
                          dim1="xc",dim2="yc",dim3="month",start=[1,1,m],count=[nx,ny,1],ncid=ncid)
            call nc_write(filename,"shf_s",dom%mon(m)%shf_s,units="W m**-2",long_name="Sensible heat flux", &
                          dim1="xc",dim2="yc",dim3="month",start=[1,1,m],count=[nx,ny,1],ncid=ncid)
            

            call nc_write(filename,"pr",dom%mon(m)%pr,units="mm d**-1",long_name="Precipitation", &
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

    ! # Convert geopotential into geopotential height, (m2/s2) => (m)
    elemental function calc_geo_height(phi) result(Z)
        implicit none 
        real(wp), intent(IN)  :: phi
        real(wp) :: Z
        real(wp), parameter :: g0 = 9.80665d0

        Z = Phi/g0

        return
    end function calc_geo_height

    elemental function calc_albedo_t2m(t2m,als_min,als_max,afac,tmid) result(al_s)
        ! Calculate the surface albedo of snow as a function of near-surface temperature
        ! Alexander Robinson, inspired from Slater et al, etc. 

        implicit none

        real(wp), intent(IN) :: t2m 
        real(wp), intent(IN) :: als_min, als_max, afac, tmid 
        real(wp) :: al_s 

        al_s = als_min + (als_max - als_min)*(0.5*tanh(afac*(t2m-tmid))+0.5)

        return 

    end function calc_albedo_t2m

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

end program rembo_test 

