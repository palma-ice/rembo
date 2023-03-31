module rembo_api
    ! Wrapper to hold all modules needed for librembo 
       
    use rembo_defs 
    use rembo_grid
    use rembo_atm 
    use rembo_physics 

    use nml 
    use ncio  
    use solvers 
    use insolation 

    use coordinates_mapping_scrip, only : map_scrip_init

    implicit none 

    ! Day count for the middle of each month of the year
    integer, parameter :: mdays(12) =[30,60,90,120,150,180,210,240,270,300,330,360]-15

    private 
    public :: rembo_update
    public :: rembo_init 
    public :: rembo_end 
    public :: rembo_write_init 
    public :: rembo_grid_write 

contains 

    subroutine rembo_update(dom,z_srf,f_ice,f_shlf,reg_mask,t2m,Z,co2_a,year)
        ! Calculate atmosphere for each month of the year 

        implicit none 

        type(rembo_class), intent(INOUT) :: dom 

        real(wp), intent(IN) :: z_srf(:,:)      ! [m]     Surface elevation
        real(wp), intent(IN) :: f_ice(:,:)      ! [--]    Fraction of land-based ice coverage in cell
        real(wp), intent(IN) :: f_shlf(:,:)     ! [--]    Fraction of floating (shelf) ice coverage in cell
        real(wp), intent(IN) :: reg_mask(:,:)   ! [--]    Maximum region of interest to model 
        real(wp), intent(IN) :: t2m(:,:,:)      ! [K]     Near-surface temperature (used for boundary)
        real(wp), intent(IN) :: Z(:,:,:)        ! [m?]    Geopotential height of 750 Mb layer
        real(wp), intent(IN) :: co2_a           ! [ppm]   Atmospheric CO2 concentration
        integer,  intent(IN) :: year            ! [yrs ago (since 1950)]

        ! Local variables 
        integer  :: day, m, nm 
        
        nm = size(dom%mon)

        ! == Store annual boundary variables =====
        dom%bnd%z_srf  = z_srf 
        dom%bnd%f_ice  = f_ice 
        dom%bnd%f_shlf = f_shlf 

        ! == Calculate annual derived boundary variables =====

        ! Calculate atmospheric density from surface elevation
        dom%now%rho_a = calc_airdens(dom%bnd%z_srf)

        ! Calculate the coriolis parameter for the current grid points
        dom%bnd%f = calc_coriolis(real(dom%grid%lat,wp),dom%par%c%omega)

        ! Calculate surface gradients and total magnitude, staggered onto ac-nodes 
        call d_dx(dom%bnd%dzsdx,dom%bnd%z_srf,dx=dom%grid%dx)
        call d_dy(dom%bnd%dzsdy,dom%bnd%z_srf,dx=dom%grid%dy)
        dom%bnd%dzsdxy = calc_magnitude_from_staggered(dom%bnd%dzsdx,dom%bnd%dzsdy)

!         ! Test calculating gradient as distance to sea level 
!         call calc_gradient_to_sealevel(dom%bnd%dzsdxy,dom%bnd%z_srf,dom%bnd%z_srf*0.0, &
!                     dom%grid%x,dom%grid%y)

        ! Calculate the rembo relaxation mask
        dom%bnd%mask = gen_relaxation(dom%bnd%z_srf,dom%grid%x,dom%grid%y,radius=16.0)  
        where(reg_mask .eq. 0.0) dom%bnd%mask = 1.0 

        ! EMB OUTPUT FOR TESTING 
        call rembo_emb_write_init(dom%emb,"test.nc",dom%par%domain,dom%par%grid_name, &
                                                    time_init=real(year,wp),units="kyr ago")

        ! Loop over each month, calculate rembo atmosphere 
        do m = 1, nm 

            ! Determine representative day of the month
            day = mdays(m)

            ! == Store monthly boundary variables =====

            ! Calculate representative insolation for the month
            dom%now%S = calc_insol_day(day,dble(dom%grid%lat),dble(year),fldr="input")

            ! Save all other boundary variables 
            dom%now%t2m_bnd = t2m(:,:,m) 
            dom%now%co2_a   = co2_a 
            dom%now%Z       = Z(:,:,m)

            ! == Calculate monthly derived boundary variables =====

            ! Calc gradient of Z: dZdx, dZdy 
            call d_dx(dom%now%dZdx,dom%now%Z,dx=dom%grid%dx)
            call d_dy(dom%now%dZdy,dom%now%Z,dx=dom%grid%dy)
            ! ajr: to do: move this inside of rembo_calc_atmosphere using local variables
            
            ! == Calculate rembo atmosphere =====
            
            call rembo_calc_atmosphere(dom%now,dom%emb,dom%bnd,dom%grid,dom%par,day,year)
            
            ! Print summary 
            call rembo_print(dom,m,day,year)

            ! Store data in monthly object 
            dom%mon(m) = dom%now

            ! EMB OUTPUT FOR TESTING 
            call rembo_emb_write_step_grid(dom%emb,"test.nc",m)

        end do 

        return 

    end subroutine rembo_update
    
    subroutine rembo_init(dom,path_par,domain,grid)

        use solvers

        implicit none
        
        type(rembo_class), intent(INOUT) :: dom
        character(len=*),  intent(IN)    :: path_par  
        character(len=*),  intent(IN), optional :: domain 
        type(rgrid_class), intent(IN), optional :: grid 

        ! Local variables         
        character(len=256) :: filename, file_boundary 
        integer            :: i, m, nm 
        real(wp)           :: dt_check(2) 
        character(len=64)  :: fmt1 

        ! Load the rembo parameters
        call rembo_par_load(dom%par,trim(path_par),domain)

        ! Define rembo grid based on input filename
        call rembo_grid_define(dom%grid,dom%par%grid_name,dom%par%grid_path,grid_in=grid)

        ! Initialize grid size variables
        dom%par%npts   = dom%grid%npts 
        dom%par%nx     = dom%grid%nx 
        dom%par%ny     = dom%grid%ny 
        
        nm = size(dom%mon)

        ! Allocate boundary variables
        call rembo_bnd_alloc(dom%bnd,dom%par%nx,dom%par%ny)
        
        ! Allocate state variables
        call rembo_alloc(dom%now,dom%par%nx,dom%par%ny)
        do m = 1, nm 
            call rembo_alloc(dom%mon(m),dom%par%nx,dom%par%ny)
        end do 
        call rembo_alloc(dom%ann,dom%par%nx,dom%par%ny)

        ! Define emb diffusion grid and variables
        call rembo_grid_define(dom%emb%grid,dom%par%grid_name_emb,dom%par%grid_path_emb)

        call rembo_emb_init(dom%emb)
        
        write(*,*) "rembo_init :: allocated rembo variables."

        ! Perform mapping between emb and rembo grids
        !call map_init(dom%emb%map_toemb,  dom%grid,dom%emb%grid,max_neighbors=10,load=.TRUE.)
        !call map_init(dom%emb%map_fromemb,dom%emb%grid,dom%grid,max_neighbors=10,load=.TRUE.)
        
        ! Load the scrip map from file (should already have been generated via cdo externally)
        call map_scrip_init(dom%emb%map_toemb,dom%grid%name,dom%emb%grid%name, &
                                    method="con",fldr="maps",load=.TRUE.)
        call map_scrip_init(dom%emb%map_fromemb,dom%emb%grid%name,dom%grid%name, &
                                    method="con",fldr="maps",load=.TRUE.)
        
        ! Check diffusion time step consistency
        dt_check(1) = diff2D_timestep(dom%emb%grid%dx,dom%emb%grid%dy, &
                                      min(dom%par%en_D_sum,dom%par%en_D_win))
        dt_check(2) = diff2D_timestep(dom%emb%grid%dx,dom%emb%grid%dy, &
                                      dom%par%ccw_D)
!         dt_check(1) = diff2Dadi_timestep(dom%emb%grid%dx,dom%emb%grid%dy,dom%par%en_D)
!         dt_check(2) = diff2Dadi_timestep(dom%emb%grid%dx,dom%emb%grid%dy,dom%par%ccw_D)
        fmt1="(a,f10.1,a,f10.1,a)"
        write(*,fmt1) "Diffusion time step,   energy: ", dom%par%en_dt, " s ( max = ",dt_check(1)," s )"
        write(*,fmt1) "Diffusion time step, moisture: ", dom%par%ccw_dt," s ( max = ",dt_check(2)," s )"

        if (dom%par%en_dt .gt. dt_check(1) .or. dom%par%ccw_dt .gt. dt_check(2)) then 
            write(*,*) "Time step too big!"
            !stop 
        end if 

        write(*,*) "rembo_init:: Initialization complete for domain: "// &
                   trim(dom%par%domain)

        return 

    end subroutine rembo_init

    subroutine rembo_emb_init(emb)
        ! Initialize the diffusion variables on
        ! the desired grid resolution.
        ! grid = original grid definition of rembo
        ! emb%grid = low resolution grid of resolution dx

        implicit none

        type(diffusion_class), intent(INOUT) :: emb 

        integer :: nx, ny 

        nx = emb%grid%nx
        ny = emb%grid%ny 

        ! Topography 
        allocate(emb%mask(nx,ny))
        allocate(emb%z_srf(nx,ny))
        allocate(emb%rho_a(nx,ny))

        allocate(emb%dzsdx(nx,ny))
        allocate(emb%dzsdy(nx,ny))
        allocate(emb%dzsdxy(nx,ny))

        ! Energy balance 
        allocate(emb%tsl(nx,ny))
        allocate(emb%tsl_bnd(nx,ny))
        allocate(emb%tsl_F(nx,ny))

        ! Moisture balance
        allocate(emb%ccw(nx,ny))
        allocate(emb%ccw_bnd(nx,ny))
        allocate(emb%ccw_F(nx,ny))
        allocate(emb%ccw_cw(nx,ny))
        allocate(emb%ccw_pr(nx,ny))
        allocate(emb%tcw(nx,ny))
        allocate(emb%tcw_bnd(nx,ny))

        ! Energy and moisture balance
        allocate(emb%ug(nx,ny))
        allocate(emb%vg(nx,ny))
        allocate(emb%uvg(nx,ny))
        allocate(emb%ww(nx,ny))

        allocate(emb%q_r(nx,ny))

        ! Diffusion
        allocate(emb%kappa(nx,ny))
        allocate(emb%kappaw(nx,ny))


        ! Initialize all variables to zero 

        emb%mask        = 0.0
        emb%z_srf       = 0.0
        emb%rho_a       = 0.0

        emb%dzsdx       = 0.0
        emb%dzsdy       = 0.0
        emb%dzsdxy      = 0.0

        ! Energy balance 
        emb%tsl         = 0.0
        emb%tsl_bnd     = 0.0
        emb%tsl_F       = 0.0

        ! Moisture balance
        emb%ccw         = 0.0
        emb%ccw_bnd     = 0.0
        emb%ccw_F       = 0.0
        emb%ccw_cw      = 0.0
        emb%ccw_pr      = 0.0
        emb%tcw         = 0.0
        emb%tcw_bnd     = 0.0

        ! Energy and moisture balance
        emb%ug          = 0.0
        emb%vg          = 0.0
        emb%uvg         = 0.0
        emb%ww          = 0.0

        emb%q_r         = 0.0

        ! Diffusion
        emb%kappa       = 0.0
        emb%kappaw      = 0.0 

        return 

    end subroutine rembo_emb_init

    subroutine rembo_end(dom)

        implicit none 

        type(rembo_class), intent(INOUT) :: dom

        ! Local variables 
        integer :: m, nm 

        nm = size(dom%mon)
        
        call rembo_bnd_dealloc(dom%bnd)

        call rembo_dealloc(dom%now)
        do m = 1, nm 
            call rembo_dealloc(dom%mon(m))
        end do
        call rembo_dealloc(dom%ann)

        write(*,*) "rembo_end :: deallocated rembo variables: "//trim(dom%par%domain)
        write(*,*)

        return 

    end subroutine rembo_end

    subroutine rembo_grid_define(grid,grid_name,grid_path,grid_in)

        implicit none 

        type(rgrid_class), intent(OUT)  :: grid  
        character(len=*), intent(INOUT) :: grid_name      ! Overwritable if grid_in is present
        character(len=*), intent(IN)    :: grid_path
        type(rgrid_class), intent(IN), optional :: grid_in 

        if (present(grid_in)) then 

            grid = grid_in 

            ! Ensure parameter grid_name is consistent with defined grid 
            grid_name = grid%name 
        
        else 
            ! Define rembo grid from predefined options (currently only from file)

            write(*,*) "grid_path: ", trim(grid_path) 

            ! Use grid_path to load grid definition from NetCDF file 
            call rembo_init_grid(grid,grid_path,grid_name)

        end if 

        return 

    end subroutine rembo_grid_define
    
    subroutine rembo_print(dom,m,d,year)

        implicit none 

        type(rembo_class) :: dom 
        integer  :: m, d, year
        real(wp) :: npts 
        character(len=256) :: fmt_head, fmt_table 

        npts = real(count(dom%bnd%mask>0),wp)

        fmt_head  = "(a6,a5,  5a8,2x,  a8,  3a8,2x,  4a8)"
        fmt_table = "(i6,i5,5f8.2,2x,f8.2,3f8.2,2x,4f8.2)"

        if (m .eq. 1) &
        write(*,fmt_head) "year","day","swd","tas", &
                           "tcw","ccw","c_w","pr"
        
        if (npts>0) then 

            write(*,fmt_table) year, d, &
                sum(dom%now%swd,        mask=dom%bnd%mask>0)/npts, &     ! [W/m^2]
                sum(dom%now%t2m,        mask=dom%bnd%mask>0)/npts, &     ! [K]
                sum(dom%now%tcw,        mask=dom%bnd%mask>0)/npts, &     ! [mm]
                sum(dom%now%ccw,        mask=dom%bnd%mask>0)/npts, &     ! [mm]
                sum(dom%now%c_w*dom%par%c%sec_day,mask=dom%bnd%mask>0)/npts, &     ! [mm/d]
                sum(dom%now%pr*dom%par%c%sec_day, mask=dom%bnd%mask>0)/npts        ! [mm/d]

        else 
            ! Print empty table 
            write(*,fmt_table) year, d, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 

        end if 

        return 

    end subroutine rembo_print
    
    subroutine rembo_par_load(par,filename,domain)

        type(rembo_param_class),    intent(INOUT) :: par 
        character(len=*),           intent(IN)    :: filename 
        character(len=*),           intent(IN), optional :: domain 

        ! Local variables 
        character(len=56) :: group0
        character(len=56) :: group1 
        character(len=56) :: group2 
        
        par%domain = trim(domain)

        group0 = "rembo_phys_const"
        group1 = "rembo" 
        group2 = "rembo1" 

        ! First load physical constants for this domain 
        call nml_read(filename,group0,"nm",             par%c%nm)
        call nml_read(filename,group0,"ndm",            par%c%ndm)
        call nml_read(filename,group0,"sec_year",       par%c%sec_year)
        call nml_read(filename,group0,"sec_day_ref",    par%c%sec_day_ref)
        call nml_read(filename,group0,"g",              par%c%g)
        call nml_read(filename,group0,"omega",          par%c%omega)
        call nml_read(filename,group0,"T0",             par%c%T0)
        call nml_read(filename,group0,"rho_ice",        par%c%rho_ice)
        call nml_read(filename,group0,"rho_w",          par%c%rho_w)
        
        call nml_read(filename,group1,"domain",         par%domain)
        call nml_read(filename,group1,"grid_name",      par%grid_name)
        call nml_read(filename,group1,"grid_name_emb",  par%grid_name_emb)
        call nml_read(filename,group1,"grid_path",      par%grid_path)
        call nml_read(filename,group1,"grid_path_emb",  par%grid_path_emb)
        call nml_read(filename,group1,"restart",        par%restart)

        call nml_read(filename,group1,"H_e",            par%H_e)
        call nml_read(filename,group1,"rembo1",         par%rembo1)
        call nml_read(filename,group1,"en_dt",          par%en_dt)
        call nml_read(filename,group1,"en_D",           par%en_D)
        call nml_read(filename,group1,"en_D_win",       par%en_D_win)
        call nml_read(filename,group1,"en_D_sum",       par%en_D_sum)
        call nml_read(filename,group1,"en_kr",          par%en_kr)
        call nml_read(filename,group1,"en_kz",          par%en_kz)
        call nml_read(filename,group1,"en_kl",          par%en_kl)
        call nml_read(filename,group1,"en_kdT",         par%en_kdT)
        call nml_read(filename,group1,"ccw_dt",         par%ccw_dt)
        call nml_read(filename,group1,"ccw_D",          par%ccw_D)
        call nml_read(filename,group1,"ccw_kr",         par%ccw_kr)
        call nml_read(filename,group1,"k_c",            par%k_c)
        call nml_read(filename,group1,"k_p",            par%k_p)
        call nml_read(filename,group1,"k_z",            par%k_z)
        call nml_read(filename,group1,"k_x",            par%k_x)
        call nml_read(filename,group1,"e0",             par%e0)
        call nml_read(filename,group1,"c1",             par%c1)
        call nml_read(filename,group1,"k_p_now",        par%k_p_now)
        call nml_read(filename,group1,"k_t",            par%k_t)
        call nml_read(filename,group1,"k_w",            par%k_w)
        call nml_read(filename,group1,"sf_a",           par%sf_a)
        call nml_read(filename,group1,"sf_b",           par%sf_b)
        call nml_read(filename,group1,"f_k",            par%f_k)
        call nml_read(filename,group1,"nk1",            par%nk1)
        call nml_read(filename,group1,"nk2",            par%nk2)
        call nml_read(filename,group1,"nk3",            par%nk3)
        call nml_read(filename,group1,"gamma",          par%gamma)
        call nml_read(filename,group1,"gamma2",         par%gamma2)
        call nml_read(filename,group1,"S0",             par%S0)
        call nml_read(filename,group1,"Lw",             par%Lw)
        call nml_read(filename,group1,"Lm",             par%Lm)
        call nml_read(filename,group1,"Ls",             par%Ls)
        call nml_read(filename,group1,"rho_w",          par%rho_w)
        call nml_read(filename,group1,"cp",             par%cp)
        call nml_read(filename,group1,"ci",             par%ci)
        call nml_read(filename,group1,"alp_a",          par%alp_a)
        call nml_read(filename,group1,"alp_b",          par%alp_b)
        
        call nml_read(filename,group1,"alp_c",          par%alp_c)
        call nml_read(filename,group1,"lwu_a",          par%lwu_a)
        call nml_read(filename,group1,"lwu_b",          par%lwu_b)
        call nml_read(filename,group1,"lwu_c",          par%lwu_c)
        call nml_read(filename,group1,"swds_a",         par%swds_a)
        call nml_read(filename,group1,"swds_b",         par%swds_b)
        call nml_read(filename,group1,"swds_c",         par%swds_c)
        call nml_read(filename,group1,"lwds_a",         par%lwds_a)
        call nml_read(filename,group1,"lwds_b",         par%lwds_b)
        call nml_read(filename,group1,"lwds_c",         par%lwds_c)
        call nml_read(filename,group1,"lwds_d",         par%lwds_d)
        
        call nml_read(filename,group1,"shfs_Cm",        par%shfs_Cm)
        call nml_read(filename,group1,"shfs_p",         par%shfs_p)

        ! === REMBO1 parameters ===
        call nml_read(filename,group2,"ta",             par%r1_ta)
        call nml_read(filename,group2,"tb",             par%r1_tb)
        ! =========================

        ! Get additional derived constants   
        par%c%nd         = par%c%nm * par%c%ndm

        par%c%day_year   = real(par%c%nd,wp)
        par%c%day_month  = real(par%c%ndm,wp)
        par%c%month_year = real(par%c%nm,wp)

        par%c%sec_day    = par%c%sec_year / par%c%day_year   ! 8.765813e4
        par%c%sec_frac   = par%c%sec_day  / par%c%sec_day_ref
        
        ! How many time steps in 1 day, aprx?
        par%en_nstep  = floor(par%c%sec_day / par%en_dt)
        par%ccw_nstep = floor(par%c%sec_day / par%ccw_dt)
        
        ! Overwrite parameter values with argument definitions if available
        if (present(domain))     par%domain    = trim(domain)
        
        ! Make sure to parse grid_path and grid_path_emb
        call rembo_parse_path(par%grid_path,par%domain,par%grid_name)
        call rembo_parse_path(par%grid_path_emb,par%domain,par%grid_name_emb)
        
        if (rembo_write_log) then 
            write(*,*) "yelmo:: loaded global constants:"
            write(*,*) "    nm       = ", par%c%nm 
            write(*,*) "    ndm      = ", par%c%ndm 
            write(*,*) "    sec_year = ", par%c%sec_year 
            write(*,*) "    sec_day  = ", par%c%sec_day 
            write(*,*) "    g        = ", par%c%g 
            write(*,*) "    omega    = ", par%c%omega
            write(*,*) "    T0       = ", par%c%T0 
            write(*,*) "    rho_ice  = ", par%c%rho_ice 
            write(*,*) "    rho_w    = ", par%c%rho_w 
        end if 

        return

    end subroutine rembo_par_load


    subroutine rembo_alloc(now,nx,ny)

        implicit none 

        type(rembo_state_class), intent(INOUT) :: now 
        integer, intent(IN) :: nx, ny 

        ! Make sure this object is fully deallocated first
        call rembo_dealloc(now)

        allocate(now%S(nx,ny))       ! Insolation top of the atmosphere (W/m2)
        allocate(now%t2m_bnd(nx,ny)) ! Near-surface temp
        allocate(now%al_s(nx,ny))    ! Surface albedo (0 - 1)
        allocate(now%co2_a(nx,ny))   ! Atmospheric CO2 (ppm)
        allocate(now%Z(nx,ny))       ! Geopotential height at input pressure level (eg 750Mb) (m)
        allocate(now%dZdx(nx,ny))    ! Geopotential height gradient (m m**-1)
        allocate(now%dZdy(nx,ny))    ! Geopotential height gradient (m m**-1)

        allocate(now%rco2_a(nx,ny))  ! Radiative forcing of CO2 (W m-2)
        allocate(now%rho_a(nx,ny))   ! Air density (kg m-3)
        allocate(now%sp(nx,ny))      ! Surface pressure (Pa)
        
        allocate(now%gamma(nx,ny))   ! Temperature lapse rate
        allocate(now%t2m(nx,ny))     ! Near-surface temp
        allocate(now%ct2m(nx,ny))    ! Near-surface temp inversion correction
        allocate(now%tsurf(nx,ny))   ! Near-surface temp
        allocate(now%pr(nx,ny))      ! Precipitation
        allocate(now%sf(nx,ny))      ! Precipitation (snow)
        allocate(now%q_s(nx,ny))     ! Specific humidity at the surface (kg/kg)
        allocate(now%q_sat(nx,ny))   ! Saturated specific humidity at the surface (kg/kg)
        allocate(now%q_r(nx,ny))     ! Relative humidity (0 - 1)
        allocate(now%tcw(nx,ny))     ! Total water content (kg/m2)
        allocate(now%tcw_sat(nx,ny)) ! Saturated total water content (kg/m2)
        allocate(now%ccw_prev(nx,ny))! Previous tcw value
        allocate(now%ccw(nx,ny))     ! Cloud water content (kg/m2)
        allocate(now%c_w(nx,ny))     ! Condensated water (kg/m2)
        allocate(now%ug(nx,ny))      ! Horizontal x-component 750Mb velocity (m/s)
        allocate(now%vg(nx,ny))      ! Horizontal y-component 750Mb velocity (m/s)
        allocate(now%uvg(nx,ny))     ! Horizontal magnitude 750Mb velocity (m/s)
        allocate(now%ww(nx,ny))      ! Vertical velocity 750Mb (m/s)
        allocate(now%cc(nx,ny))      ! Cloud fraction (0 - 1)        
        allocate(now%swd(nx,ny))     ! Short-wave downward at toa (W/m2)
        allocate(now%lwu(nx,ny))     ! Long-wave upward radiation at toa (W/m2)
        allocate(now%al_p(nx,ny))    ! Planetary albedo (0 - 1)
        allocate(now%at(nx,ny))      ! Atmospheric transmissivity (0 - 1)
        
        allocate(now%u_k(nx,ny))     ! x-component katabatic wind velocity (m/s)
        allocate(now%v_k(nx,ny))     ! y-component katabatic wind velocity (m/s)
        allocate(now%uv_k(nx,ny))    ! magnitude katabatic wind velocity (m/s)
        
        allocate(now%dtsldx(nx,ny))     ! temperature gradient [K m-1]
        allocate(now%dtsldy(nx,ny))     ! temperature gradient [K m-1]
        allocate(now%dtsldxy(nx,ny))    ! temperature gradient [K m-1]
        
        allocate(now%swd_s(nx,ny))   ! Short-wave downward at surface (W/m2)
        allocate(now%lwd_s(nx,ny))   ! Long-wave downward at surface (W/m2)
        allocate(now%lwu_s(nx,ny))   ! Long-wave upward at surface (W/m2)
        allocate(now%shf_s(nx,ny))   ! Sensible heat flux at surface (W/m2)
        allocate(now%lhf_s(nx,ny))   ! Latent heat flux at surface (W/m2)
        allocate(now%u_s(nx,ny))     ! Horizontal x-component surface velocity (m/s)
        allocate(now%v_s(nx,ny))     ! Horizontal y-component surface velocity (m/s)
        allocate(now%uv_s(nx,ny))    ! Horizontal magnitude surface velocity (m/s)

        now%S           = 0.0 
        now%t2m_bnd     = 0.0 
        now%al_s        = 0.0 
        now%co2_a       = 0.0 
        now%Z           = 0.0 
        now%dZdx        = 0.0   
        now%dZdy        = 0.0  

        now%rco2_a      = 0.0 
        now%rho_a       = 0.0 
        now%sp          = 0.0 

        now%gamma       = 0.0
        now%t2m         = 0.0
        now%ct2m        = 0.0
        now%tsurf       = 0.0
        now%pr          = 0.0        
        now%sf          = 0.0
        now%q_s         = 0.0 
        now%q_sat       = 0.0 
        now%q_r         = 0.0 
        now%tcw         = 0.0 
        now%tcw_sat     = 0.0 
        now%ccw_prev    = 0.0
        now%ccw         = 0.0 
        now%c_w         = 0.0  
        now%ug          = 0.0 
        now%vg          = 0.0   
        now%uvg         = 0.0
        now%ww          = 0.0 
        now%cc          = 0.0           
        now%swd         = 0.0    
        now%lwu         = 0.0 
        now%al_p        = 0.0  
        now%at          = 0.0  
        
        now%u_k         = 0.0    
        now%v_k         = 0.0 
        now%uv_k        = 0.0   
        
        now%dtsldx      = 0.0  
        now%dtsldy      = 0.0 
        now%dtsldxy     = 0.0  
          
        now%swd_s       = 0.0  
        now%lwd_s       = 0.0  
        now%lwu_s       = 0.0  
        now%shf_s       = 0.0  
        now%lhf_s       = 0.0 
        now%u_s         = 0.0  
        now%v_s         = 0.0 
        now%uv_s        = 0.0 

        return

    end subroutine rembo_alloc

    subroutine rembo_dealloc(now)

        implicit none 

        type(rembo_state_class), intent(INOUT) :: now

        ! Deallocate state variables 

        if (allocated(now%S   ))        deallocate(now%S)       ! Insolation top of the atmosphere (W/m2)
        if (allocated(now%t2m_bnd ))    deallocate(now%t2m_bnd) ! Near-surface temp
        if (allocated(now%al_s ))       deallocate(now%al_s)    ! Surface albedo (0 - 1)
        if (allocated(now%co2_a ))      deallocate(now%co2_a)   ! Atmospheric CO2 (ppm)
        if (allocated(now%Z ))          deallocate(now%Z)       ! Geopotential height at input pressure level (eg 750Mb) (m)
        if (allocated(now%dZdx ))       deallocate(now%dZdx)    ! Geopotential height gradient (m m**-1)
        if (allocated(now%dZdy ))       deallocate(now%dZdy)    ! Geopotential height gradient (m m**-1)
        
        if (allocated(now%rco2_a) )     deallocate(now%rco2_a)  ! Radiative forcing of CO2 (W m-2)
        if (allocated(now%rho_a)  )     deallocate(now%rho_a)   ! Air density (kg m-3)
        if (allocated(now%sp)     )     deallocate(now%sp)      ! Surface pressure (Pa)

        if (allocated(now%gamma) )      deallocate(now%gamma)   ! Temperature lapse rate 
        if (allocated(now%t2m) )        deallocate(now%t2m)     ! Near-surface temp
        if (allocated(now%ct2m) )       deallocate(now%ct2m)    ! Near-surface temp inversion correction
        if (allocated(now%tsurf) )      deallocate(now%tsurf)   ! Surface temp
        if (allocated(now%pr) )         deallocate(now%pr)      ! Precipitation
        if (allocated(now%sf) )         deallocate(now%sf)      ! Precipitation (snow)
        if (allocated(now%q_s) )        deallocate(now%q_s)     ! Specific humidity at the surface (kg/kg)
        if (allocated(now%q_sat) )      deallocate(now%q_sat)   ! Saturated specific humidity at the surface (kg/kg)
        if (allocated(now%q_r) )        deallocate(now%q_r)     ! Relative humidity (0 - 1)
        if (allocated(now%tcw) )        deallocate(now%tcw)     ! Total water content (kg/m2)
        if (allocated(now%tcw_sat) )    deallocate(now%tcw_sat) ! Saturated total water content (kg/m2)
        if (allocated(now%ccw_prev) )   deallocate(now%ccw_prev)! Previous tcw value
        if (allocated(now%ccw) )        deallocate(now%ccw)     ! Cloud water content (kg/m2)
        if (allocated(now%c_w) )        deallocate(now%c_w)     ! Condensated water (kg/m2)
        if (allocated(now%ug) )         deallocate(now%ug)      ! Horizontal x-component 750Mb velocity (m/s)
        if (allocated(now%vg) )         deallocate(now%vg)      ! Horizontal y-component 750Mb velocity (m/s)
        if (allocated(now%uvg) )        deallocate(now%uvg)     ! Horizontal magnitude 750Mb velocity (m/s)
        if (allocated(now%ww) )         deallocate(now%ww)      ! Vertical velocity 750Mb (m/s)
        if (allocated(now%cc) )         deallocate(now%cc)      ! Cloud fraction (0 - 1)        
        if (allocated(now%swd) )        deallocate(now%swd)     ! Short-wave downward at toa (W/m2)
        if (allocated(now%lwu) )        deallocate(now%lwu)     ! Long-wave upward radiation at toa (W/m2)
        if (allocated(now%al_p) )       deallocate(now%al_p)    ! Planetary albedo (0 - 1)
        if (allocated(now%at) )         deallocate(now%at)      ! Atmospheric transmissivity (0 - 1)

        if (allocated(now%u_k) )        deallocate(now%u_k)     ! x-component katabatic wind velocity (m/s)
        if (allocated(now%v_k) )        deallocate(now%v_k)     ! y-component katabatic wind velocity (m/s)
        if (allocated(now%uv_k) )       deallocate(now%uv_k)    ! magnitude katabatic wind velocity (m/s)
        
        if (allocated(now%dtsldx) )     deallocate(now%dtsldx)     ! temperature gradient [K m-1]
        if (allocated(now%dtsldy) )     deallocate(now%dtsldy)     ! temperature gradient [K m-1]
        if (allocated(now%dtsldxy) )    deallocate(now%dtsldxy)    ! temperature gradient [K m-1]
        
        if (allocated(now%swd_s) )      deallocate(now%swd_s)   ! Short-wave downward at surface (W/m2)
        if (allocated(now%lwd_s) )      deallocate(now%lwd_s)   ! Long-wave downward at surface (W/m2)
        if (allocated(now%lwu_s) )      deallocate(now%lwu_s)   ! Long-wave upward at surface (W/m2)
        if (allocated(now%shf_s) )      deallocate(now%shf_s)   ! Sensible heat flux at surface (W/m2)
        if (allocated(now%lhf_s) )      deallocate(now%lhf_s)   ! Latent heat flux at surface (W/m2)
        if (allocated(now%u_s) )        deallocate(now%u_s)     ! Horizontal x-component surface velocity (m/s)
        if (allocated(now%v_s) )        deallocate(now%v_s)     ! Horizontal y-component surface velocity (m/s)
        if (allocated(now%uv_s) )       deallocate(now%uv_s)    ! Horizontal magnitude surface velocity (m/s)

        return 

    end subroutine rembo_dealloc

    subroutine rembo_bnd_alloc(bnd,nx,ny)

        implicit none 
         
        type(rembo_boundary_class), intent(INOUT) :: bnd
        integer, intent(IN) :: nx, ny 

        call rembo_bnd_dealloc(bnd)

        allocate(bnd%z_srf(nx,ny))   ! Surface elevation 
        allocate(bnd%f_ice(nx,ny))   ! Ice thickness (grounded)
        allocate(bnd%f_shlf(nx,ny))  ! Ice thickness (floating)
        
        allocate(bnd%mask(nx,ny))    ! Ocean-land-ice mask 
        allocate(bnd%f(nx,ny))       ! Coriolis parameter (1/s)
        allocate(bnd%dzsdx(nx,ny))   ! Surface gradient (x-dir)
        allocate(bnd%dzsdy(nx,ny))   ! Surface gradient (y-dir)
        allocate(bnd%dzsdxy(nx,ny))  ! Surface gradient (magnitude)
        
        bnd%z_srf       = 0.0
        bnd%f_ice       = 0.0 
        bnd%f_shlf      = 0.0 

        bnd%mask        = 0.0
        bnd%f           = 0.0 
        bnd%dzsdx       = 0.0 
        bnd%dzsdy       = 0.0 
        bnd%dzsdxy      = 0.0    
        
        return 

    end subroutine rembo_bnd_alloc

    subroutine rembo_bnd_dealloc(bnd)

        implicit none 
         
        type(rembo_boundary_class), intent(INOUT) :: bnd

        if (allocated(bnd%z_srf ))  deallocate(bnd%z_srf)   ! Surface elevation 
        if (allocated(bnd%f_ice ))  deallocate(bnd%f_ice)   ! Ice thickness (grounded)
        if (allocated(bnd%f_shlf )) deallocate(bnd%f_shlf)  ! Ice thickness (floating)
        
        if (allocated(bnd%mask ))   deallocate(bnd%mask)    ! Ocean-land-ice mask 
        if (allocated(bnd%f ))      deallocate(bnd%f)       ! Coriolis parameter (1/s)
        if (allocated(bnd%dzsdx ))  deallocate(bnd%dzsdx)   ! Surface gradient (x-dir)
        if (allocated(bnd%dzsdy ))  deallocate(bnd%dzsdy)   ! Surface gradient (y-dir)
        if (allocated(bnd%dzsdxy )) deallocate(bnd%dzsdxy)  ! Surface gradient (magnitude)

        return 

    end subroutine rembo_bnd_dealloc

    subroutine rembo_write_init(dom,filename,time_init,units)

        implicit none 

        type(rembo_class), intent(IN) :: dom 
        character(len=*),  intent(IN) :: filename 
        real(wp),          intent(IN) :: time_init
        character(len=*),  intent(IN) :: units

        ! Local variables 

        ! Initialize netcdf file and dimensions
        call nc_create(filename)
        call nc_write_dim(filename,"xc",    x=dom%grid%xc,  units="kilometers")
        call nc_write_dim(filename,"yc",    x=dom%grid%yc,  units="kilometers")
        call nc_write_dim(filename,"month", x=1,dx=1,nx=12, units="month")
        call nc_write_dim(filename,"time",  x=time_init,dx=1.0_wp,nx=1,units=trim(units),unlimited=.TRUE.)

        ! Write grid information
        call rembo_grid_write(dom%grid,filename,dom%par%domain,dom%par%grid_name,create=.FALSE.)

        return

    end subroutine rembo_write_init
    
    subroutine rembo_emb_write_init(emb,filename,domain,grid_name,time_init,units)

        implicit none 

        type(diffusion_class), intent(IN) :: emb 
        character(len=*),  intent(IN) :: filename 
        character(len=*),  intent(IN) :: domain
        character(len=*),  intent(IN) :: grid_name
        real(wp),          intent(IN) :: time_init
        character(len=*),  intent(IN) :: units

        ! Local variables 

        ! Initialize netcdf file and dimensions
        call nc_create(filename)
        call nc_write_dim(filename,"xc",    x=emb%grid%xc,  units="kilometers")
        call nc_write_dim(filename,"yc",    x=emb%grid%yc,  units="kilometers")
        call nc_write_dim(filename,"month", x=1,dx=1,nx=12, units="month")
        call nc_write_dim(filename,"time",  x=time_init,dx=1.0_wp,nx=1,units=trim(units),unlimited=.TRUE.)

        ! Write grid information
        call rembo_grid_write(emb%grid,filename,domain,grid_name,create=.FALSE.)

        return

    end subroutine rembo_emb_write_init
    
    subroutine rembo_emb_write_step_grid(emb,filename,m)

        implicit none 

        type(diffusion_class) :: emb 
        character(len=*)  :: filename 
        integer :: nx, ny, m 

        nx = emb%grid%nx
        ny = emb%grid%ny
        
        if (m .eq. 1) then 
            call nc_write(filename,"z_srf",emb%z_srf,dim1="xc",dim2="yc")
            call nc_write(filename,"mask", emb%mask, dim1="xc",dim2="yc")
        end if 

        call nc_write(filename,"kappa", real(emb%kappa),dim1="xc",dim2="yc",dim3="month", &
                      start=[1,1,m],count=[nx,ny,1])
        call nc_write(filename,"tsl_bnd", real(emb%tsl_bnd),dim1="xc",dim2="yc",dim3="month",&
                      start=[1,1,m],count=[nx,ny,1])
        call nc_write(filename,"tsl",real(emb%tsl),dim1="xc",dim2="yc",dim3="month",&
                      start=[1,1,m],count=[nx,ny,1])
        call nc_write(filename,"tsl_F",real(emb%tsl_F),dim1="xc",dim2="yc",dim3="month",&
                      start=[1,1,m],count=[nx,ny,1])
        call nc_write(filename,"kappaw", real(emb%kappaw),dim1="xc",dim2="yc",dim3="month", &
                      start=[1,1,m],count=[nx,ny,1])
        call nc_write(filename,"ccw_bnd",real(emb%ccw_bnd),dim1="xc",dim2="yc",dim3="month",&
                      start=[1,1,m],count=[nx,ny,1])
        call nc_write(filename,"ccw",real(emb%ccw),dim1="xc",dim2="yc",dim3="month",&
                      start=[1,1,m],count=[nx,ny,1])
        call nc_write(filename,"ccw_F",real(emb%ccw_F),dim1="xc",dim2="yc",dim3="month",&
                      start=[1,1,m],count=[nx,ny,1])
        call nc_write(filename,"tcw_bnd",real(emb%tcw_bnd),dim1="xc",dim2="yc",dim3="month",&
                      start=[1,1,m],count=[nx,ny,1])
        call nc_write(filename,"tcw",real(emb%tcw),dim1="xc",dim2="yc",dim3="month",&
                      start=[1,1,m],count=[nx,ny,1])
        
        call nc_write(filename,"ug",real(emb%ug),dim1="xc",dim2="yc",dim3="month",&
                      start=[1,1,m],count=[nx,ny,1])
        call nc_write(filename,"vg",real(emb%vg),dim1="xc",dim2="yc",dim3="month",&
                      start=[1,1,m],count=[nx,ny,1])
        
        return 

    end subroutine rembo_emb_write_step_grid

end module rembo_api