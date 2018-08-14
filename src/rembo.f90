module rembo
    ! Wrapper to hold all modules needed for librembo 
       
    use rembo_defs 
    use rembo_atm 
    use rembo_physics 

    use nml 
    use ncio 
    use solvers 
    use insolation 

    implicit none 

    ! Day count for the middle of each month of the year
    integer, parameter :: mdays(12) =[30,60,90,120,150,180,210,240,270,300,330,360]-15

    private 
    public :: rembo_update
    public :: rembo_init 
    public :: rembo_end 
    public :: rembo_write_init 

contains 

    subroutine rembo_update(dom,z_srf,f_ice,f_shlf,t2m,Z,co2_a,year)
        ! Calculate atmosphere for each month of the year 

        implicit none 

        type(rembo_class), intent(INOUT) :: dom 

        real(wp), intent(IN) :: z_srf(:,:)      ! [m]     Surface elevation
        real(wp), intent(IN) :: f_ice(:,:)      ! [--]    Fraction of land-based ice coverage in cell
        real(wp), intent(IN) :: f_shlf(:,:)     ! [--]    Fraction of floating (shelf) ice coverage in cell
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

        ! Calculate the coriolis parameter for the current gridpoints
        dom%bnd%f = calc_coriolis(real(dom%grid%lat,wp))

        ! Calculate surface gradients and total magnitude
        call d_dx(dom%bnd%dzsdx,dom%bnd%z_srf,dx=real(dom%grid%G%dx*dom%grid%xy_conv,wp))
        call d_dy(dom%bnd%dzsdy,dom%bnd%z_srf,dx=real(dom%grid%G%dy*dom%grid%xy_conv,wp))
        dom%bnd%dzsdxy = calc_magnitude(dom%bnd%dzsdx,dom%bnd%dzsdy)

        ! Calculate the rembo relaxation mask
        dom%bnd%mask = 0
        where(dom%bnd%z_srf .gt. 0.0) dom%bnd%mask = 1 

        dom%bnd%mask = gen_relaxation(dom%bnd%z_srf,real(dom%grid%x,wp),real(dom%grid%y,wp),radius=20.0)  
        
        ! EMB OUTPUT FOR TESTING 
        call rembo_emb_write_init(dom%emb,"test.nc",time_init=real(year,wp),units="kyr ago")

        ! Loop over each month, calculate rembo atmosphere 
        do m = 1, nm 

            ! Determine representative day of the month
            day = mdays(m)

            ! == Store monthly boundary variables =====

            ! Calculate representative insolation for the month
            dom%now%S = calc_insol_day(day,dom%grid%lat,dble(year),fldr="libs/insol/input")

            ! Save all other boundary variables 
            dom%now%t2m_bnd = t2m(:,:,m) 
            dom%now%co2_a   = co2_a 
            dom%now%Z       = Z(:,:,m)

            ! == Calculate monthly derived boundary variables =====

            ! Calc gradient of Z: dZdx, dZdy 
            call d_dx(dom%now%dZdx,dom%now%Z,dx=real(dom%grid%G%dx*dom%grid%xy_conv,wp))
            call d_dy(dom%now%dZdy,dom%now%Z,dx=real(dom%grid%G%dy*dom%grid%xy_conv,wp))
            ! ajr: to do: move this inside of rembo_calc_atmosphere using local variables
            
            ! == Calculate rembo atmosphere =====

            !dom%now%t2m = dom%now%t2m_bnd
            
            call rembo_calc_atmosphere(dom%now,dom%emb,dom%bnd,dom%par,day,year)
            
            ! Print summary 
            call rembo_print(dom,m,day,year)

            ! Store data in monthly object 
            dom%mon(m) = dom%now

            ! EMB OUTPUT FOR TESTING 
            call rembo_emb_write_step_grid(dom%emb,"test.nc",m)

        end do 

        return 

    end subroutine rembo_update 
    
    subroutine rembo_init(dom,par_path,domain,grid)

        use solvers

        implicit none
        
        type(rembo_class), intent(INOUT) :: dom
        character(len=*),  intent(IN)    :: par_path  
        character(len=*),  intent(IN)    :: domain 
        type(grid_class),  intent(IN)    :: grid 

        ! Local variables         
        character(len=256) :: filename, file_boundary 
        integer            :: i, m, nm 
        real(wp)           :: dt_check(2) 
        character(len=64)  :: fmt1 

        ! Load the rembo parameters
        call rembo_par_load(dom%par,trim(par_path),domain)

        ! Initialize rembo domain and grid
        dom%grid       = grid 
        dom%par%npts   = grid%npts 
        dom%par%nx     = grid%G%nx 
        dom%par%ny     = grid%G%ny 
        
        nm = size(dom%mon)

        ! Allocate boundary variables
        call rembo_bnd_alloc(dom%bnd,dom%par%nx,dom%par%ny)
        
        ! Allocate state variables
        call rembo_alloc(dom%now,dom%par%nx,dom%par%ny)
        do m = 1, nm 
            call rembo_alloc(dom%mon(m),dom%par%nx,dom%par%ny)
        end do 
        call rembo_alloc(dom%ann,dom%par%nx,dom%par%ny)

        ! Diffusion grid and variables
        call rembo_emb_init(dom%emb,grid,dx=dom%par%emb_dx)
        
        write(*,*) "rembo_init :: allocated rembo variables."

        ! Perform mapping between emb and rembo grids
        call map_init(dom%emb%map_toemb,  dom%grid,dom%emb%grid,max_neighbors=10,load=.TRUE.)
        call map_init(dom%emb%map_fromemb,dom%emb%grid,dom%grid,max_neighbors=10,load=.TRUE.)
        
        ! Check diffusion time step consistency
        dt_check(1) = diff2D_timestep(real(dom%emb%grid%G%dx*dom%emb%grid%xy_conv,wp), &
                                      real(dom%emb%grid%G%dy*dom%emb%grid%xy_conv,wp), &
                                      dom%par%en_D)
        dt_check(2) = diff2D_timestep(real(dom%emb%grid%G%dx*dom%emb%grid%xy_conv,wp), &
                                      real(dom%emb%grid%G%dy*dom%emb%grid%xy_conv,wp), &
                                      dom%par%ccw_D)
!         dt_check(1) = diff2Dadi_timestep(dom%emb%grid%G%dx*dom%emb%grid%xy_conv, &
!                                          dom%emb%grid%G%dy*dom%emb%grid%xy_conv, &
!                                          dom%par%en_D)
!         dt_check(2) = diff2Dadi_timestep(dom%emb%grid%G%dx*dom%emb%grid%xy_conv, &
!                                          dom%emb%grid%G%dy*dom%emb%grid%xy_conv, &
!                                          dom%par%ccw_D)
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

    subroutine rembo_emb_init(emb,grid,dx)
        ! Initialize the diffusion variables on
        ! the desired grid resolution.
        ! grid = original grid definition of rembo
        ! emb%grid = low resolution grid of resolution dx

        implicit none

        type(diffusion_class), intent(INOUT) :: emb 
        type(grid_class),      intent(IN)    :: grid 
        real(wp),              intent(IN)    :: dx 
        
        integer :: nx, ny 
        character(len=256) :: name 
        real(wp), allocatable :: x(:)
        real(wp), allocatable :: y(:) 

        if (dx .ge. 100.d0) then 
            write(name,"(a,i3,a)") trim(grid%name)//"-emb", int(dx), "km"
        else if (dx .ge. 10.d0) then 
            write(name,"(a,i2,a)") trim(grid%name)//"-emb", int(dx), "km"
        else
            write(*,*) "rembo_emb_init:: Error: ", &
            "Such a high resolution grid for diffusion is not supported."
            write(*,*) "emb dx =",dx 
            stop 
        end if 

        ! Initialize the diffusion grid
        call grid_init(emb%grid,grid,name, &
                       x0=minval(grid%x)-dx,dx=real(dx,dp),nx=ceiling((maxval(grid%x)-minval(grid%x))/dx)+2, &
                       y0=minval(grid%y)-dx,dy=real(dx,dp),ny=ceiling((maxval(grid%y)-minval(grid%y))/dx)+2)

        nx = emb%grid%G%nx
        ny = emb%grid%G%ny 

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
        allocate(emb%en(nx,ny))
        allocate(emb%en_bnd(nx,ny))
        allocate(emb%en_F(nx,ny))

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
        emb%en          = 0.0
        emb%en_bnd      = 0.0
        emb%en_F        = 0.0

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
                sum(dom%now%c_w*sec_day,mask=dom%bnd%mask>0)/npts, &     ! [mm/d]
                sum(dom%now%pr*sec_day, mask=dom%bnd%mask>0)/npts        ! [mm/d]

        else 
            ! Print empty table 
            write(*,fmt_table) year, d, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 

        end if 

        return 

    end subroutine rembo_print 
    
    subroutine rembo_par_load(par,filename,domain)

        type(rembo_param_class), intent(INOUT) :: par 
        character(len=*),        intent(IN)    :: filename 
        character(len=*),        intent(IN)    :: domain 

        par%domain = trim(domain)
        call nml_read(filename,"rembo_params","restart",    par%restart    )

        call nml_read(filename,"rembo_params","H_e",  par%H_e   )
        call nml_read(filename,"rembo_params","emb_dx",  par%emb_dx  )
        call nml_read(filename,"rembo_params","rembo1",  par%rembo1  )
        call nml_read(filename,"rembo_params","dist_rel",  par%dist_rel  )
        call nml_read(filename,"rembo_params","en_dt",  par%en_dt  )
        call nml_read(filename,"rembo_params","en_D",  par%en_D  )
        call nml_read(filename,"rembo_params","en_D_win",  par%en_D_win  )
        call nml_read(filename,"rembo_params","en_D_sum",  par%en_D_sum  )
        call nml_read(filename,"rembo_params","en_kr",  par%en_kr  )
        call nml_read(filename,"rembo_params","en_kz",  par%en_kz  )
        call nml_read(filename,"rembo_params","en_kl",  par%en_kl  )
        call nml_read(filename,"rembo_params","en_kdT",  par%en_kdT  )
        call nml_read(filename,"rembo_params","en_Ha",  par%en_Ha  )
        call nml_read(filename,"rembo_params","ccw_dt",  par%ccw_dt  )
        call nml_read(filename,"rembo_params","ccw_D",  par%ccw_D  )
        call nml_read(filename,"rembo_params","ccw_kr",  par%ccw_kr  )
        call nml_read(filename,"rembo_params","k_c",  par%k_c  )
        call nml_read(filename,"rembo_params","k_p",  par%k_p  )
        call nml_read(filename,"rembo_params","k_z",  par%k_z  )
        call nml_read(filename,"rembo_params","k_x",  par%k_x  )
        call nml_read(filename,"rembo_params","e0",  par%e0  )
        call nml_read(filename,"rembo_params","c1",  par%c1  )
        call nml_read(filename,"rembo_params","k_p_now",  par%k_p_now  )
        call nml_read(filename,"rembo_params","k_t",  par%k_t  )
        call nml_read(filename,"rembo_params","k_w",  par%k_w  )
        call nml_read(filename,"rembo_params","sf_a",  par%sf_a  )
        call nml_read(filename,"rembo_params","sf_b",  par%sf_b  )
        call nml_read(filename,"rembo_params","f_k",  par%f_k  )
        call nml_read(filename,"rembo_params","nk1",  par%nk1  )
        call nml_read(filename,"rembo_params","nk2",  par%nk2  )
        call nml_read(filename,"rembo_params","nk3",  par%nk3  )
        call nml_read(filename,"rembo_params","gamma",  par%gamma  )
        call nml_read(filename,"rembo_params","gamma2",  par%gamma2  )
        call nml_read(filename,"rembo_params","S0",  par%S0  )
        call nml_read(filename,"rembo_params","Lw",  par%Lw  )
        call nml_read(filename,"rembo_params","Lm",  par%Lm  )
        call nml_read(filename,"rembo_params","Ls",  par%Ls  )
        call nml_read(filename,"rembo_params","rho_w",  par%rho_w  )
        call nml_read(filename,"rembo_params","cp",  par%cp  )
        call nml_read(filename,"rembo_params","ci",  par%ci  )
        call nml_read(filename,"rembo_params","alp_a",  par%alp_a  )
        call nml_read(filename,"rembo_params","alp_b",  par%alp_b  )
        
        call nml_read(filename,"rembo_params","alp_c",  par%alp_c  )
        call nml_read(filename,"rembo_params","lwu_a",  par%lwu_a  )
        call nml_read(filename,"rembo_params","lwu_b",  par%lwu_b  )
        call nml_read(filename,"rembo_params","lwu_c",  par%lwu_c  )
        call nml_read(filename,"rembo_params","swds_a",  par%swds_a  )
        call nml_read(filename,"rembo_params","swds_b",  par%swds_b  )
        call nml_read(filename,"rembo_params","swds_c",  par%swds_c  )
        call nml_read(filename,"rembo_params","lwds_a",  par%lwds_a  )
        call nml_read(filename,"rembo_params","lwds_b",  par%lwds_b  )
        call nml_read(filename,"rembo_params","lwds_c",  par%lwds_c  )
        call nml_read(filename,"rembo_params","lwds_d",  par%lwds_d  )
        
        call nml_read(filename,"rembo_params","shfs_Cm",  par%shfs_Cm)
        call nml_read(filename,"rembo_params","shfs_p",  par%shfs_p  )

        call nml_read(filename,"rembo1_params","ta",  par%r1_ta  )
        call nml_read(filename,"rembo1_params","tb",  par%r1_tb  )

        ! How many time steps in 1 day, aprx?
        par%en_nstep  = floor(sec_day / par%en_dt)
        par%ccw_nstep = floor(sec_day / par%ccw_dt)
        
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
        call nc_write_dim(filename,"xc",    x=dom%grid%G%x,  units="kilometers")
        call nc_write_dim(filename,"yc",    x=dom%grid%G%y,  units="kilometers")
        call nc_write_dim(filename,"month", x=1,dx=1,nx=12,   units="month")
        call nc_write_dim(filename,"time",  x=time_init,dx=1.0_wp,nx=1,units=trim(units),unlimited=.TRUE.)

        ! Write grid information
        call grid_write(dom%grid,filename,xnm="xc",ynm="yc",create=.FALSE.)

        return

    end subroutine rembo_write_init 
    
    subroutine rembo_emb_write_init(emb,filename,time_init,units)

        implicit none 

        type(diffusion_class), intent(IN) :: emb 
        character(len=*),  intent(IN) :: filename 
        real(wp),          intent(IN) :: time_init
        character(len=*),  intent(IN) :: units

        ! Local variables 

        ! Initialize netcdf file and dimensions
        call nc_create(filename)
        call nc_write_dim(filename,"xc",    x=emb%grid%G%x,  units="kilometers")
        call nc_write_dim(filename,"yc",    x=emb%grid%G%y,  units="kilometers")
        call nc_write_dim(filename,"month", x=1,dx=1,nx=12,   units="month")
        call nc_write_dim(filename,"time",  x=time_init,dx=1.0_wp,nx=1,units=trim(units),unlimited=.TRUE.)

        ! Write grid information
        call grid_write(emb%grid,filename,xnm="xc",ynm="yc",create=.FALSE.)

        return

    end subroutine rembo_emb_write_init 
    
    subroutine rembo_emb_write_step_grid(emb,filename,m)

        implicit none 

        type(diffusion_class) :: emb 
        character(len=*)  :: filename 
        integer :: nx, ny, m 

        nx = emb%grid%G%nx
        ny = emb%grid%G%ny
        
        if (m .eq. 1) then 
            call nc_write(filename,"z_srf",emb%z_srf,dim1="xc",dim2="yc")
            call nc_write(filename,"mask", emb%mask, dim1="xc",dim2="yc")
        end if 

        call nc_write(filename,"kappa", real(emb%kappa),dim1="xc",dim2="yc",dim3="month", &
                      start=[1,1,m],count=[nx,ny,1])
        call nc_write(filename,"kappaw", real(emb%kappaw),dim1="xc",dim2="yc",dim3="month", &
                      start=[1,1,m],count=[nx,ny,1])
        call nc_write(filename,"tsl_bnd", real(emb%tsl_bnd),dim1="xc",dim2="yc",dim3="month",&
                      start=[1,1,m],count=[nx,ny,1])
        call nc_write(filename,"ccw_bnd",real(emb%ccw_bnd),dim1="xc",dim2="yc",dim3="month",&
                      start=[1,1,m],count=[nx,ny,1])
        call nc_write(filename,"tcw_bnd",real(emb%tcw_bnd),dim1="xc",dim2="yc",dim3="month",&
                      start=[1,1,m],count=[nx,ny,1])
        call nc_write(filename,"tsl",real(emb%tsl),dim1="xc",dim2="yc",dim3="month",&
                      start=[1,1,m],count=[nx,ny,1])
        call nc_write(filename,"ccw",real(emb%ccw),dim1="xc",dim2="yc",dim3="month",&
                      start=[1,1,m],count=[nx,ny,1])
        call nc_write(filename,"ccw_F",real(emb%ccw_F),dim1="xc",dim2="yc",dim3="month",&
                      start=[1,1,m],count=[nx,ny,1])
        call nc_write(filename,"tcw",real(emb%tcw),dim1="xc",dim2="yc",dim3="month",&
                      start=[1,1,m],count=[nx,ny,1])
        call nc_write(filename,"en",real(emb%en),dim1="xc",dim2="yc",dim3="month",&
                      start=[1,1,m],count=[nx,ny,1])
        call nc_write(filename,"en_bnd",real(emb%en_bnd),dim1="xc",dim2="yc",dim3="month",&
                      start=[1,1,m],count=[nx,ny,1])
        call nc_write(filename,"en_F",real(emb%en_F),dim1="xc",dim2="yc",dim3="month",&
                      start=[1,1,m],count=[nx,ny,1])
            
        return 

    end subroutine rembo_emb_write_step_grid

end module rembo 