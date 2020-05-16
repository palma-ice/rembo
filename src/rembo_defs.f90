module rembo_defs 

    use nml 
    use coord 

    implicit none 

    ! =========================================================================
    !
    ! CONSTANTS (program precision, global constants)
    !
    ! =========================================================================

    ! Internal constants
    integer,  parameter :: dp  = kind(1.d0)
    integer,  parameter :: sp  = kind(1.0)

    ! Choose the precision of the library (sp,dp)
    integer,  parameter :: wp = sp 

    ! Write flags 
    logical, parameter :: rembo_write_log = .TRUE. 

    ! Missing value and aliases
    real(wp), parameter :: MISSING_VALUE_DEFAULT = real(-9999.0,wp)
    real(wp), parameter :: MISSING_VALUE = MISSING_VALUE_DEFAULT
    real(wp), parameter :: MV = MISSING_VALUE_DEFAULT
    
    ! Error distance (very large), error index, and smallest number epsilon 
    real(wp), parameter :: ERR_DIST = real(1E8,wp) 
    integer,    parameter :: ERR_IND  = -1 
    real(wp), parameter :: eps      = real(1E-8,wp) 
    
    ! Mathematical constants
    real(wp), parameter :: pi  = real(2._dp*acos(0.0_dp),wp)
    real(wp), parameter :: degrees_to_radians = real(pi / 180._dp,wp)  ! Conversion factor between radians and degrees
    real(wp), parameter :: radians_to_degrees = real(180._dp / pi,wp)  ! Conversion factor between degrees and radians
    
    logical :: rembo_use_omp 

    integer,  parameter :: nd  = 360
    integer,  parameter :: nm  = 12
    integer,  parameter :: ndm = 30

    ! The constants below should be loaded using the global subroutine
    ! defined below `rembo_constants_load`.
    ! Note: The key limitation imposed by defining the parameters defined 
    ! globally is that these constants must be the same for all domains 
    ! being run in the same program. 

    ! Physical constants
    real(wp) :: sec_year       ! [s] seconds per year 
    real(wp) :: g              ! Gravitational accel.  [m s-2]
    real(wp) :: omega          ! Coriolis constant [omega = 7.2921d-5]
    real(wp) :: T0             ! Reference freezing temperature [K] 
    real(wp) :: rho_ice        ! Density ice           [kg m-3] 
    real(wp) :: rho_w          ! Density water         [kg m-3] 

    real(wp), parameter :: sec_day0 = 8.64e4_wp 
    real(wp) :: sec_day 
    real(wp) :: sec_frac 

    ! Time conversions
    real(wp), parameter :: day_year   = real(nd,wp)
    real(wp), parameter :: day_month  = real(ndm,wp)
    real(wp), parameter :: month_year = real(nm,wp)


    ! First define all parameters needed to represent a given domain
    type rembo_param_class

        character(len=256)  :: domain
        character(len=256)  :: grid_name
        character(len=256)  :: restart 
        integer             :: npts, nx, ny
        real(wp)            :: dx 

        ! Physics
        real(wp)    :: H_e   ! Precip vapor scale height (m)
        logical     :: rembo1  
        real(wp)    :: emb_dx 
        real(wp)    :: dist_rel
        real(wp)    :: en_dt, en_D, en_kr, en_kz, en_kl, en_kdT, en_Ha 
        real(wp)    :: en_D_win, en_D_sum 
        real(wp)    :: ccw_dt, ccw_D, ccw_kr 
        integer     :: en_nstep, ccw_nstep 
        real(wp)    :: k_c, k_p, k_z, k_x, e0, c1
        real(wp)    :: k_c_sec 
        real(wp)    :: k_p_now, k_t, k_w, sf_a, sf_b 
        real(wp)    :: f_k 
        real(wp)    :: nk1, nk2, nk3 
        real(wp)    :: gamma, gamma2 
        real(wp)    :: S0, Lw, Lm, Ls
        real(wp)    :: rho_w   
        real(wp)    :: cp, ci
        real(wp)    :: alp_a, alp_b, alp_c, lwu_a, lwu_b, lwu_c
        real(wp)    :: lwds_a, lwds_b, lwds_c, lwds_d
        real(wp)    :: swds_a, swds_b, swds_c
        real(wp)    :: shfs_Cm, shfs_p   
        real(wp)    :: r1_ta, r1_tb 

    end type

    type rembo_boundary_class

        ! Annual boundary variables
        real(wp), allocatable :: z_srf(:,:)    ! [m]     Surface elevation
        real(wp), allocatable :: f_ice(:,:)    ! [--]    Fraction of land-based ice coverage in cell
        real(wp), allocatable :: f_shlf(:,:)   ! [--]    Fraction of floating (shelf) ice coverage in cell
        
        ! Derived annual boundary variables
        integer,  allocatable :: mask(:,:)     ! [--]    0: Ocean; 1: Land, 2: Grounded ice, 3: Floating ice
        real(wp), allocatable :: f(:,:)        ! [--]    Coriolis parameter
        real(wp), allocatable :: dzsdx(:,:)
        real(wp), allocatable :: dzsdy(:,:)
        real(wp), allocatable :: dzsdxy(:,:)

    end type 

    ! Now define all variables of the domain
    type rembo_state_class

        ! Monthly forcing variables 
        real(wp), allocatable :: S(:,:)        ! [W m-2] Insolation top-of-atmosphere
        real(wp), allocatable :: t2m_bnd(:,:)  ! [K]     Near-surface temperature (used for boundary)
        real(wp), allocatable :: al_s(:,:)     ! [--]    Surface albedo 
        real(wp), allocatable :: co2_a(:,:)    ! [ppm]   Atmospheric CO2 concentration
        real(wp), allocatable :: Z(:,:)        ! [m?]    Geopotential height of 750 Mb layer
        real(wp), allocatable :: dZdx(:,:)
        real(wp), allocatable :: dZdy(:,:)
        
        ! Annual variables 
        real(wp), allocatable :: rco2_a(:,:)
        real(wp), allocatable :: rho_a(:,:)
        real(wp), allocatable :: sp(:,:)
        
        ! Seasonal variables
        real(wp), allocatable :: gamma(:,:)
        real(wp), allocatable :: t2m(:,:)   
        real(wp), allocatable :: ct2m(:,:)
        real(wp), allocatable :: tsurf(:,:)   
        real(wp), allocatable :: pr(:,:)
        real(wp), allocatable :: sf(:,:)
        real(wp), allocatable :: q_s(:,:)
        real(wp), allocatable :: q_sat(:,:)
        real(wp), allocatable :: q_r(:,:)
        real(wp), allocatable :: tcw(:,:), tcw_sat(:,:)
        real(wp), allocatable :: ccw(:,:), c_w(:,:), ccw_prev(:,:) 
        real(wp), allocatable :: ug(:,:), vg(:,:), uvg(:,:), ww(:,:), cc(:,:)
        real(wp), allocatable :: swd(:,:), lwu(:,:), al_p(:,:), at(:,:)
        real(wp), allocatable :: swd_s(:,:), lwd_s(:,:), shf_s(:,:), lhf_s(:,:), lwu_s(:,:)
        real(wp), allocatable :: u_s(:,:), v_s(:,:), uv_s(:,:)  
        
        real(wp), allocatable :: u_k(:,:), v_k(:,:), uv_k(:,:), dtsldx(:,:), dtsldy(:,:), dtsldxy(:,:)
    end type 

    ! Define all variables needed for diffusion on lo-res grid
    type diffusion_class

        type(grid_class)      :: grid ! EMB diffusion resolution grid
        type(map_class)       :: map_toemb, map_fromemb  ! map EMB => rembo grid

        ! Relaxation mask, topography  
        integer,  allocatable :: mask(:,:)
        real(wp), allocatable :: z_srf(:,:), rho_a(:,:) 
        real(wp), allocatable :: dzsdx(:,:), dzsdy(:,:), dzsdxy(:,:) 

        ! Energy and moisture balance variables
        real(wp), allocatable :: tsl(:,:), tsl_bnd(:,:), tsl_F(:,:)  
        real(wp), allocatable :: ccw(:,:), ccw_bnd(:,:), ccw_F(:,:) 
        real(wp), allocatable :: ccw_cw(:,:), ccw_pr(:,:) 
        real(wp), allocatable :: tcw(:,:), tcw_bnd(:,:)
        real(wp), allocatable :: ug(:,:), vg(:,:), uvg(:,:), ww(:,:), q_r(:,:)  

        ! Diffusion 
        real(wp), allocatable :: kappa(:,:), kappaw(:,:) 
        real(wp) :: en_dt, tsl_fac, en_kr, en_kz
        integer :: en_nstep
        logical :: bnd_pr 
    end type 

    type rembo_class

        type(rembo_param_class)   :: par        ! physical parameters
        type(grid_class)          :: grid       ! Grid definition   (from coordinates module)
        
        ! Boundary variables
        type(rembo_boundary_class) :: bnd 

        ! current variables, month and various averages
        type(rembo_state_class) :: now, mon(12), ann

        ! Variables and grid definitions for energy-moisture balance calculations
        type(diffusion_class) :: emb

    end type


    public   ! All rembo defs are public

contains 

    function rembo_get_working_precision() result(rembo_prec)

        implicit none 

        integer :: rembo_prec 

        rembo_prec = kind(wp)

        return 

    end function rembo_get_working_precision

        
    subroutine rembo_parse_path(path,domain,grid_name)

        implicit none 

        character(len=*), intent(INOUT) :: path 
        character(len=*), intent(IN)    :: domain, grid_name 

        call nml_replace(path,"{domain}",   trim(domain))
        call nml_replace(path,"{grid_name}",trim(grid_name))
        
        return 

    end subroutine rembo_parse_path

    subroutine rembo_global_init(filename)

        !$ use omp_lib 

        character(len=*), intent(IN)  :: filename
        
        ! Local variables
        logical :: init_pars 
        integer :: n_threads 
        character(len=10) :: n_threads_str 

        init_pars = .TRUE. 
        
        ! Check openmp status - set global variable to use as a switch 
        rembo_use_omp = .FALSE. 
        !$ rembo_use_omp = .TRUE.

        ! Output some information about openmp status 
        if (rembo_use_omp) then 
            
            n_threads = 1
            !$ n_threads = omp_get_max_threads() 

            write(n_threads_str,"(i10)") n_threads 
            n_threads_str = adjustl(n_threads_str)

            write(*,*) "rembo_global_init:: openmp is active, rembo will run on "//trim(n_threads_str)//" thread(s)."
            
        else 
            
            n_threads = 1
            write(*,*) "rembo_global_init:: openmp is not active, rembo will run on 1 thread."

        end if 

        ! Store parameter values in output object
        call nml_read(filename,"rembo_constants","sec_year",    sec_year,   init=init_pars)
        call nml_read(filename,"rembo_constants","g",           g,          init=init_pars)
        call nml_read(filename,"rembo_constants","omega",       omega,      init=init_pars)
        call nml_read(filename,"rembo_constants","T0",          T0,         init=init_pars)
        
        call nml_read(filename,"rembo_constants","rho_ice",     rho_ice,    init=init_pars)
        call nml_read(filename,"rembo_constants","rho_w",       rho_w,      init=init_pars)

        if (rembo_write_log) then 
            write(*,*) "yelmo:: loaded global constants:"
            write(*,*) "    sec_year = ", sec_year 
            write(*,*) "    g        = ", g 
            write(*,*) "    omega    = ", omega
            write(*,*) "    T0       = ", T0 
            write(*,*) "    rho_ice  = ", rho_ice 
            write(*,*) "    rho_w    = ", rho_w 

        end if 

        sec_day    = sec_year / day_year   ! 8.765813e4
        sec_frac   = sec_day / sec_day0

        return

    end subroutine rembo_global_init
    
    subroutine rembo_grid_define(grid,grid_name,grid_in)

        implicit none 

        type(grid_class), intent(OUT)   :: grid  
        character(len=*), intent(INOUT) :: grid_name      ! Overwritable if grid_in is present
        type(grid_class), intent(IN), optional :: grid_in 

        if (present(grid_in)) then 

            grid = grid_in 

            ! Ensure parameter grid_name is consistent with defined grid 
            grid_name = grid%name 
        
        else 
            ! Define rembo grid from predefined options 

            select case(trim(grid_name))

                ! Note - all North projections now use the ESPG-3413
                ! polar stereographic projection with (lambda=-45.d0,phi=70.d0)
                ! Smaller Northern domains like Eurasia and Greenland use
                ! the same projection for consistency. 
                ! ESPG-3413 (lambda=-45.d0,phi=70.d0) is used for Greenland in 
                ! model intercomparison exercises, eg ISMIP6. 

                ! NORTH DOMAINS ======================= 

                case("NH-40KM")
                    call grid_init(grid,name="NH-40KM",mtype="polar_stereographic",units="kilometers", &
                            lon180=.TRUE.,x0=-4900.d0,dx=40.0d0,nx=221,y0=-5400.d0,dy=40.0d0,ny=221, &
                            lambda=-45.d0,phi=70.d0)
                
                case("NH-20KM")
                    call grid_init(grid,name="NH-20KM",mtype="polar_stereographic",units="kilometers", &
                            lon180=.TRUE.,x0=-4900.d0,dx=20.0d0,nx=441,y0=-5400.d0,dy=20.0d0,ny=441, &
                            lambda=-45.d0,phi=70.d0)
                
                case("NH-10KM")
                    call grid_init(grid,name="NH-10KM",mtype="polar_stereographic",units="kilometers", &
                            lon180=.TRUE.,x0=-4900.d0,dx=10.0d0,nx=881,y0=-5400.d0,dy=10.0d0,ny=881, &
                            lambda=-45.d0,phi=70.d0)
                
                case("NH-5KM")
                    call grid_init(grid,name="NH-5KM",mtype="polar_stereographic",units="kilometers", &
                            lon180=.TRUE.,x0=-4900.d0,dx=5.0d0,nx=1761,y0=-5400.d0,dy=5.0d0,ny=1761, &
                            lambda=-45.d0,phi=70.d0)
            
                ! EURASIA DOMAINS ======================= 

                case("EIS-40KM")
                    call grid_init(grid,name="EIS-40KM",mtype="polar_stereographic",units="kilometers", &
                            lon180=.TRUE.,x0=380.d0,dx=40.0d0,nx=89,y0=-5000.d0,dy=40.0d0,ny=161, &
                            lambda=-45.d0,phi=70.d0)
                
                case("EIS-20KM")
                    call grid_init(grid,name="EIS-20KM",mtype="polar_stereographic",units="kilometers", &
                            lon180=.TRUE.,x0=380.d0,dx=20.0d0,nx=177,y0=-5000.d0,dy=20.0d0,ny=321, &
                            lambda=-45.d0,phi=70.d0)
                
                case("EIS-10KM")
                    call grid_init(grid,name="EIS-10KM",mtype="polar_stereographic",units="kilometers", &
                            lon180=.TRUE.,x0=380.d0,dx=10.0d0,nx=353,y0=-5000.d0,dy=10.0d0,ny=641, &
                            lambda=-45.d0,phi=70.d0)
                
                case("EIS-5KM")
                    call grid_init(grid,name="EIS-5KM",mtype="polar_stereographic",units="kilometers", &
                            lon180=.TRUE.,x0=380.d0,dx=5.0d0,nx=705,y0=-5000.d0,dy=5.0d0,ny=1281, &
                            lambda=-45.d0,phi=70.d0)
                    
                ! GREENLAND DOMAINS =======================

                case("GRL-40KM")
                    call grid_init(grid,name="GRL-40KM",mtype="polar_stereographic",units="kilometers", &
                            lon180=.TRUE.,x0=-720.d0,dx=40.0d0,nx=43,y0=-3450.d0,dy=40.0d0,ny=73, &
                            lambda=-45.d0,phi=70.d0)
                
                case("GRL-20KM")
                    call grid_init(grid,name="GRL-20KM",mtype="polar_stereographic",units="kilometers", &
                            lon180=.TRUE.,x0=-720.d0,dx=20.0d0,nx=85,y0=-3450.d0,dy=20.0d0,ny=145, &
                            lambda=-45.d0,phi=70.d0)
                
                case("GRL-10KM")
                    call grid_init(grid,name="GRL-10KM",mtype="polar_stereographic",units="kilometers", &
                            lon180=.TRUE.,x0=-720.d0,dx=10.0d0,nx=169,y0=-3450.d0,dy=10.0d0,ny=289, &
                            lambda=-45.d0,phi=70.d0)
                
                case("GRL-5KM")
                    call grid_init(grid,name="GRL-5KM",mtype="polar_stereographic",units="kilometers", &
                            lon180=.TRUE.,x0=-720.d0,dx=5.0d0,nx=337,y0=-3450.d0,dy=5.0d0,ny=577, &
                            lambda=-45.d0,phi=70.d0)
                
                case("GRL-2KM")
                    call grid_init(grid,name="GRL-2KM",mtype="polar_stereographic",units="kilometers", &
                            lon180=.TRUE.,x0=-720.d0,dx=2.0d0,nx=841,y0=-3450.d0,dy=2.0d0,ny=1441, &
                            lambda=-45.d0,phi=70.d0)
                
                case("GRL-1KM")
                    call grid_init(grid,name="GRL-1KM",mtype="polar_stereographic",units="kilometers", &
                            lon180=.TRUE.,x0=-720.d0,dx=1.0d0,nx=1681,y0=-3450.d0,dy=1.0d0,ny=2881, &
                            lambda=-45.d0,phi=70.d0)

                case("Bamber01-20KM")
                    call grid_init(grid,name="Bamber01-20KM",mtype="polar_stereographic",units="kilometers", &
                            lon180=.TRUE.,x0=-800.d0,dx=20.d0,nx=76,y0=-3400.d0,dy=20.d0,ny=141, &
                            lambda=-39.d0,phi=90.d0)

                case("Bamber01-10KM")
                    call grid_init(grid,name="Bamber01-10KM",mtype="polar_stereographic",units="kilometers", &
                            lon180=.TRUE.,x0=-800.d0,dx=10.d0,nx=151,y0=-3400.d0,dy=10.d0,ny=281, &
                            lambda=-39.d0,phi=90.d0)

                ! ANTARCTICA DOMAINS ======================= 

                case("ANT-80KM")
                    call grid_init(grid,name="ANT-80KM",mtype="polar_stereographic",units="kilometers", &
                           lon180=.TRUE.,dx=80.d0,nx=79,dy=80.d0,ny=74,lambda=0.d0,phi=-71.d0)

                case("ANT-40KM")
                    call grid_init(grid,name="ANT-40KM",mtype="polar_stereographic",units="kilometers", &
                           lon180=.TRUE.,dx=40.d0,nx=157,dy=40.d0,ny=147,lambda=0.d0,phi=-71.d0)

                case("ANT-20KM")
                    call grid_init(grid,name="ANT-20KM",mtype="polar_stereographic",units="kilometers", &
                           lon180=.TRUE.,dx=20.d0,nx=313,dy=20.d0,ny=293,lambda=0.d0,phi=-71.d0)

                case("ANT-10KM")
                    call grid_init(grid,name="ANT-10KM",mtype="polar_stereographic",units="kilometers", &
                           lon180=.TRUE.,dx=10.d0,nx=625,dy=10.d0,ny=585,lambda=0.d0,phi=-71.d0)

                case("ANT-5KM")
                    call grid_init(grid,name="ANT-5KM",mtype="polar_stereographic",units="kilometers", &
                           lon180=.TRUE.,dx=5.d0,nx=1249,dy=5.d0,ny=1169,lambda=0.d0,phi=-71.d0)

                case("ANT-1KM")
                    call grid_init(grid,name="ANT-1KM",mtype="polar_stereographic",units="kilometers", &
                           lon180=.TRUE.,dx=1.d0,nx=6241,dy=1.d0,ny=5841,lambda=0.d0,phi=-71.d0)

                case DEFAULT
                    write(*,*) "rembo_grid_define:: error: grid name not recognized: "//trim(grid_name)
                    stop 

            end select

        end if 

        return 

    end subroutine rembo_grid_define
    
    subroutine rembo_cpu_time(time,time0,dtime)
        ! Calculate time intervals using system_clock.

        ! Note: for mulithreading, cpu_time() won't work properly.
        ! Instead, system_clock() should be used as it is here, 
        ! unless use_cpu_time=.TRUE. 

        !$ use omp_lib

        implicit none 

        real(8), intent(OUT) :: time 
        real(8), intent(IN),  optional :: time0 
        real(8), intent(OUT), optional :: dtime 

        ! Local variables
        logical    :: using_omp 
        integer(4) :: clock 
        integer(4) :: clock_rate
        integer(4) :: clock_max 
        real(8)    :: wtime 

        ! Check openmp status - set global variable to use as a switch 
        using_omp = .FALSE. 
        !$ using_omp = .TRUE.

        if (using_omp) then 
            ! --------------------------------------
            ! omp_get_wtime must be used for multithread openmp execution to get timing on master thread 
            ! The following lines will overwrite time with the result from omp_get_wtime on the master thread 

            !$ time = omp_get_wtime()

            ! --------------------------------------
            
        else 

            ! cpu_time can be used for serial execution to get timing on 1 processor
            call cpu_time(time)

        end if 

        if (present(dtime)) then 
            ! Calculate time interval 

            if (.not. present(time0)) then  
                write(*,*) "rembo_cpu_time:: Error: time0 argument is missing, but necessary."
                stop
            end if 
            
            ! Calculate the difference between current time and time0 in [s]
            dtime = time - time0

            ! Limit dtime to non-zero number 
            if (dtime .eq. 0.0d0) then
                write(*,*) "rembo_cpu_time:: Error: dtime cannot equal zero - check precision of timing variables, &
                            &which should be real(kind=8) to maintain precision."
                write(*,*) "clock", time, time0, dtime  
                stop  
            end if 

        end if 

        return 

    end subroutine rembo_cpu_time 

end module rembo_defs 
