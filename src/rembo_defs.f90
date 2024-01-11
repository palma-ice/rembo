module rembo_defs 

    use precision
    use nml
    use monthlydaily, only : monthlydaily_class

    use coordinates_mapping_scrip, only : map_scrip_class

    implicit none 

    ! =========================================================================
    !
    ! CONSTANTS (program precision, global constants)
    !
    ! =========================================================================

    ! Write flags 
    logical, parameter :: rembo_write_log = .TRUE. 

    ! Missing value and aliases
    real(wp), parameter :: MISSING_VALUE_DEFAULT = real(-9999.0,wp)
    real(wp), parameter :: MISSING_VALUE = MISSING_VALUE_DEFAULT
    real(wp), parameter :: MV = MISSING_VALUE_DEFAULT
    
    ! Error distance (very large), error index, and smallest number epsilon 
    real(wp), parameter :: ERR_DIST = real(1E8,wp) 
    integer,  parameter :: ERR_IND  = -1 
    real(wp), parameter :: eps      = real(1E-8,wp) 
    
    ! Mathematical constants
    real(wp), parameter :: pi  = real(2._dp*acos(0.0_dp),wp)
    real(wp), parameter :: degrees_to_radians = real(pi / 180._dp,wp)  ! Conversion factor between radians and degrees
    real(wp), parameter :: radians_to_degrees = real(180._dp / pi,wp)  ! Conversion factor between degrees and radians
    
    logical :: rembo_use_omp 

    ! The constants below should be loaded using the global subroutine
    ! defined below `rembo_constants_load`.
    ! Note: The key limitation imposed by defining the parameters defined 
    ! globally is that these constants must be the same for all domains 
    ! being run in the same program. 

    type rembo_phys_const_class

        ! Physical constants
        integer  :: nm              ! [--] Number of months in a year
        integer  :: ndm             ! [--] Number of days in a month
        real(wp) :: sec_year        ! [s] seconds per year
        real(wp) :: sec_day_ref     ! [s] seconds per day, normally 
        real(wp) :: g               ! Gravitational accel.  [m s-2]
        real(wp) :: omega           ! Coriolis constant [omega = 7.2921d-5]
        real(wp) :: T0              ! Reference freezing temperature [K] 
        real(wp) :: rho_ice         ! Density ice           [kg m-3] 
        real(wp) :: rho_w           ! Density water         [kg m-3] 

        
        ! Internal parameters
        integer  :: nd
        real(wp) :: day_year
        real(wp) :: day_month
        real(wp) :: month_year
        real(wp) :: sec_day 
        real(wp) :: sec_frac 

    end type

    ! First define all parameters needed to represent a given domain
    type rembo_param_class

        type(rembo_phys_const_class) :: c           ! Physical constants
        
        character(len=256)  :: domain
        character(len=256)  :: grid_name
        character(len=256)  :: grid_name_emb
        character(len=256)  :: grid_name_hi
        character(len=256)  :: grid_path
        character(len=256)  :: grid_path_emb
        character(len=256)  :: grid_path_hi
        character(len=256)  :: restart 
        
        integer             :: npts, nx, ny
        real(wp)            :: dx 
        
        logical             :: rembo1  
        character(len=56)   :: solver 
        character(len=56)   :: step 
        real(wp)            :: mask_zs_min
        real(wp)            :: mask_radius

        ! Physics
        real(wp)    :: H_e   ! Precip vapor scale height (m)
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

    type rgrid_class 

        ! Grid name 
        character(len=256) :: name 
        
        ! Parameters
        integer    :: nx, ny, npts
        real(wp)   :: dx, dy

        ! Projection parameters (optional)
        character(len=256) :: mtype 
        real(wp)   :: lambda
        real(wp)   :: phi
        real(wp)   :: alpha
        real(wp)   :: scale
        real(wp)   :: x_e
        real(wp)   :: y_n
        real(wp)   :: semi_major_axis
        real(wp)   :: inverse_flattening
        logical    :: is_sphere 
        logical    :: is_projection 

        ! Axes
        real(wp), allocatable :: xc(:)    
        real(wp), allocatable :: yc(:) 

        ! Grid arrays 
        real(wp), allocatable :: x(:,:)
        real(wp), allocatable :: y(:,:)
        real(wp), allocatable :: lon(:,:)
        real(wp), allocatable :: lat(:,:)
        real(wp), allocatable :: area(:,:)
        
    end type

    type rembo_boundary_class

        ! Annual boundary variables
        real(wp), allocatable :: z_srf(:,:)    ! [m]     Surface elevation
        real(wp), allocatable :: f_ice(:,:)    ! [--]    Fraction of land-based ice coverage in cell
        real(wp), allocatable :: f_shlf(:,:)   ! [--]    Fraction of floating (shelf) ice coverage in cell
        
        ! Derived annual boundary variables
        real(wp), allocatable :: f(:,:)        ! [--]    Coriolis parameter
        real(wp), allocatable :: dzsdx(:,:)
        real(wp), allocatable :: dzsdy(:,:)
        real(wp), allocatable :: dzsdxy(:,:)

        integer,  allocatable :: mask(:,:)     ! [--]    1: solve model, -1: fix values to boundary
        
    end type 

    ! Now define all variables of the domain
    type rembo_state_class

        ! Monthly forcing variables 
        real(wp), allocatable :: S(:,:)        ! [W m-2] Insolation top-of-atmosphere
        real(wp), allocatable :: t2m_bnd(:,:)  ! [K]     Near-surface temperature (used for boundary)
        real(wp), allocatable :: tsl_bnd(:,:)  ! [K]     Near-surface temperature at sea level (used for boundary)
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
        real(wp), allocatable :: tce(:,:)
        real(wp), allocatable :: tcm(:,:)
        real(wp), allocatable :: gamma(:,:)
        real(wp), allocatable :: tsl(:,:)
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
        real(wp), allocatable :: ccw_bnd(:,:)
        real(wp), allocatable :: ug(:,:), vg(:,:), uvg(:,:), ww(:,:), cc(:,:)
        real(wp), allocatable :: swd(:,:), lwu(:,:), al_p(:,:), at(:,:)
        real(wp), allocatable :: swd_s(:,:), lwd_s(:,:), shf_s(:,:), lhf_s(:,:), lwu_s(:,:)
        real(wp), allocatable :: u_s(:,:), v_s(:,:), uv_s(:,:)  
        
        real(wp), allocatable :: u_k(:,:), v_k(:,:), uv_k(:,:), dtsldx(:,:), dtsldy(:,:), dtsldxy(:,:)
    end type 

    ! Define all variables needed for diffusion on lo-res grid
    type diffusion_class

        type(rgrid_class)     :: grid ! EMB diffusion resolution grid
        type(map_scrip_class) :: map_toemb, map_fromemb  ! map EMB <=> rembo grid

        ! Relaxation mask, topography  
        integer,  allocatable :: mask(:,:)                  ! 1: solve model, -1: fix values to boundary
        real(wp), allocatable :: z_srf(:,:), rho_a(:,:) 
        real(wp), allocatable :: dzsdx(:,:), dzsdy(:,:), dzsdxy(:,:) 

        ! Energy and moisture balance variables
        real(wp), allocatable :: tsl(:,:), tsl_bnd(:,:), tsl_F(:,:)  
        real(wp), allocatable :: ccw(:,:), ccw_bnd(:,:), ccw_F(:,:) 
        real(wp), allocatable :: ccw_cw(:,:), ccw_pr(:,:) 
        real(wp), allocatable :: tcw(:,:)
        real(wp), allocatable :: ug(:,:), vg(:,:), uvg(:,:), ww(:,:), q_r(:,:)  

        ! Diffusion 
        real(wp), allocatable :: kappa_t(:,:), kappa_w(:,:) 
        real(wp) :: en_dt, tsl_fac, en_kr, en_kz
        integer :: en_nstep
        logical :: bnd_pr 
    end type 

    type rembo_class

        type(rembo_param_class)      :: par         ! Model parameters
        type(rgrid_class)            :: grid        ! Grid definition   (from coordinates module)
        type(rgrid_class)            :: gridhi      ! High-resolution grid definition   (from coordinates module)
        type(monthlydaily_class)     :: mm          ! Monthly to daily interpolation weights

        ! Boundary variables
        type(rembo_boundary_class) :: bnd 

        ! current variables, month and various averages
        type(rembo_state_class) :: now, mon(12), ann

        ! Variables and grid definitions for energy-moisture balance calculations
        type(diffusion_class) :: emb

    end type

    type rembo_forcing_class
        ! Climatology and forcing data for a whole year 
        ! This is just a convenient data holding class, which can be useful for
        ! holding boundary forcing data. But it is not necessary
        ! to use it in a program, as the fields are passed individually to rembo_update.

        real(wp), allocatable :: z_srf(:,:)      ! [m]     Surface elevation
        real(wp), allocatable :: t2m(:,:,:)      ! [K]     Near-surface temperature (used for boundary)
        real(wp), allocatable :: tsl(:,:,:)      ! [K]     Sea-level temperature (used for boundary)
        real(wp), allocatable :: al_s(:,:,:)     ! [--]    Surface albedo 
        real(wp), allocatable :: co2_a           ! [ppm]   Atmospheric CO2 concentration
        real(wp), allocatable :: Z(:,:,:)        ! [m?]    Geopotential height of 750 Mb layer

        real(wp), allocatable :: tcwv(:,:,:)     ! [kg m-2] Total column water vapor (used for boundary)
        
        ! Extra variables that can be helpful
        real(wp), allocatable :: pr(:,:,:)       ! [km m-2 s-1] Precipitation rate

    end type 

    !public   ! All rembo defs are public
    
    ! safer ==>
    private 

    public :: dp, sp, wp, rembo_write_log, MISSING_VALUE_DEFAULT, MISSING_VALUE, MV
    public :: ERR_DIST, ERR_IND, eps, pi, degrees_to_radians, radians_to_degrees
    public :: rembo_use_omp 
    public :: rembo_phys_const_class
    public :: rembo_param_class, rembo_boundary_class
    public :: rgrid_class
    public :: rembo_state_class 
    public :: diffusion_class 
    public :: rembo_class 
    public :: rembo_forcing_class
    
    public :: rembo_get_working_precision
    public :: rembo_parse_path
    public :: rembo_global_init
    public :: rembo_cpu_time

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

    subroutine rembo_global_init()
        ! This routine will ensure that the rembo global variable rembo_use_omp is initialized.
        ! The routine should be called internally via rembo_init(), which means it will
        ! be called multiple times if multiple domains are intialized. But this shouldn't
        ! be a problem. However, rembo_use_omp will have the same value for all domains. 

        !$ use omp_lib 

        ! Local variables
        integer :: n_threads 
        character(len=10) :: n_threads_str 

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

        return

    end subroutine rembo_global_init
    
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
