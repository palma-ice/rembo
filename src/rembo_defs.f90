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
        real(wp), allocatable :: tsl(:,:), tsl_bnd(:,:) 
        real(wp), allocatable :: en(:,:), en_bnd(:,:), en_F(:,:) 
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

        character(len=*), intent(IN)  :: filename
        
        ! Local variables
        logical :: init_pars 

        init_pars = .TRUE. 
        
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
    

end module rembo_defs 
