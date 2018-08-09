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

        integer             :: npts, nx, ny  
        character(len=256)  :: boundary(30)
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
        real(wp)    :: k_p_now, k_t, k_w, ps_a, ps_b 
        real(wp)    :: f_k 
        real(wp)    :: nk1, nk2, nk3
        real(wp)    :: teff_sigma 
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

    ! Now define all variables of the domain
    type rembo_state_class

        ! Static (per year or greater) variables
        integer,  allocatable :: mask(:,:), basin(:,:) 
        real(wp), allocatable :: zs(:,:), zb(:,:), H(:,:) 
        real(wp), allocatable :: co2_a(:,:), rco2_a(:,:), rho_a(:,:), f(:,:), sp(:,:)
        real(wp), allocatable :: dzsdx(:,:), dzsdy(:,:), dzsdxy(:,:)

        ! Seasonal (daily) variables that require averaging
        real(wp), allocatable :: t2m(:,:), teff(:,:), ct2m(:,:), pp(:,:), ps(:,:)  
        real(wp), allocatable :: q_s(:,:), q_sat(:,:), q_r(:,:)
        real(wp), allocatable :: tcw(:,:), tcw_sat(:,:)
        real(wp), allocatable :: ccw(:,:), c_w(:,:), ccw_prev(:,:) 
        real(wp), allocatable :: ug(:,:), vg(:,:), uvg(:,:), ww(:,:), cc(:,:)
        real(wp), allocatable :: S(:,:), dS(:,:), swd(:,:), lwu(:,:), al_s(:,:), alp(:,:), at(:,:)
        real(wp), allocatable :: tsurf(:,:), swd_s(:,:), lwd_s(:,:), shf_s(:,:), lhf_s(:,:), lwu_s(:,:)
        real(wp), allocatable :: u_s(:,:), v_s(:,:), uv_s(:,:)  
        real(wp), allocatable :: Z(:,:), dZdx(:,:), dZdy(:,:)
        real(wp), allocatable :: u_k(:,:), v_k(:,:), uv_k(:,:), dtsldx(:,:), dtsldy(:,:), dtsldxy(:,:)
    end type 

    type insol_class

        integer               :: year, nd 
        real(wp)              :: year_ref 
        character(len=256)    :: fldr
        logical               :: calendar_time  
        real(wp), allocatable :: lat(:,:) 
        real(wp), allocatable :: S(:,:,:), S0(:,:,:), dS(:,:,:)

    end type 

    ! Define all variables needed for diffusion on lo-res grid
    type diffusion_class

        type(grid_class)     :: grid ! EMB diffusion resolution grid
        type(map_class)      :: map_toemb, map_fromemb  ! map EMB => rembo grid

        ! Relaxation mask, topography  
        integer,  allocatable :: mask(:,:)
        real(wp), allocatable :: zs(:,:), rho_a(:,:) 
        real(wp), allocatable :: dzsdx(:,:), dzsdy(:,:), dzsdxy(:,:) 

        ! Energy and moisture balance variables
        real(wp), allocatable :: tsl(:,:), tsl_bnd(:,:), tsl_bnd0(:,:) 
        real(wp), allocatable :: en(:,:), en_bnd(:,:), en_F(:,:) 
        real(wp), allocatable :: ccw(:,:), ccw_bnd(:,:), ccw_F(:,:) 
        real(wp), allocatable :: ccw_cw(:,:), ccw_pp(:,:) 
        real(wp), allocatable :: tcw(:,:), tcw_bnd(:,:)
        real(wp), allocatable :: ug(:,:), vg(:,:), uvg(:,:), ww(:,:), qr(:,:)  

        ! Diffusion 
        real(wp), allocatable :: kappa(:,:), kappaw(:,:) 
        real(wp) :: en_dt, tsl_fac, en_kr, en_kz
        integer :: en_nstep
        logical :: bnd_pp 
    end type 

    type boundary_opt_class 
        logical :: t2m, ct2m, pp, q_r, tcw, ccw, ug, vg, Z, cc
        logical :: swd_s, lwd_s, u_s, v_s
        logical :: swd, lwu, alp 
    end type

    type rembo_class

        type(rembo_param_class)   :: par        ! physical parameters
        type(boundary_opt_class)  :: bnd, bnd0  ! boundary switches (bnd0 for equilibration)

        ! Daily variables, month and annual averages, forcing variables
        type(rembo_state_class) :: now, mon, ann 

        ! Variables and grid definitions for energy-moisture balance calculations
        type(diffusion_class) :: emb 

        ! Variables holding the daily insolation (present-day and anomalies)
        type(insol_class) :: ins 

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
