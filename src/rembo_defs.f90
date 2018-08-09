module rembo_defs 

    use nml 
    !use coord 

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
