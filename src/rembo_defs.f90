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
    integer,  parameter :: prec = sp 

    ! Write flags 
    logical, parameter :: rembo_write_log = .TRUE. 

    ! Missing value and aliases
    real(prec), parameter :: MISSING_VALUE_DEFAULT = real(-9999.0,prec)
    real(prec), parameter :: MISSING_VALUE = MISSING_VALUE_DEFAULT
    real(prec), parameter :: MV = MISSING_VALUE_DEFAULT
    
    ! Error distance (very large), error index, and smallest number epsilon 
    real(prec), parameter :: ERR_DIST = real(1E8,prec) 
    integer,    parameter :: ERR_IND  = -1 
    real(prec), parameter :: eps      = real(1E-8,prec) 
    
    ! Mathematical constants
    real(prec), parameter :: pi  = real(2._dp*acos(0.0_dp),prec)
    real(prec), parameter :: degrees_to_radians = real(pi / 180._dp,prec)  ! Conversion factor between radians and degrees
    real(prec), parameter :: radians_to_degrees = real(180._dp / pi,prec)  ! Conversion factor between degrees and radians
    
    ! The constants below should be loaded using the global subroutine
    ! defined below `rembo_constants_load`.
    ! Note: The key limitation imposed by defining the parameters defined 
    ! globally is that these constants must be the same for all domains 
    ! being run in the same program. 

    ! Physical constants 
    real(prec) :: sec_year       ! [s] seconds per year 
    real(prec) :: g              ! Gravitational accel.  [m s-2]
    real(prec) :: T0             ! Reference freezing temperature [K] 
    real(prec) :: rho_ice        ! Density ice           [kg m-3] 
    real(prec) :: rho_w          ! Density water         [kg m-3] 
    real(prec) :: rho_sw         ! Density seawater      [kg m-3] 
    real(prec) :: rho_a          ! Density asthenosphere [kg m-3] 
    real(prec) :: rho_m          ! Density mantle (lith) [kg m-3]

public   ! All yelmo defs are public

contains 

    function rembo_get_precision() result(rembo_prec)

        implicit none 

        integer :: rembo_prec 

        rembo_prec = kind(prec)

        return 

    end function rembo_get_precision

        
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
        call nml_read(filename,"rembo_constants","T0",          T0,         init=init_pars)
        
        call nml_read(filename,"rembo_constants","rho_ice",     rho_ice,    init=init_pars)
        call nml_read(filename,"rembo_constants","rho_w",       rho_w,      init=init_pars)
        call nml_read(filename,"rembo_constants","rho_sw",      rho_sw,     init=init_pars)
        call nml_read(filename,"rembo_constants","rho_a",       rho_a,      init=init_pars)
        call nml_read(filename,"rembo_constants","rho_m",       rho_m,      init=init_pars)
        
        if (rembo_write_log) then 
            write(*,*) "yelmo:: loaded global constants:"
            write(*,*) "    sec_year = ", sec_year 
            write(*,*) "    g        = ", g 
            write(*,*) "    T0       = ", T0 
            write(*,*) "    rho_ice  = ", rho_ice 
            write(*,*) "    rho_w    = ", rho_w 
            write(*,*) "    rho_sw   = ", rho_sw 
            write(*,*) "    rho_a    = ", rho_a 
            write(*,*) "    rho_m    = ", rho_m 
            
        end if 

        return

    end subroutine rembo_global_init
    

end module rembo_defs 
