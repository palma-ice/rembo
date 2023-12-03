module rembo1

    use precision, only : wp
    use rembo_defs, only : pi 

    implicit none


    private
    public :: rembo1_calc_atm


contains

    subroutine rembo1_calc_atm(tt,qq,pp,rain,snow,t0,rhum,S,ap,dh,rrco2,zs,dzs,mask,dx,dt,sec_day)

        implicit none

        real(wp), intent(INOUT) :: tt(:,:)
        real(wp), intent(INOUT) :: qq(:,:)
        real(wp), intent(INOUT) :: pp(:,:)
        real(wp), intent(IN) :: rain(:,:)
        real(wp), intent(IN) :: snow(:,:)
        real(wp), intent(IN) :: t0(:,:)
        real(wp), intent(IN) :: rhum(:,:)
        real(wp), intent(IN) :: S(:,:)
        real(wp), intent(IN) :: ap(:,:)
        real(wp), intent(IN) :: dh(:,:)
        real(wp), intent(IN) :: rrco2(:,:)
        real(wp), intent(IN) :: zs(:,:)
        real(wp), intent(IN) :: dzs(:,:)
        integer,  intent(IN) :: mask(:,:)
        real(wp), intent(IN) :: dx
        real(wp), intent(IN) :: dt
        real(wp), intent(IN) :: sec_day

        ! Local variables
        integer :: nx, ny 

        real(wp), allocatable :: kappa(:,:)

        integer, parameter :: solver = 2        ! 1: adi

        real(wp), parameter :: Lm       = 3.35e5
        real(wp), parameter :: Ls       = 2.84e6
        real(wp), parameter :: Lw       = 2.50e6
        real(wp), parameter :: rho_w    = 1000.0
        real(wp), parameter :: ta       = 224.58
        real(wp), parameter :: tb       = 1.94
        real(wp), parameter :: tkappa   = 1.0e12
        real(wp), parameter :: tce      = 8.37e6
        real(wp), parameter :: tlfac    = 0.0065
        real(wp), parameter :: trfac    = -1.0e3
        
        real(wp), parameter :: tau      = 5.0
        real(wp), parameter :: lapse    = 0.0
        real(wp), parameter :: h_e      = 2000.0
        real(wp), parameter :: scaleT   = 0.0
        real(wp), parameter :: p_k      = 50.0
        !real(wp), parameter :: p_k_eastfrac =  1.0
        real(wp), parameter :: pkappa   = 35e4
        real(wp), parameter :: pce      = 1.0
        real(wp), parameter :: prfac    = -1e-3
        real(wp), parameter :: ppfac    = -0.03

        nx = size(tt,1)
        ny = size(tt,2)
        
        allocate(kappa(nx,ny))

        ! Update energy balance
        kappa = tkappa
        call emb_temper(tt,t0,rain,snow,S,ap,dh,rrco2,zs,real(mask,wp),kappa, &
                        Lm,Ls,Lw,rho_w,ta,tb,tce,tlfac,trfac,dx,dt,solver)

        ! Update moisture balance
        !kappa = pkappa
        !call emb_precip(pp,qq,tt,rhum,zs,dzs,real(mask,wp),kappa,tau, &
        !                lapse,h_e,scaleT,p_k,pce,prfac,sec_day,dx,dt,solver)

        return

    end subroutine rembo1_calc_atm

    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Subroutine :  e m b _ t e m p e r
    ! Author     :  Alex Robinson
    ! Purpose    :  Determine the temperature over land through 
    !               diffusion, relaxing to data values over ocean
    ! dt [days]
    ! dxl [m]
    ! pp  [kg/m2-s]
    ! tt, t0  [°C]   - Sea-level T
    ! rhum [0:1]
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine emb_temper(tt,t0,rain,snow,S,ap,dh,rrco2,zs,mr,kappa, &
                            Lm,Ls,Lw,rho_w,ta,tb,tce,tlfac,trfac,dx,dt,solver)
    
        implicit none

        real(wp), intent(INOUT) :: tt(:,:)
        real(wp), intent(IN)  :: t0(:,:)
        real(wp), intent(IN)  :: rain(:,:)
        real(wp), intent(IN)  :: snow(:,:)
        real(wp), intent(IN)  :: S(:,:)
        real(wp), intent(IN)  :: ap(:,:)
        real(wp), intent(IN)  :: dh(:,:)
        real(wp), intent(IN)  :: rrco2(:,:)
        real(wp), intent(IN)  :: zs(:,:)
        real(wp), intent(IN)  :: mr(:,:)
        real(wp), intent(IN)  :: kappa(:,:)
        real(wp), intent(IN)  :: Lm
        real(wp), intent(IN)  :: Ls
        real(wp), intent(IN)  :: Lw
        real(wp), intent(IN)  :: rho_w
        real(wp), intent(IN)  :: ta
        real(wp), intent(IN)  :: tb
        real(wp), intent(IN)  :: tce
        real(wp), intent(IN)  :: tlfac
        real(wp), intent(IN)  :: trfac
        real(wp), intent(IN)  :: dx
        real(wp), intent(IN)  :: dt
        integer,  intent(IN)  :: solver

        ! Local variables
        integer  :: i, j, nx, ny 
        real(wp) :: gama, beta, dxl, dt_sec, chck
        real(wp), allocatable :: alpha(:,:)
        real(wp), allocatable :: F(:,:)
        real(wp), allocatable :: relax(:,:)
        
        nx = size(tt,1)
        ny = size(tt,2)

        allocate(alpha(nx,ny))
        allocate(F(nx,ny))
        allocate(relax(nx,ny))

        ! Convert dt from [sec] to [sec]
        dt_sec = dt 

        ! Constants, in time units of sec
        alpha = kappa * dt_sec / (dx*dx*tce)
        gama = dt_sec / tce 
        beta = tb

        ! Scale the relaxation mask
        relax = mr * trfac
        
        ! Get temperature forcing 
        F =  S*(1-ap) + rrco2 + (Lw*rain + Ls*snow) &
            - (ta - tb*tlfac*zs) - dh*rho_w*Lm

        if (solver .eq. 1) then
            call adi(tt,t0,F,alpha,relax,gama,beta)
        else      
            call explicit(tt,t0,F,alpha,relax,gama,beta)
        end if

        return
    
    end subroutine emb_temper

    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Subroutine :  e m b _ p r e c i p
    ! Author     :  Alex Robinson
    ! Purpose    :  Determine an amount of precipitation based
    !               on the water content in the atmosphere for current
    !               temperature
    ! dt [days]
    ! dxl [m]
    ! P  [kg/m2-s]
    ! T  [°C]   - 'sea-level' temperature, to be corrected by lapse
    ! rhum [0:1]
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine emb_precip(pp,qq,tt,rhum,zs,dzs,mr,kappa,tau, &
                            lapse,h_e,scaleT,p_k,pce,prfac,sec_day,dx,dt,solver)
    
        implicit none

        real(wp), intent(INOUT)  :: pp(:,:)
        real(wp), intent(INOUT)  :: qq(:,:)
        real(wp), intent(IN)  :: tt(:,:)
        real(wp), intent(IN)  :: rhum(:,:)
        real(wp), intent(IN)  :: zs(:,:)
        real(wp), intent(IN)  :: dzs(:,:)
        real(wp), intent(IN)  :: mr(:,:)
        real(wp), intent(IN)  :: kappa(:,:)
        real(wp), intent(IN)  :: tau
        real(wp), intent(IN)  :: lapse
        real(wp), intent(IN)  :: h_e            ! Characteristic elevation (usually 2000m)
        real(wp), intent(IN)  :: scaleT 
        real(wp), intent(IN)  :: p_k
        real(wp), intent(IN)  :: pce
        real(wp), intent(IN)  :: prfac
        real(wp), intent(IN)  :: sec_day
        real(wp), intent(IN)  :: dx 
        real(wp), intent(IN)  :: dt
        integer,  intent(IN)  :: solver

        ! Local variables
        integer :: i, j, nx, ny 
        real(wp) :: dt_sec, gama, chck
        real(wp) :: tau_sec

        real(wp), allocatable :: p0(:,:)
        real(wp), allocatable :: q0(:,:)
        real(wp), allocatable :: q_sat(:,:)
        real(wp), allocatable :: relax(:,:)
        real(wp), allocatable :: alpha(:,:)
        real(wp), allocatable :: tmp(:,:)
        real(wp), allocatable :: fr(:,:)

        nx = size(pp,1)
        ny = size(pp,2)

        allocate(p0(nx,ny))
        allocate(q0(nx,ny))
        allocate(q_sat(nx,ny))
        allocate(relax(nx,ny))
        allocate(alpha(nx,ny))
        allocate(tmp(nx,ny))
        allocate(fr(nx,ny))

        ! Convert time [sec] to [sec]
        dt_sec = dt 

        ! Constants, in time units of sec
        alpha = kappa*dt_sec/(dx*dx*pce)
        gama  = dt_sec / pce   
        
        tau_sec = tau * sec_day ! usually 5 days
        
        ! Scale the relaxation mask
        relax = mr * prfac
        
        ! scale the temperature to surface T, 
        ! and increase by scaling factor where no relaxation occurs
        tmp = tt - lapse*zs + scaleT*(1.d0-mr)       
         
        do j = 1, ny
        do i = 1, nx
            q_sat(j,i) = emb_FQSAT(tmp(i,j)) * h_e * airdens(zs(i,j))
        end do
        end do
            
        ! Get current Q based on relative humidity
        q0 = q_sat*rhum
    
        ! Adjust by qfactor (only increases precip outside Greenland)
        !q0 = q0 * qqfac
        ! ajr (2023-11-29): disabled for now, probably was not being used!! 

        ! Use diffusion to determine the new Q, forced by current P 
        if (solver .eq. 1) then     
            call adi(qq,q0,-pp,alpha,relax,gama,0.0_wp)
        else
            call explicit(qq,q0,-pp,alpha,relax,gama,0.0_wp)
        end if
       
        ! Hack in case q comes out negative (which shouldn't happen, but does in some cases!)
        ! ajr: find out why it goes negative, emb_PRECIP=1,emb_TEMPER=2
        ! ...or if that happens anymore?? ajr, 10.01.2009
        qq = max(qq,0.d0)
    
        ! Using the new value of Q, determine P
        ! Precip is either the total Q spread over eg tau=5days or
        ! the amount of Q over the saturation level
        !pp = max(qq/tau,(qq - q_sat)/dt_sec) 
        
        ! Only p1
        !pp = qq/tau
        
        ! Only p1b
        !n = 4.d0
        !fr = (qq/q_sat)**n
        !pp = (qq/tau)* fr
        
        ! Only p1c
        !fr = 1.d0 + p_k*dzs   !! now this is calculated outside in subroutine precip_grad...
        pp = (qq/tau_sec) * (1.d0 + p_k*dzs)
            
        return
    
    end subroutine emb_precip    
    
    subroutine kappa(kapfac,lats,lons,zs,mask,day,klat,klon,kamp,kzs,kdir,dx,days_in_year)
    
        implicit none

        real(wp), intent(OUT) :: kapfac(:,:)
        real(wp), intent(IN)  :: lats(:,:)
        real(wp), intent(IN)  :: lons(:,:)
        real(wp), intent(IN)  :: zs(:,:)
        real(wp), intent(IN)  :: mask(:,:)
        integer,  intent(IN)  :: day
        real(wp), intent(IN)  :: klat
        real(wp), intent(IN)  :: klon
        real(wp), intent(IN)  :: kamp
        real(wp), intent(IN)  :: kzs
        real(wp), intent(IN)  :: kdir
        real(wp), intent(IN)  :: dx
        real(wp), intent(IN)  :: days_in_year

        ! Local variables
        integer  :: nx, ny, i, j
        integer  :: f_lat,f_season,f_zs, f_lon
        real(wp) :: latmin, latmax, lonmin, lonmax
        real(wp) :: dudx, dudy, inv_2dx
        
        real(wp), allocatable :: uu_dir(:,:)
        
        nx = size(kapfac,2)
        ny = size(kapfac,1)
        
        allocate(uu_dir(ny,nx))
      
        inv_2dx = 1.0 / (2.0 * dx)
      
!       klat      = kappalat
!       klon      = kappalon
!       kamp      = kappaamp   ! dimensionless
!       kzs1      = kappazs    ! dimensionless
!       if (present(kzs)) kzs1 = kzs
      
        latmin = minval(lats); latmax = maxval(lats)
        lonmin = minval(lons); lonmax = maxval(lons)
        
        kapfac = 1.0

        ! Adjust kappa with latitude
        ! Make a scale that goes from 1 at low lat to the fraction desired at top of grid
        ! kappalat (-1:1), negative means lower latitudes have higher kappa,
        ! positive kappa lat means higher latitudes have higher kappa
        kapfac = 1.0 + klat*(2.0*(lats-latmin)/(latmax-latmin) - 1.0) 
      
        ! Adjust kappa with latitude
        ! Make a scale that goes from 1 at low lat to the fraction desired at top of grid
        ! kappalat (-1:1), negative means lower latitudes have higher kappa,
        ! positive kappa lat means higher latitudes have higher kappa
        kapfac = kapfac * ( 1.0 + klon*(2.0*(lons-lonmin)/(lonmax-lonmin) - 1.0) )
        
        
        ! Adjust kappa seasonally
        kapfac = kapfac * (1.0 + kamp*cos((day-15)*2.0*pi/days_in_year))

        ! Adjust kappa with elevation
        kapfac = kapfac* (1.0 + kzs*zs)
        
        ! Adjust kappa with direction of elevation
        ! Assign fraction as 0, assuming all points are East-facing
        uu_dir = 0.0
      
        do j = 2, ny-1
        do i = 2, nx-1
        
            dudx = (zs(i+1,j) - zs(i-1,j)) * inv_2dx
            dudy = (zs(i,j+1) - zs(i,j-1)) * inv_2dx

            ! Correct the fraction for western facing points
            if (dudx .gt. 0.0) uu_dir(i,j) = kdir

        end do
        end do
      
        where ( mask .ne. 2.0 ) kapfac = kapfac* (1.0 + uu_dir)
  
        return
    
    end subroutine kappa
    
      ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  R c o 2
  ! Author     :  Alex Robinson
  ! Purpose    :  Determine the radiative forcing for certain 
  !               co2 concentration
  !               (Based on paper by ? - taken from Andrey)
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
  function Rco2(co2)
    
    real(wp) :: Rco2, co2
    
    real(wp), parameter :: co2_0    = 280.d0
    real(wp), parameter :: Rco2_fac = 5.35d0 
    
    !if ( co2_0 .ne. 280.d0 ) write(*,*) "co2_0: ",co2_0
    !write(*,*) "co2: ",co2,", co2_0: ",co2_0
    !write(*,*) "dlog( co2/co2_0 ): ",dlog( co2 / co2_0 )
    Rco2 = Rco2_fac * log( co2 / co2_0 )
  
    return

  end function Rco2
  
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Subroutine :  a i r d e n s i t y
    ! Author     :  Alex Robinson
    ! Purpose    :  Determine the air density depending on elevation
    ! zs [m]
    ! rho [kg/m3]
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
    function airdens(zs)
    
        real(wp), intent(IN) :: zs
        real(wp) :: airdens
    
        airdens = 1.3d0 * exp(-zs/8600.d0)
  
        return
    end function airdens
        

    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Subroutine :  e m b _ F Q S A T
    ! Purpose    :  Same as FQSAT in climber, except real (8),
    !               returns sat. specific humidity, as a function 
    !               of °C
    ! Author     :  Alex Robinson (24. June 2008)
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    function emb_FQSAT(T)  

        implicit none
            
        real(wp), intent(IN) :: T
        real(wp) :: emb_FQSAT

        emb_FQSAT=3.8d-3*EXP(17.67d0*T/(T+243.5d0))

        return

    end function emb_FQSAT
    

    ! ======== PDE SOLVERS ==========
    
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Subroutine :  e x p l i c i t
    ! Purpose    :  Perform 1 timestep of explicit routine on 2d grid
    ! Author     :  Alex Robinson (24. June 2008)
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine explicit(uu,u0,F,alpha,relax,gama,beta)

        implicit none
        
        real(wp), intent(OUT) :: uu(:,:)
        real(wp), intent(IN)  :: u0(:,:)
        real(wp), intent(IN)  :: F(:,:)
        real(wp), intent(IN)  :: alpha(:,:)
        real(wp), intent(IN)  :: relax(:,:)
        real(wp), intent(IN)  :: gama
        real(wp), intent(IN)  :: beta
        
        ! Local variables
        integer :: i, j, k, nx, ny 
        real(wp) :: gam
        real(wp), allocatable :: u(:,:)
        real(wp), allocatable :: bet(:,:)
        real(wp), allocatable :: alph(:,:)
        real(wp), allocatable :: rel(:,:)
        
        nx = size(uu,1)
        ny = size(uu,2)

        allocate(u(nx,ny))
        allocate(bet(nx,ny))
        allocate(alph(nx,ny))
        allocate(rel(nx,ny))
        
        ! Constants
        alph = alpha
        gam  = gama
        bet  = 1.d0-4.d0*alph - gam*beta + gam*relax
    
        rel = gam*relax*u0
    
        ! Update boundary conditions
        uu(1,1:ny)  = u0(1,1:ny)
        uu(nx,1:ny) = u0(nx,1:ny)
        uu(1:nx,1)  = u0(1:nx,1)
        uu(1:nx,ny) = u0(1:nx,ny)
    
        ! Define matrix (to fill in boundary regions)
        u = uu
    
        ! explicit scheme
        do j= 2, ny-1
        do i= 2, nx-1
            u(i,j)= bet(i,j) * uu(i,j)   &
                    + alph(i,j) * (uu(i+1,j)+uu(i-1,j)+uu(i,j+1)+uu(i,j-1)) &
                    + gam*F(i,j) - rel(i,j)
        end do
        end do
      
        ! Send to output matrix
        uu = u
    
        return
    
    end subroutine explicit

    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Subroutine :  a d i
    ! Purpose    :  Perform 1 timestep of adi routine on 2d grid
    !               kappa varying in space (prescribed)
    ! Author     :  Alex Robinson (01. Aug 2008)
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine adi(uu,u0,F,alpha,relax,gama,beta)

        implicit none
    
        real(wp), intent(OUT) :: uu(:,:)
        real(wp), intent(IN)  :: u0(:,:)
        real(wp), intent(IN)  :: F(:,:)
        real(wp), intent(IN)  :: alpha(:,:)
        real(wp), intent(IN)  :: relax(:,:)
        real(wp), intent(IN)  :: gama
        real(wp), intent(IN)  :: beta
        
        ! Local variables
        integer  :: i, j, k, nx, ny 
        real(wp) :: gam, bet0, bet
    
        real(wp), allocatable :: alph(:,:)
        real(wp), allocatable :: coeffs(:,:)
        real(wp), allocatable :: rel(:,:)
        
        real(wp), allocatable :: ax(:), bx(:), cx(:), uprimex(:), cprimex(:), mx(:), rx(:)
        real(wp), allocatable :: ay(:), by(:), cy(:), uprimey(:), cprimey(:), my(:), ry(:)
        
        nx = size(uu,1)
        ny = size(uu,2)

        allocate(alph(nx,ny))
        allocate(coeffs(nx,ny))
        allocate(rel(nx,ny))

        allocate(ax(nx),bx(nx),cx(nx),uprimex(nx),cprimex(nx),mx(nx),rx(nx))
        allocate(ay(ny),by(ny),cy(ny),uprimey(ny),cprimey(ny),my(ny),ry(ny))
        
        ! Constants, divide by two because it is 1/2 timestep
        alph = alpha / 2.d0
        gam  = gama / 2.d0
        bet  = 0.d0             ! Explicit coefficient, 0 if using implicit
        bet0 = beta             ! Implicit coefficient, 0 if using explicit
        
        rel = gam*relax*u0     ! Assign 'boundary' values to relaxation matrix
        coeffs = 1.d0 - 2.d0*alph - gam*bet + gam*relax
    
        ! Define the three diagonal rows of the matrix: x
        ax = 0.d0
        bx = 1.d0
        cx = 0.d0

        ! Define the three diagonal rows of the matrix: y
        ay = 0.d0
        by = 1.d0
        cy = 0.d0
    
        ! Update boundary conditions
        uu(1:nx,1)  = u0(1:nx,1)
        uu(1:nx,ny) = u0(1:nx,ny)
        uu(1,1:ny)  = u0(1,1:ny)
        uu(nx,1:ny) = u0(nx,1:ny)
    
        ! First implicit x, explicit y - column by column
        do j = 2, ny-1
      
            cprimex(1) = cx(1) / bx(1)     
      
            ! RHS
            uprimex(1) = ( coeffs(1,j) * uu(1,j)     &
                        + alph(1,j) * (uu(1,j+1) + uu(1,j-1))    &
                        + gam*F(1,j) - rel(1,j) ) / bx(1)     
      
            ! Define the three diagonal rows of the matrix: x
            ax(2:nx-1) = -alph(2:nx-1,j)
            bx(2:nx-1) = 1 + 2.d0 * alph(2:nx-1,j) + gam*bet0
            cx(2:nx-1) = -alph(2:nx-1,j) 
            
            do i = 2, nx

                cprimex(i) = cx(i) / (bx(i) - ax(i) * cprimex(i-1))
        
                ! RHS
                uprimex(i) = ( ( coeffs(i,j) * uu(i,j)      &
                            + alph(i,j) * (uu(i,j+1) + uu(i,j-1))      & 
                            + gam*F(i,j) - rel(i,j) ) - ax(i)    &
                            * uprimex(i-1)) / (bx(i) - ax(i) * cprimex(i-1))        
        
            end do 
      
            do i = nx-1, 2, -1
                uu(i,j) = uprimex(i) - cprimex(i) * uu(i+1,j)
            end do
      
        end do ! End of j loop
    
        ! Now explicit x, implicit y, row by row
        do i = 2, nx-1
      
            cprimey(1) = cy(1) / by(1)
      
            ! RHS
            uprimey(1) = ( coeffs(i,1) * uu(i,1)      &
                            + alph(i,1) * (uu(i+1,1) + uu(i-1,1))   &  
                            + gam*F(i,1) - rel(i,1) ) / by(1)
      
            ! Define the three diagonal rows of the matrix: y
            ay(2:ny-1) = -alph(i,2:ny-1)
            by(2:ny-1) = 1 + 2.d0 * alph(i,2:ny-1) + gam*bet0
            cy(2:ny-1) = -alph(i,2:ny-1)
    
            do j = 2, ny

                cprimey(j) = cy(j) / ( by(j) - ay(j) * cprimey(j-1) )
        
                ! RHS               
                uprimey(j) = ( ( coeffs(i,j) * uu(i,j)      &
                            + alph(i,j) * (uu(i+1,j) + uu(i-1,j))     &  
                            + gam*F(i,j) - rel(i,j) ) - ay(j)    &
                            * uprimey(j-1)) / (by(j) - ay(j) * cprimey(j-1))                  

            end do
      
            do j = ny-1, 2, -1
                uu(j,i) = uprimey(j) - cprimey(j) * uu(i,j+1)
            end do
      
        end do  ! End of i-loop
    
        return
    
    end subroutine adi
      
end module rembo1