module rembo_physics 
    ! This module holds all physics functions related to the 
    ! 2D energy-moisture balance (emb) atmospheric model of REMBO 

    use rembo_defs 

    implicit none 

    private 
    public :: calc_total_column_energy
    public :: calc_albedo_t2m
    public :: calc_condensation
    public :: calc_precip
    public :: calc_snowfrac
    public :: calc_coriolis
    public :: calc_geo_height
    public :: calc_u_geo
    public :: calc_v_geo
    public :: calc_w
    public :: calc_w_omega
    public :: calc_u_kata
    public :: calc_v_kata
    public :: calc_ttcorr1
    public :: calc_sp
    public :: calc_esat
    public :: calc_airdens
    public :: calc_elevation
    public :: calc_qs
    public :: calc_qsat
    public :: calc_rad_co2
    public :: calc_cloudfrac
    public :: calc_rad_surf_up
    public :: calc_rad_surf_down
    public :: calc_radshort_surf_down
    public :: calc_ebs_sky
    public :: calc_sensible_heat_flux
    
    public :: d_dx 
    public :: d_dy 
    public :: calc_magnitude
    public :: calc_magnitude_from_staggered
    public :: calc_gradient_to_sealevel

contains
    
    elemental subroutine calc_total_column_energy(tce,t2m,z_srf,gamma,cv,rho_0,H_toa)
        ! Calculate the total column energy from near-surface air temperature
        ! following Willeit et al. (2022), Eq. A40 and subsequent text.
        ! tce = integral_{z_srf}^{H_toa} cv * rho * T dz

        implicit none

        real(wp), intent(OUT) :: tce            ! [J/m^2] Total column energy
        real(wp), intent(IN)  :: t2m            ! [K] Near-surface air temperature
        real(wp), intent(IN)  :: z_srf          ! [m] Surface elevation
        real(wp), intent(IN)  :: gamma          ! [K/m] Atmospheric lapse rate
        real(wp), intent(IN)  :: cv             ! [J / kg-K] Air specific heat capacity at constant volume
        real(wp), intent(IN)  :: rho_0          ! [kg/m^3] Air density at sea level
        real(wp), intent(IN)  :: H_toa          ! [m] Maximum considered height of atmosphere
        ! Local variables
        integer :: k, nz
        real(wp) :: rho_srf

        nz = 10 

        ! Get air density at surface elevation
        rho_srf = calc_airdens(z_srf)


        return

    end subroutine calc_total_column_energy

    elemental function calc_albedo_t2m(t2m,als_min,als_max,afac,tmid) result(al_s)
        ! Calculate the surface albedo of snow as a function of near-surface temperature
        ! Alexander Robinson, inspired from Slater et al, etc. 

        implicit none

        real(wp), intent(IN) :: t2m 
        real(wp), intent(IN) :: als_min, als_max, afac, tmid 
        real(wp) :: al_s 

        al_s = als_min + (als_max - als_min)*(0.5*tanh(afac*(t2m-tmid))+0.5)
        !al_s = als_min+(als_max-als_min)*afac*(t2m-tmid)

        return 

    end function calc_albedo_t2m

    elemental function calc_condensation(tcw,qr,ww,k_c,k_x,sec_day) result(c_w)

        implicit none 

        real(wp), intent(IN) :: tcw, qr, ww
        real(wp), intent(IN) :: k_c, k_x
        real(wp), intent(IN) :: sec_day
        real(wp) :: c_w 

        c_w = (tcw/(k_c*sec_day))*qr * (1.0 + k_x*ww)
        c_w = max(c_w,0.0)

        return

    end function calc_condensation

    elemental function calc_precip(ccw,ww,tt,zs,k_w,k_z,k_t) result(pp)

        implicit none 

        real(wp), intent(IN) :: ccw, ww, tt, zs
        real(wp), intent(IN) :: k_w, k_z, k_t  
        real(wp) :: pp 

        real(wp) :: k_tmp 

!         pp = (ccw/k_w) * (1.0 - k_t*tt) * 1.0/rho_a
        
        k_tmp = k_t 
!         if (zs > 0.0) k_tmp = 0.0 

        pp = (ccw/k_w) * (1.0 + k_z*ww - k_tmp*tt)

        pp = max(pp,0.0)

        return

    end function calc_precip

    elemental function calc_snowfrac(t2m,a,b) result(f)
        ! Return the fraction of snow from total precipitation
        ! expected for a given temperature
        
        implicit none 

        real(wp), intent(IN) :: t2m, a, b 
        real(wp)             :: f 

        f = -0.50*tanh(a*(t2m-b))+0.50 
!         f = 1.0 - 1.0 / (1.0+exp(-a*2.0*(t2m-b)))
        
        return 

    end function calc_snowfrac 

    ! Get coriolis parameter (1/s)
    elemental function calc_coriolis(lat,omega) result(f)

        implicit none

        real(wp), intent(IN) :: lat
        real(wp), intent(IN) :: omega
        real(wp) :: f

        f = 2.0*omega*sin(lat*degrees_to_radians)

        return

    end function calc_coriolis

    elemental function calc_geo_height(phi,g) result(Z)
        ! # Convert geopotential into geopotential height, (m2/s2) => (m)

        implicit none 
        real(wp), intent(IN)  :: phi
        real(wp), intent(IN)  :: g
        real(wp) :: Z

        Z = phi/g

        return

    end function calc_geo_height

    elemental function calc_u_geo(dZdy,f,g) result(ug)
        ! Get horizontal geostrophic wind component, u

        implicit none 

        real(wp), intent(IN)  :: dZdy
        real(wp), intent(IN)  :: f
        real(wp), intent(IN)  :: g
        real(wp) :: ug

        ug = -(g / f) * dZdy

        return

    end function calc_u_geo

    elemental function calc_v_geo(dZdx,f,g) result(vg)
        ! Get horizontal geostrophic wind component, v
        
        implicit none 
        
        real(wp), intent(IN)  :: dZdx
        real(wp), intent(IN)  :: f
        real(wp), intent(IN)  :: g
        real(wp) :: vg

        vg = (g / f) * dZdx

        return

    end function calc_v_geo

    subroutine calc_w(w,u,v,dzdx,dzdy)
        ! Get vertical wind
        
        implicit none 
        
        real(wp), intent(OUT) :: w(:,:)     ! aa-nodes
        real(wp), intent(IN)  :: u(:,:)     ! ac-nodes
        real(wp), intent(IN)  :: v(:,:)     ! ac-nodes 
        real(wp), intent(IN)  :: dzdx(:,:)  ! ac-nodes
        real(wp), intent(IN)  :: dzdy(:,:)  ! ac-nodes

        ! Local variables 
        integer :: i, j, nx, ny 
        integer :: im1, jm1 
        real(wp) :: u_aa, v_aa, dzdx_aa, dzdy_aa  

        nx = size(w,1)
        ny = size(w,2) 


        do j = 1, ny 
        do i = 1, nx 

            im1 = max(i-1,1)
            jm1 = max(j-1,1) 

            u_aa = 0.5*(u(i,j)+u(im1,j)) 
            v_aa = 0.5*(v(i,j)+v(im1,j)) 
            dzdx_aa = 0.5*(dzdx(i,j)+dzdx(im1,j))
            dzdy_aa = 0.5*(dzdy(i,j)+dzdy(i,jm1))
            
            w(i,j) = u_aa*dzdx_aa + v_aa*dzdy_aa

        end do 
        end do  

        return
    end subroutine calc_w

    elemental function calc_w_omega(omega,T,p,g) result(w)
        ! Get vertical wind from omega (Pa/s)
        
        implicit none 
        
        real(wp), intent(IN) :: omega
        real(wp), intent(IN) :: T
        real(wp), intent(IN) :: p
        real(wp), intent(IN) :: g
        real(wp) :: w
        real(wp) :: rho
        real(wp), parameter :: rgas = 287.058_wp
        
        rho  = p/(rgas*T)         ! density => kg/m3
        w    = -1.0*omega/(rho*g) 
        
        return

    end function calc_w_omega

    elemental function calc_u_kata(dTdx,dzsdx,f_k) result(uk)
        ! Get katabatic wind component, u 
        
        implicit none 
        
        real(wp), intent(IN) :: dTdx
        real(wp), intent(IN) :: dzsdx
        real(wp), intent(IN) :: f_k
        real(wp) :: uk
        
        uk = f_k* min(0.0,dTdx*dzsdx)
        uk = sign(uk,-dzsdx)

        return

    end function calc_u_kata

    elemental function calc_v_kata(dTdy,dzsdy,f_k) result(vk)
        ! Get katabatic wind component, v 
        
        implicit none 
        
        real(wp), intent(IN) :: dTdy
        real(wp), intent(IN) :: dzsdy
        real(wp), intent(IN) :: f_k
        real(wp) :: vk
        
        vk = f_k* min(0.0,dTdy*dzsdy)
        vk = sign(vk,-dzsdy)

        return

    end function calc_v_kata

    elemental function calc_ttcorr1(T2m,zs,a,b,c) result(ttcorr)
        ! Calculate a correction factor to adjust 2-m temperature
        ! for inversion effects
        
        implicit none 

        real(wp), intent(IN) :: T2m, zs, a, b, c
        real(wp) :: ttcorr 

        ! Get the difference with the actual 2m temperature
        ttcorr = a*zs + b*T2m + c

        return

    end function calc_ttcorr1 

    elemental function calc_sp(zs,g) result(sp)
        ! Surface pressure ( Pa=> kg/(m s2) )
        
        implicit none 

        real(wp), intent(IN) :: zs
        real(wp), intent(IN) :: g
        real(wp) :: sp

        real(wp), parameter :: p0 = 101325_wp
        real(wp), parameter :: M  = 0.0289644_wp
        real(wp), parameter :: R  = 8.31447_wp
        real(wp), parameter :: T0 = 298.15_wp 

        sp = p0*exp(-g*M*zs/(R*T0))

        return

    end function calc_sp

    elemental function calc_esat(Ts,T0) result(esat)

        implicit none 
        
        real(wp), intent(IN) :: Ts
        real(wp), intent(IN) :: T0
        real(wp) :: esat
        real(wp), parameter :: e0 = 6.112_wp
        real(wp), parameter :: c1 = 17.67_wp
        real(wp), parameter :: c2 = 243.5_wp

        esat = e0*exp((c1*(Ts-T0))/(c2+(Ts-T0)))*100_wp

        return

    end function calc_esat

    elemental function calc_airdens(zs) result(rho_a)
        ! Air density (kg/m3) for given elevation
    
        implicit none 

        real(wp), intent(IN) :: zs
        real(wp) :: rho_a 

        rho_a = 1.3_wp * exp(-zs/8.6d3)

        return

    end function calc_airdens

    elemental function calc_elevation(p,g) result(zs)
        ! Elevation corresponding to a given air pressure (m)
    
        implicit none 

        real(wp), intent(IN) :: p 
        real(wp), intent(IN) :: g
        real(wp) :: zs

        real(wp), parameter :: p0 = 1013.25_wp
        real(wp), parameter ::  M = 0.0289644_wp
        real(wp), parameter ::  R = 8.31447_wp
        
        !real(wp), parameter ::  g = 9.80665_wp
        real(wp), parameter :: T0 = 288.15_wp 

        zs = -log(p/p0)*R*T0/(g*M)

        return

    end function calc_elevation

    elemental function calc_qs(Ts,zs,e0,c1,g,T0) result(qs)
        ! Specific water content parameterization(kg/kg)
        ! Ts [K], zs [m] 
        
        implicit none

        real(wp), intent(IN) :: Ts
        real(wp), intent(IN) :: zs
        real(wp), intent(IN) :: e0
        real(wp), intent(IN) :: c1 
        real(wp), intent(IN) :: g
        real(wp), intent(IN) :: T0
        real(wp) :: qs
        
        real(wp) :: esat, p
        real(wp), parameter :: ebs = 0.62198_wp
        real(wp), parameter :: c2  = 243.5_wp
        
        ! First, calculate the saturation vapor pressure
        ! ( hPa *100 => Pa => kg/(m s2) )
        esat = e0*exp((c1*(Ts-T0))/(c2+(Ts-T0)))*100_wp

        ! and the surface pressure
        p    = calc_sp(zs,g)

        ! Then calculate the specific humidity
        qs = ebs * esat / (p-(1.0-ebs)*esat)

        return

    end function calc_qs

    elemental function calc_qsat(Ts,zs,g,T0) result(qsat)
        ! Saturation specific water content (kg/kg)
        ! Ts [K], zs [m] 
        
        implicit none 
        
        real(wp), intent(IN) :: Ts
        real(wp), intent(IN) :: zs
        real(wp), intent(IN) :: g
        real(wp), intent(IN) :: T0
        real(wp) :: qsat
        
        real(wp), parameter :: e0  = 6.112_wp
        real(wp), parameter :: c1  = 17.67_wp
 
        ! Calculate the specific humidity with 
        ! the correct saturation parameters
        qsat = calc_qs(Ts,zs,e0,c1,g,T0)

        return

    end function calc_qsat

    elemental function calc_rad_co2(CO2) result(RCO2)
        ! Calculate radiative forcing of CO2 from the CO2 concentration
    
        real(wp), intent(IN) :: CO2
        real(wp)             :: RCO2

        real(wp), parameter :: CO2_0    = 280.d0
        real(wp), parameter :: RCO2_fac = 5.35d0 

        RCO2 = RCO2_fac * log( CO2 / CO2_0 )

        return

    end function calc_rad_co2
    
    elemental function calc_cloudfrac(tt,ccw,rho_a,k1,k2,k3) result(cc)
        ! Cloud cover fraction (1)
        ! tt [K], ccw [kg m2], rho_a [Pa]
        
        implicit none 
        
        real(wp), intent(IN)  :: tt, ccw, rho_a, k1, k2, k3 
        real(wp) :: cc
 
        ! Calculate the specific humidity with 
        ! the correct saturation parameters
        cc = (k1+k2*ccw+k3*tt)  / rho_a

        ! Make sure to limit to appropriate range (probably unnecessary)
        cc = min(cc,1.d0)
        cc = max(cc,0.d0)

        return

    end function calc_cloudfrac

    
    elemental function calc_rad_surf_up(tsurf) result(lwu)
        ! Stefan-Boltzmann's law for outgoing long-wave radiation flux (W m-2)

        implicit none 

        real(wp), intent(IN) :: tsurf
        real(wp), parameter :: sigm = 5.67d-8 ! W m−2 K−4
        real(wp) :: lwu 

        lwu = sigm*(tsurf*tsurf*tsurf*tsurf)
        
        return 

    end function calc_rad_surf_up 

    elemental function calc_rad_surf_down(t2m,cc,tcw,a,b,c,d) result(lwd)
        ! Approximation for downward long-wave radiation at the surface (W m-2)

        implicit none 

        real(wp), intent(IN) :: t2m, tcw, cc, a, b, c, d 
        real(wp) :: lwd

        lwd = a + b*t2m + c*cc + d*tcw

        return 

    end function calc_rad_surf_down

    elemental function calc_radshort_surf_down(S,cc,zs,a,b,c) result(swd)
        ! Approximation for downward short-wave radiation at the surface (W m-2)
        ! Modified form of Eq. 3 of Konzelmann et al (1994)

        implicit none 

        real(wp), intent(IN) :: S, cc, zs, a, b, c 
        real(wp) :: swd 

        swd = S* (a - b*cc*exp(-c*zs))

        return

    end function calc_radshort_surf_down

    elemental function calc_ebs_sky(uv,alpha,beta,m) result(ebs)
        ! For Trigo et al lwd_s parameterization (not used)

        implicit none 

        real(wp), intent(IN) :: uv, alpha, beta, m 
        real(wp) :: ebs 

        ebs = 1.d0 - (1.d0+0.1d0*uv)*exp(-(alpha+beta*0.1d0*uv)**m)

        return 

    end function calc_ebs_sky 

    elemental function calc_sensible_heat_flux(ts,ta,wind,rhoa,csh_pos,csh_neg,cap) result(shf)
        ! Bulk formulation of sensible heat flux to the atmosphere
        ! Following Krapp et al. (2017)

        real(wp), intent(in) :: ts      !< surface temperature
        real(wp), intent(in) :: ta      !< air temperature
        real(wp), intent(in) :: wind    !< wind speed
        real(wp), intent(in) :: rhoa    !< air density
        real(wp), intent(in) :: csh_pos !< sensible heat exchange coefficient (for positive shf, ie summer)
        real(wp), intent(in) :: csh_neg !< sensible heat exchange coefficient (for positive shf)
        real(wp), intent(in) :: cap     !< air specific heat capacity
        real(wp) :: shf                 !< sensible heat flux [W/m2]
        
        ! Local variables
        real(wp) :: coeff

        ! enhancement over land
        if ( wind*(ts-ta) > 0.0_dp ) then
            coeff = csh_pos
        else
            coeff = csh_neg
        end if

        shf = coeff*cap*rhoa*wind*(ts-ta)

        return

    end function calc_sensible_heat_flux

    subroutine d_dx(du,u0,dx)  
        ! Horizontal gradient, x-component
    
        implicit none

        integer :: i, j, nx, ny
        real(wp), intent(IN) :: u0(:,:)
        real(wp) :: du(:,:)
        real(wp) :: dx, inv_dx, inv_2dx

        inv_dx  = 1.d0 / (dx)
        inv_2dx = 1.d0 / (2.d0 * dx)
        nx = size(u0,1)
        ny = size(u0,2)

        ! Assign boundary values
        du = 0.d0

if (.FALSE.) then 
        ! Centered, second order
        do i = 2, nx-1
            du(i,:) = (u0(i+1,:) - u0(i-1,:)) * inv_2dx
        end do

else 
        ! Staggered, second order 
        do i = 1, nx-1
            du(i,:) = (u0(i+1,:) - u0(i,:)) * inv_dx
        end do
end if 

        return

    end subroutine d_dx

    subroutine d_dy(du,u0,dx)  
        ! Horizontal gradient, y-component
    
        implicit none

        integer :: i, j, nx, ny
        real(wp), intent(IN) :: u0(:,:)
        real(wp) :: du(:,:)
        real(wp) :: dx, inv_dx, inv_2dx

        inv_dx  = 1.d0 / (dx)
        inv_2dx = 1.d0 / (2.d0 * dx)
        nx = size(u0,1)
        ny = size(u0,2)

        ! Assign boundary values
        du = 0.d0

if (.FALSE.) then 
        ! Centered, second order
        do j = 2, ny-1
            du(:,j) = (u0(:,j+1) - u0(:,j-1)) * inv_2dx
        end do

else 
        ! Staggered, second order 
        do j = 1, ny-1
            du(:,j) = (u0(:,j+1) - u0(:,j)) * inv_dx
        end do
end if 

        return

    end subroutine d_dy

    elemental function calc_magnitude(u,v) result(umag)
        ! Get the vector magnitude from two components
        
        implicit none 
        
        real(wp), intent(IN)  :: u, v 
        real(wp) :: umag 

        umag = sqrt(u*u+v*v)

        return

    end function calc_magnitude 

    function calc_magnitude_from_staggered(u,v) result(umag)
        ! Calculate the centered (aa-nodes) magnitude of a vector 
        ! from the staggered (ac-nodes) components

        implicit none 
        
        real(wp), intent(IN)  :: u(:,:), v(:,:)  
        real(wp) :: umag(size(u,1),size(u,2)) 

        ! Local variables 
        integer :: i, j, nx, ny 
        integer :: ip1, jp1, im1, jm1 
        real(wp) :: unow, vnow 

        nx = size(u,1)
        ny = size(u,2) 

        umag = 0.0_wp 

        do j = 1, ny 
        do i = 1, nx

            im1 = max(i-1,1)
            jm1 = max(j-1,1)

            unow = 0.5_wp*(u(i,j)+u(im1,j))
            vnow = 0.5_wp*(v(i,j)+v(i,jm1))
            umag(i,j) = sqrt(unow*unow+vnow*vnow)

        end do 
        end do 

        return

    end function calc_magnitude_from_staggered 
    
    subroutine calc_gradient_to_sealevel(dzsdxy,z_srf,z_sl,xx,yy)

        implicit none 

        real(wp), intent(OUT) :: dzsdxy(:,:)    ! aa-nodes 
        real(wp), intent(IN)  :: z_srf(:,:) 
        real(wp), intent(IN)  :: z_sl(:,:) 
        real(wp), intent(IN)  :: xx(:,:) 
        real(wp), intent(IN)  :: yy(:,:) 
        
        ! Local variables
        integer :: i, j, nx, ny
        real(wp), allocatable :: dists(:,:) 
        real(wp) :: dist_min 

        nx = size(z_srf,1)
        ny = size(z_srf,2)

        allocate(dists(nx,ny))

        dists = 0.0 

        do j = 1, ny 
        do i = 1, nx 

            ! Calculate distances to all sea-level points
            where (z_srf .eq. z_sl) dists = sqrt((xx-xx(i,j))**2+(yy-yy(i,j))**2)

            ! Determine minimum distance and calculate slope 
            dist_min = minval(dists,mask=dists .gt. 0.0)

            dzsdxy(i,j) = (z_srf(i,j)-z_sl(i,j)) / dist_min 

        end do 
        end do 


        return 

    end subroutine calc_gradient_to_sealevel

end module rembo_physics 
