module emb 
    ! This module holds all physics functions related to the 
    ! 2D energy-moisture balance (emb) atmospheric model of REMBO 

    use rembo_defs 

    implicit none 
    
    private 
    public :: calc_condensation
    public :: calc_precip
    public :: calc_snowfrac
    public :: calc_magnitude
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
    
contains
    
    elemental function calc_condensation(tcw,qr,ww,k_c,k_x) result(c_w)

        implicit none 

        real(wp), intent(IN) :: tcw, qr, ww
        real(wp), intent(IN) :: k_c, k_x
        real(wp) :: c_w 

        c_w = (tcw/(k_c*sec_day))*qr * (1.d0 + k_x*ww)
        c_w = max(c_w,0.d0)

        return

    end function calc_condensation

    elemental function calc_precip(ccw,ww,tt,zs,k_w,k_z,k_t) result(pp)

        implicit none 

        real(wp), intent(IN) :: ccw, ww, tt, zs
        real(wp), intent(IN) :: k_w, k_z, k_t  
        real(wp) :: pp 

        real(wp) :: k_tmp 

!         pp = (ccw/k_w) * (1.d0 - k_t*tt) * 1.d0/rho_a
        
        k_tmp = k_t 
!         if (zs > 0.d0) k_tmp = 0.d0 

        pp = (ccw/k_w) * (1.d0 + k_z*ww - k_tmp*tt)

        pp = max(pp,0.d0)

        return

    end function calc_precip

    elemental function calc_snowfrac(t2m,a,b) result(f)
        ! Return the fraction of snow from total precipitation
        ! expected for a given temperature
        
        implicit none 

        real(wp), intent(IN) :: t2m, a, b 
        real(wp)             :: f 

        f = -0.5d0*tanh(a*(t2m-b))+0.5d0 
!         f = 1.d0 - 1.d0 / (1.d0+exp(-a*2.d0*(t2m-b)))
        
        return 

    end function calc_snowfrac 

    ! Get the vector magnitude from two components
    elemental function calc_magnitude(u,v) result(umag)
        implicit none 
        real(wp), intent(IN)  :: u, v 
        real(wp) :: umag 

        umag = sqrt(u*u+v*v)

        return
    end function calc_magnitude 

    ! Get coriolis parameter (1/s)
    elemental function calc_coriolis(lat) result(f)
        implicit none 
        real(wp), intent(IN)  :: lat
        real(wp) :: f

        f = 2.d0*omega*sin(lat*degrees_to_radians)

        return
    end function calc_coriolis

! # Convert geopotential into geopotential height, (m2/s2) => (m)
    elemental function calc_geo_height(phi) result(Z)
        implicit none 
        real(wp), intent(IN)  :: phi
        real(wp) :: Z

        Z = phi/g

        return
    end function calc_geo_height

    ! Get horizontal geostrophic wind component, u
    elemental function calc_u_geo(dZdy,f) result(ug)
        implicit none 
        real(wp), intent(IN)  :: dZdy, f
        real(wp) :: ug

        ug = -(g / f) * dZdy

        return
    end function calc_u_geo

    ! Get horizontal geostrophic wind component, v
    elemental function calc_v_geo(dZdx,f) result(vg)
        implicit none 
        real(wp), intent(IN)  :: dZdx, f
        real(wp) :: vg

        vg = (g / f) * dZdx

        return
    end function calc_v_geo

    ! Get vertical wind
    elemental function calc_w(u,v,dzdx,dzdy) result(w)
        implicit none 
        real(wp), intent(IN)  :: u, v, dzdx, dzdy
        real(wp) :: w

        w = u*dzdx + v*dzdy 

        return
    end function calc_w

    ! Get vertical wind from omega (Pa/s)
    elemental function calc_w_omega(omega,T,p) result(w)
        implicit none 
        real(wp), intent(IN)  :: omega,T,p
        real(wp) :: w
        real(wp) :: rho
        real(wp), parameter :: rgas = 287.058d0
        
        rho  = p/(rgas*T)         ! density => kg/m3
        w    = -1*omega/(rho*g) 
        return
    end function calc_w_omega

    ! Get katabatic wind component, u 
    elemental function calc_u_kata(dTdx,dzsdx,f_k) result(uk)
        implicit none 
        real(wp), intent(IN)  :: dTdx,dzsdx, f_k
        real(wp) :: uk
        
        uk = f_k* min(0.d0,dTdx*dzsdx)
        uk = sign(uk,-dzsdx)

        return
    end function calc_u_kata

    ! Get katabatic wind component, v 
    elemental function calc_v_kata(dTdy,dzsdy,f_k) result(vk)
        implicit none 
        real(wp), intent(IN)  :: dTdy,dzsdy, f_k
        real(wp) :: vk
        
        vk = f_k* min(0.d0,dTdy*dzsdy)
        vk = sign(vk,-dzsdy)

        return
    end function calc_v_kata

    ! Calculate a correction factor to adjust 2-m temperature
    ! for inversion effects
    elemental function calc_ttcorr1(T2m,zs,a,b,c) result(ttcorr)
        implicit none 

        real(wp), intent(IN) :: T2m, zs, a, b, c
        real(wp) :: ttcorr 

        ! Get the difference with the actual 2m temperature
        ttcorr = a*zs + b*T2m + c

        return
    end function calc_ttcorr1 

    ! Surface pressure ( Pa=> kg/(m s2) )
    elemental function calc_sp(zs) result(sp)
        implicit none 
        real(wp), intent(IN)  :: zs
        real(wp) :: sp
        real(wp), parameter :: p0 = 101325d0
        real(wp), parameter :: M  = 0.0289644d0
        real(wp), parameter :: R  = 8.31447
        real(wp), parameter :: T0 = 298.15d0 

        sp = p0*exp(-g*M*zs/(R*T0))
        return
    end function calc_sp

    
    elemental function calc_esat(Ts) result(esat)
        implicit none 
        real(wp), intent(IN)  :: Ts
        real(wp) :: esat
        real(wp), parameter :: e0 = 6.112d0
        real(wp), parameter :: c1 = 17.67d0
        real(wp), parameter :: c2 = 243.5d0
        real(wp), parameter :: T0 = 273.15d0 

        esat = e0*exp((c1*(Ts-T0))/(c2+(Ts-T0)))*100d0

        return
    end function calc_esat

    ! Air density (kg/m3) for given elevation
    elemental function calc_airdens(zs) result(rho_a)
        implicit none 
        real(wp), intent(IN)  :: zs
        real(wp) :: rho_a 

        rho_a = 1.3d0 * exp(-zs/8.6d3)

        return
    end function calc_airdens

    ! Elevation corresponding to a given air pressure (m)
    elemental function calc_elevation(p) result(zs)
        implicit none 
        real(wp), intent(IN) :: p 
        real(wp) :: zs 
        real(wp), parameter :: p0 = 1013.25d0
        real(wp), parameter ::  g = 9.80665d0
        real(wp), parameter ::  M = 0.0289644d0
        real(wp), parameter ::  R = 8.31447d0
        real(wp), parameter :: T0 = 288.15d0 

        zs = -log(p/p0)*R*T0/(g*M)

        return
    end function calc_elevation

    ! Specific water content parameterization(kg/kg)
    ! Ts [K], zs [m] 
    elemental function calc_qs(Ts,zs,e0,c1) result(qs)
        implicit none 
        real(wp), intent(IN)  :: Ts, zs
        real(wp), intent(IN) :: e0, c1 
        real(wp) :: qs, esat, p
        real(wp), parameter :: ebs = 0.62198
        real(wp), parameter :: c2  = 243.5d0
        real(wp), parameter :: T0  = 273.15d0 

        ! First, calculate the saturation vapor pressure
        ! ( hPa *100 => Pa => kg/(m s2) )
        esat = e0*exp((c1*(Ts-T0))/(c2+(Ts-T0)))*100d0

        ! and the surface pressure
        p    = calc_sp(zs)

        ! Then calculate the specific humidity
        qs = ebs * esat / (p-(1.d0-ebs)*esat)

        return
    end function calc_qs

    ! Saturation specific water content (kg/kg)
    ! Ts [K], zs [m] 
    elemental function calc_qsat(Ts,zs) result(qsat)
        implicit none 
        real(wp), intent(IN)  :: Ts, zs
        real(wp) :: qsat 
        real(wp), parameter :: e0  = 6.112d0
        real(wp), parameter :: c1  = 17.67d0
 
        ! Calculate the specific humidity with 
        ! the correct saturation parameters
        qsat = calc_qs(Ts,zs,e0,c1)

        return
    end function calc_qsat

    ! Calculate radiative forcing of CO2 from the CO2 concentration
    elemental function calc_rad_co2(CO2) result(RCO2)

    real(wp), intent(IN) :: CO2
    real(wp)             :: RCO2

    real(wp), parameter :: CO2_0    = 280.d0
    real(wp), parameter :: RCO2_fac = 5.35d0 

    RCO2 = RCO2_fac * log( CO2 / CO2_0 )

    return

    end function calc_rad_co2
    
    ! Cloud cover fraction (1)
    ! tt [K], ccw [kg m2], rho_a [Pa]
    elemental function calc_cloudfrac(tt,ccw,rho_a,k1,k2,k3) result(cc)
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

end module emb 
