module rembo_atm 

    use coord 
    use solvers 

    use rembo_defs 
    use rembo_physics 

    implicit none 

    private
    public rembo_calc_atmosphere

contains 

    subroutine rembo_calc_atmosphere(now,emb,bnd,par,day,year)
        ! Calculate the equilibrium rembo atmosphere solution 
        ! for a given day of a given year (can be representative of a month)

        implicit none 

        type(rembo_state_class),    intent(INOUT) :: now     ! rembo atmosphere variables
        type(diffusion_class),      intent(INOUT) :: emb     ! Diffusion variables
        type(rembo_boundary_class), intent(IN)    :: bnd     ! rembo boundary variables
        type(rembo_param_class),    intent(IN)    :: par     ! rembo parameters
        integer, intent(IN) :: day 
        integer, intent(IN) :: year 

        ! Local variables
        real(wp) :: als_max, als_min, afac, tmid    ! ajr: to do: move to parameters!!

        ! Get the surface pressure 
        now%sp = calc_sp(bnd%z_srf)

        ! Get geostrophic components and magnitude of velocity
        now%ug  = calc_u_geo(now%dZdy,bnd%f)
        now%vg  = calc_v_geo(now%dZdx,bnd%f)
        now%uvg = calc_magnitude(now%ug,now%vg)

        ! Get horizontal gradient of temperature
        !call gradient_1D(now%dtsldx,now%dtsldy,now%dtsldxy,now%t2m+bnd%z_srf*par%gamma,par%nx,par%ny,par%dx)
        
        ! Calculate katabatic wind field
        now%u_k  = calc_u_kata(now%dtsldx,bnd%dzsdx,f_k=par%f_k)
        now%v_k  = calc_v_kata(now%dtsldy,bnd%dzsdy,f_k=par%f_k)
        now%uv_k = calc_magnitude(now%u_k,now%v_k)

!         ! Combine wind fields 
!         now%ug = now%ug + now%u_k 
!         now%vg = now%vg + now%v_k 
!         now%uvg = calc_magnitude(now%ug,now%vg)

        ! Calculate vertical velocity
        now%ww  = calc_w(now%ug,now%vg,bnd%dzsdx,bnd%dzsdy)

        ! Calculate the surface velocity 
!         now%u_s = calc_u_surf()
!         now%v_s = calc_v_surf()
        now%u_s = now%ug 
        now%v_s = now%vg 
        now%uv_s = calc_magnitude(now%u_s,now%v_s)

        ! Calculate the surface albedo 
        als_max =   0.80
        als_min =   0.69
        afac     =  -0.18
        tmid     = 275.35
        now%al_s = calc_albedo_t2m(now%t2m,als_min,als_max,afac,tmid)

        ! ajr: TO do: implement loop updating t2m and al_s !!!

        ! Calculate the planetary albedo 
        now%al_p = par%alp_a + par%alp_b*now%al_s - par%alp_c*now%tcw 
        where(now%al_p .lt. 0.0) now%al_p = 0.0
        where(now%al_p .gt. 1.0) now%al_p = 1.0  
        
        ! Calculate the outgoing long-wave radiation at toa
        now%lwu = par%lwu_a + par%lwu_b*(now%t2m-273.15) + par%lwu_c*now%S

        ! Calculate the incoming short-wave radiation at toa
        now%swd = (1.0-now%al_p)*now%S 

        ! Calculate radiative forcing of CO2
        now%rco2_a  = calc_rad_co2(now%co2_a)

        ! Calculate cloud fraction
        now%cc   = calc_cloudfrac(now%t2m,now%ccw,now%rho_a, &
                                         par%nk1,par%nk2,par%nk3)
        
        ! Calculate latent heat flux at the surface

        ! ajr: TO DO 

        ! Calculate sensible heat flux at the surface

        ! ajr: TO DO 


        ! Calculate energy balance on low-resolution grid
        call rembo_calc_en(now%t2m,emb,par,day,bnd%z_srf,now%t2m_bnd, &
                          swn=(now%swd - now%swd_s*(1.0-now%al_s)), &
                          lwn=(- now%lwu + now%lwu_s - now%lwd_s), &
                          shf=now%shf_s,lhf=now%lhf_s, &
                          lhp=(par%Lw*(now%pr-now%sf) + par%Ls*now%sf)*1.0, &
                          rco2=par%en_kdT + now%rco2_a) 

        ! Calculate inversion correction for moisture balance
        now%ct2m = calc_ttcorr1(now%t2m,bnd%z_srf,-2.4e-3,-3.7e-1,106.0)

        ! Get saturated specific humidity and sat. total water content
        now%q_sat   = calc_qsat(now%t2m+now%ct2m,bnd%z_srf) 
        now%tcw_sat = now%q_sat * now%rho_a * par%H_e

        ! Calculate the current total water content
        now%q_s = calc_qs(now%t2m+now%ct2m,bnd%z_srf,par%e0,par%c1) 
        now%tcw = now%q_s * now%rho_a * par%H_e
        now%q_r = now%tcw / now%tcw_sat 

!         ! Now calculate the condensation rate (kg m**-2 s**-1)
!         now%c_w = calc_condensation(now%tcw,now%q_r,bnd%dzsdxy,now%ww, &
!                                     par%k_c,par%k_x)

        ! Calculate the current cloud water content
        call rembo_calc_ccw(now%pr,now%ccw,now%c_w,emb,par,ww=now%ww)   


        ! Calculate snowfall (kg m**-2 s**-1)
        now%sf = calc_snowfrac(now%t2m,par%sf_a,par%sf_b) * now%pr 

        ! ## Calculate surface fluxes ###

        ! Shortwave radiation down (surface)
        now%swd_s = calc_radshort_surf_down(now%S,now%cc,bnd%z_srf,&
                                par%swds_a,par%swds_b,par%swds_c) 

        ! Longwave radiation down (surface)
        now%lwd_s = calc_rad_surf_down(now%t2m,now%cc,now%tcw, &
                                       par%lwds_a,par%lwds_b,par%lwds_c,par%lwds_d)

        return 

    end subroutine rembo_calc_atmosphere

    subroutine rembo_calc_en(t2m,emb,par,day,z_srf,t2m_bnd,swn,lwn,shf,lhf,lhp,rco2)

        implicit none 

        real(wp),                intent(INOUT) :: t2m(:,:)
        type(diffusion_class),   intent(INOUT) :: emb  
        real(wp),                intent(IN)    :: z_srf(:,:)
        real(wp),                intent(IN)    :: t2m_bnd(:,:)
        real(wp),                intent(IN)    :: swn(:,:)
        real(wp),                intent(IN)    :: lwn(:,:)
        real(wp),                intent(IN)    :: shf(:,:)
        real(wp),                intent(IN)    :: lhf(:,:)
        real(wp),                intent(IN)    :: lhp(:,:)
        real(wp),                intent(IN)    :: rco2(:,:)
        type(rembo_param_class), intent(IN)    :: par 
        integer,                 intent(IN)    :: day 
        
        ! Local variables 
        real(wp), allocatable :: gamma(:,:)  
        real(wp) :: tsl_fac 
        integer :: q 

        ! Allocate working arrays
        allocate(gamma(size(t2m,1),size(t2m,2))) 

        ! Get the current lapse rate (gamma=winter, gamma2=summer)  <= CHECK !!
        gamma = par%gamma  &
           + (0.5*cos((day-15)*2.0*pi/day_year)-0.5)*(par%gamma-par%gamma2)

        ! Get the tsl => column energy conversion
        ! tsl_fac = H_a[m] c_p[J kg-1 K-1] rho_a[kg m-3] = [J m-2 K-1]
        tsl_fac = par%en_Ha *1000.0 *1.225 !* 1.225 ! = 8.6e6
        
        ! Get the 2D energy diffusion coefficient
        ! emb%kappa = par%en_D 
        emb%kappa = par%en_D_win + (par%en_D_sum-par%en_D_win)*(0.5-0.5*cos((day-15)*2.0*pi/day_year))
        
        ! Boundary sea-level temperature, tsl
        call map_field(emb%map_toemb,"tsl_bnd",t2m_bnd+gamma*z_srf,emb%tsl_bnd,method="radius",fill=.TRUE.,missing_value=dble(missing_value))
          
        ! Sea level temperature, tsl
        call map_field(emb%map_toemb,"tsl",t2m+gamma*z_srf,emb%tsl,method="radius",fill=.TRUE.,missing_value=dble(missing_value))
            
        ! Radiative forcing, en_F
        call map_field(emb%map_toemb,"en_F",swn + lwn + (shf+lhf) + lhp + rco2, &
                       emb%en_F,method="radius",fill=.TRUE.,missing_value=dble(missing_value))

        ! Radiation
        emb%en     = emb%tsl     *tsl_fac
        emb%en_bnd = emb%tsl_bnd *tsl_fac

        ! Calculate radiative balance over the day
        do q = 1, par%en_nstep
            where (emb%mask .eq. 1) emb%en = emb%en_bnd
            call adv_diff_2D(emb%en,emb%en_bnd,emb%en_F,relax=emb%mask, &
                             dx=real(emb%grid%G%dx*emb%grid%xy_conv,wp), &
                             dy=real(emb%grid%G%dx*emb%grid%xy_conv,wp), &
                             dt=par%en_dt,kappa=emb%kappa,k_relax=par%en_kr) !, &
!                              v_x=emb%ug,v_y=emb%vg)
!             call solve_diff_2D_adi(emb%en,emb%en_bnd,emb%en_F,relax=emb%mask, &
!                                    dx=emb%grid%G%dx*emb%grid%xy_conv, &
!                                    dy=emb%grid%G%dx*emb%grid%xy_conv, &
!                                    dt=par%en_dt,kappa=emb%kappa,k_relax=par%en_kr)

        end do 

        ! Re-calculate temperature 
        emb%tsl = emb%en /tsl_fac

        call map_field(emb%map_fromemb,"tsl",emb%tsl,t2m,method="nn",fill=.TRUE.,missing_value=dble(missing_value))

        ! ajr: TO DO !!
!         t2m =  interp_bilinear(is_points=.TRUE.,x=emb%grid%G%x,y=emb%grid%G%y,z=emb%tsl, &
!                                xout=emb%map_fromemb%x,yout=emb%map_fromemb%y, &
!                                missing_value=missing_value)

        t2m = t2m - gamma*z_srf

        return 

    end subroutine rembo_calc_en 

    subroutine rembo_calc_ccw(pr,ccw,c_w,emb,par,ww)

        implicit none 

        real(wp),                intent(OUT)   :: pr(:,:)
        real(wp),                intent(OUT)   :: ccw(:,:)
        real(wp),                intent(OUT)   :: c_w(:,:)
        type(diffusion_class),   intent(INOUT) :: emb 
        type(rembo_param_class), intent(IN)    :: par 
        real(wp),                intent(IN)    :: ww(:,:)
        
        ! Local variables
        real(wp)  :: dtsl_mean, lat0, lat1
        integer   :: q  
        !real(wp)  :: dtsl_mean    ! Average temp anomaly away from boundary temps 

        !emb%tcw     = calc_qs(emb%tsl_bnd+dtsl_mean-par%gamma*emb%z_srf,emb%z_srf,par%e0,par%c1) *emb%rho_a *par%H_e
        emb%tcw_bnd = calc_qs(emb%tsl_bnd-par%gamma*emb%z_srf,emb%z_srf,par%e0,par%c1) *emb%rho_a *par%H_e
        emb%tcw     = emb%tcw_bnd

        emb%qr = emb%tcw / (calc_qsat(emb%tsl_bnd+dtsl_mean-par%gamma*emb%z_srf,emb%z_srf) *emb%rho_a *par%H_e)

        ! Get 2D diffusion constant
        emb%kappaw = par%ccw_D

        ! Note: removed any lat/elevation dependence of diffusion coefficient (ajr, 2015-08-31)
!         lat0 = minval(emb%grid%lat)
!         lat1 = maxval(emb%grid%lat)
!         emb%kappaw = emb%kappaw* & 
!             (1.d0 + par%en_kl*(2.d0*(emb%grid%lat-lat0)/(lat1-lat0)-1.d0))

!         emb%kappaw = emb%kappaw* (1.d0 + par%en_kz*emb%z_srf)

        ! Get vertical wind for precipitation 
        ! ajr: TO DO 
!         call map_field(emb%map_toemb,"ww",ww,emb%ww,method="nng",sigma=240.0,fill=.TRUE.,missing_value=missing_value)
!         call map_field(emb%map_toemb,"ww",ww,emb%ccw_pr,method="radius",sigma=240.d0,fill=.TRUE.,missing_value=missing_value)
!         emb%ww = (emb%ww+emb%ccw_pr)/2.d0 

!         ! Get current moisture content [kg m-2]
!         call map_field(emb%map_toemb,"ccw",ccw,emb%ccw,method="radius",fill=.TRUE.,missing_value=missing_value)
            
        ! Calculate the initial lo-res precipitation rate [kg m-2 s-1]
        !call map_field(emb%map_toemb,"pr", pr, emb%ccw_pr,method="nn",fill=.TRUE.,missing_value=missing_value)
        ! ajr: TO DO : calculate a basic rate somehow avoiding history 


        ! Now calculate the condensation rate [kg m-2 s-1]
!         call map_field(emb%map_toemb,"ww",ww,emb%ww,method="nng",sigma=240.d0,fill=.TRUE.,missing_value=missing_value)
        emb%ccw_cw = calc_condensation(emb%tcw,emb%qr,emb%ww,par%k_c,par%k_x)

        ! Get moisture balance forcing [kg m-2 s-1]
        emb%ccw_F = emb%ccw_cw - emb%ccw_pr 
        where(emb%mask==1) emb%ccw_F = 0.d0 

        ! ajr: Why was this called twice?
!         call map_field(emb%map_toemb,"ww",ww,emb%ww,method="radius",sigma=240.0,fill=.TRUE.,missing_value=missing_value)
        
        ! Calculate moisture balance to equilibrium
        do q = 1, par%ccw_nstep
            where (emb%mask .eq. 1) emb%ccw = emb%ccw_bnd
            call adv_diff_2D(emb%ccw,emb%ccw_bnd*emb%tcw/emb%tcw_bnd, &
                             emb%ccw_F,relax=emb%mask, &
                             dx=real(emb%grid%G%dx*emb%grid%xy_conv,wp), &
                             dy=real(emb%grid%G%dx*emb%grid%xy_conv,wp), &
                             dt=par%ccw_dt,kappa=emb%kappaw,k_relax=par%ccw_kr) !, &
!                              v_x=emb%ug,v_y=emb%vg)

!             call solve_diff_2D_adi(emb%ccw,emb%ccw_bnd*emb%tcw/emb%tcw_bnd, &
!                                    emb%ccw_F,relax=emb%mask, &
!                                    dx=emb%grid%G%dx*emb%grid%xy_conv, &
!                                    dy=emb%grid%G%dx*emb%grid%xy_conv, &
!                                    dt=par%ccw_dt,kappa=emb%kappaw,k_relax=par%ccw_kr)
            where (emb%ccw .lt. 0.d0) emb%ccw = 0.d0

            ! Now calculate the precipitation rate on emb grid (kg m**-2 s**-1)
            ! *Cheaper than reinterpolating, but necessary to update moisture balance
            emb%ccw_pr = calc_precip(emb%ccw,emb%ww,emb%tsl-par%gamma*emb%z_srf,emb%z_srf,par%k_w,par%k_z,par%k_t) 

            ! ajr: Add equil. condition!!!

        end do 

        ! Send cloud moisture content back to main domain pts
!         call map_field(emb%map_fromemb,"ccw",emb%ccw,ccw,method="quadrant",fill=.TRUE.,missing_value=missing_value)
        ! ajr: TO DO 
!         ccw =  interp_bilinear(is_points=.TRUE.,x=emb%grid%G%x,y=emb%grid%G%y,z=emb%ccw, &
!                                xout=emb%map_fromemb%x,yout=emb%map_fromemb%y, &
!                                missing_value=missing_value)

        ! ajr: TO DO 
!         c_w = interp_bilinear(is_points=.TRUE.,x=emb%grid%G%x,y=emb%grid%G%y,z=emb%ccw_cw, &
!                                xout=emb%map_fromemb%x,yout=emb%map_fromemb%y, &
!                                missing_value=missing_value)

        ! Now calculate the high resolution precipitation rate (kg m**-2 s**-1)
!         now%pr = calc_precip(now%ccw,now%ww,now%t2m,par%k_w,par%k_z,par%k_t) 
        ! ajr: TO DO 
!         pr = interp_bilinear(is_points=.TRUE.,x=emb%grid%G%x,y=emb%grid%G%y,z=emb%ccw_pr, &
!                                xout=emb%map_fromemb%x,yout=emb%map_fromemb%y, &
!                                missing_value=missing_value)

        return 

    end subroutine rembo_calc_ccw 

end module rembo_atm 
