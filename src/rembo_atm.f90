module rembo_atm 

    use coord 
    use solvers 

    use rembo_defs 
    use rembo_physics 

    implicit none 

    real(wp), parameter :: cap  = 1000.0_wp  !< specific heat capacitiy of air [J/(kg K)]
    ! ajr: move to rembo_defs later...


    private
    public :: rembo_calc_atmosphere

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
        real(wp) :: beta                            ! ajr: to do: move to parameters!!
        integer  :: iter, n_iter 
        integer  :: iter_surf, n_iter_surf 

        ! Get the current lapse rate (gamma=winter, gamma2=summer)  <= CHECK !!
        now%gamma = par%gamma  &
           + (0.5*cos((day-15)*2.0*pi/day_year)-0.5)*(par%gamma-par%gamma2)
        !write(*,*) "gamma: ", day, minval(gamma), maxval(gamma)
        
        ! Get the surface pressure 
        now%sp = calc_sp(bnd%z_srf)

        ! Get geostrophic components and magnitude of velocity
        now%ug  = calc_u_geo(now%dZdy,bnd%f)
        now%vg  = calc_v_geo(now%dZdx,bnd%f)
        now%uvg = calc_magnitude(now%ug,now%vg)

        ! Initialize t2m solution with boundary solution 
        now%t2m = now%t2m_bnd 
        
        ! Initialize emb object for this day 
        call rembo_calc_iterinit(emb,bnd%mask,bnd%z_srf,now%t2m_bnd,now%gamma,par,day)

        n_iter      = 2
        n_iter_surf = 2 

        do iter = 1, n_iter 

            ! Get horizontal gradient of temperature
            !call gradient_1D(now%dtsldx,now%dtsldy,now%dtsldxy,now%t2m+bnd%z_srf*par%gamma,par%nx,par%ny,par%dx)
        
            ! Calculate katabatic wind field (ac-nodes)
            now%u_k  = calc_u_kata(now%dtsldx,bnd%dzsdx,f_k=par%f_k)
            now%v_k  = calc_v_kata(now%dtsldy,bnd%dzsdy,f_k=par%f_k)
            now%uv_k = calc_magnitude_from_staggered(now%u_k,now%v_k)

!             ! Combine wind fields 
!             now%ug = now%ug + now%u_k 
!             now%vg = now%vg + now%v_k 
!             now%uvg = calc_magnitude(now%ug,now%vg)

            ! Calculate vertical velocity
            call calc_w(now%ww,now%ug,now%vg,bnd%dzsdx,bnd%dzsdy)

            ! Calculate the surface velocity 
!             now%u_s = calc_u_surf()
!             now%v_s = calc_v_surf()
            now%u_s  = now%ug 
            now%v_s  = now%vg 
            now%uv_s = calc_magnitude_from_staggered(now%u_s,now%v_s)

!             ! Calculate the surface albedo (ice surfaces)
!             als_min  =   0.24
!             als_max  =   0.81
!             afac     =  -0.21
!             tmid     = 274.16
!             where (bnd%f_ice .ge. 0.5) &
!                 now%al_s = calc_albedo_t2m(now%t2m,als_min,als_max,afac,tmid)

!             ! Calculate the surface albedo (non-ice surfaces)
!             als_min  =   0.13
!             als_max  =   1.22
!             afac     =  -0.10
!             tmid     = 270.23
!             where (bnd%f_ice .lt. 0.5) &
!                 now%al_s = calc_albedo_t2m(now%t2m,als_min,als_max,afac,tmid)
        
            als_min  =   0.10
            als_max  =   2.53
            afac     =  -0.021
            tmid     = 278.86
            where (bnd%f_ice .ge. 0.5) &
                now%al_s = calc_albedo_t2m(now%t2m,als_min,als_max,afac,tmid)

            als_min  =   0.10 
            where (bnd%f_ice .lt. 0.5) &
                now%al_s = calc_albedo_t2m(now%t2m,als_min,als_max,afac,tmid)
           
            ! Make sure no strange limits are broken
            where(now%al_s .gt. 0.81) now%al_s = 0.81 
            where(now%al_s .lt. 0.10) now%al_s = 0.10
        
        
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
        
            ! ## Calculate surface fluxes to the atmosphere ###

            ! Shortwave radiation down (surface)
            now%swd_s = calc_radshort_surf_down(now%S,now%cc,bnd%z_srf,&
                                    par%swds_a,par%swds_b,par%swds_c) 

            ! Longwave radiation down (surface)
            now%lwd_s = calc_rad_surf_down(now%t2m,now%cc,now%tcw, &
                                           par%lwds_a,par%lwds_b,par%lwds_c,par%lwds_d)

            ! Calculate latent heat flux at the surface

            ! ajr: TO DO - assume zero for now 
            now%lhf_s = 0.0 

            ! Iterate surface temperature and sensible heat flux calcs
            do iter_surf = 1, n_iter_surf
                
                ! Calculate surface temperature 

                ! To do: calculate from energy balance equation (roughly)

                now%tsurf = now%t2m 
                where (bnd%f_ice .gt. 0.0 .and. now%tsurf .gt. T0) now%tsurf = T0 

                ! Calculate sensible heat flux at the surface
!                 now%shf_s  = 0.0

!                 beta       = 0.0      ! 7-20 W m−2 K−1 (following Krebs-Kanzow etal, tc, 2018 and Braithwaite 1995 et al)
!                 now%shf_s  = beta*(T0-now%t2m)
        
                ! Following Krapp et al., 2017
                now%shf_s  = calc_sensible_heat_flux(now%tsurf,now%t2m,now%uv_s,now%rho_a, &
                                                        csh_pos=1.65e-3*2.4,csh_neg=1.65e-3,cap=cap)

            end do 

            ! Calculate energy balance on low-resolution grid
            ! ajr: need to test whether it's better just to use slightly lower
            ! resolution for the whole atmosphere, but one grid. (High resolution
            ! can happen externally for smb model, where it is needed).
            call rembo_calc_en(now%t2m,emb,par,day,bnd%z_srf,now%t2m_bnd,now%gamma, &
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

!             ! Now calculate the condensation rate (kg m**-2 s**-1)
!             now%c_w = calc_condensation(now%tcw,now%q_r,bnd%dzsdxy,now%ww, &
!                                         par%k_c,par%k_x)
        
            ! Calculate the current cloud water content  (kg m**-2)
            call rembo_calc_ccw(now%pr,now%ccw,now%c_w,emb,par,now%tcw,now%q_r,now%ww)   

            ! Now calculate the high resolution precipitation rate (kg m**-2 s**-1)
            now%pr = calc_precip(now%ccw,now%ww,now%t2m,bnd%z_srf,par%k_w,par%k_z,par%k_t)

            ! Calculate snowfall (kg m**-2 s**-1)
            now%sf = calc_snowfrac(now%t2m,par%sf_a,par%sf_b) * now%pr 

        end do 

        return 

    end subroutine rembo_calc_atmosphere

    subroutine rembo_calc_iterinit(emb,mask,z_srf,t2m_bnd,gamma,par,day)

        implicit none 

        type(diffusion_class),   intent(INOUT) :: emb  
        integer,                 intent(IN)    :: mask(:,:)
        real(wp),                intent(IN)    :: z_srf(:,:)
        real(wp),                intent(IN)    :: t2m_bnd(:,:)
        real(wp),                intent(IN)    :: gamma(:,:)
        type(rembo_param_class), intent(IN)    :: par 
        integer,                 intent(IN)    :: day 
        
        ! Local variables  
        real(wp) :: tsl_fac 
        integer :: q, nx, ny 
        real(dp), allocatable :: tmp8(:,:) 

        nx = emb%grid%G%nx 
        ny = emb%grid%G%ny 

        allocate(tmp8(nx,ny)) 

        ! == emb relaxation mask (used for rembo_calc_en and rembo_calc_ccw)

        ! Relaxation mask (ensure borders are relaxing too)
        call map_field(emb%map_toemb,"mask",mask,emb%mask,method="nn",fill=.TRUE.,missing_value=dble(mv))
        emb%mask(1,:)             = 1 
        emb%mask(emb%grid%G%nx,:) = 1 
        emb%mask(:,1)             = 1 
        emb%mask(:,emb%grid%G%ny) = 1 
        call remove_islands(emb%mask)

        !call map_field(emb%map_toemb,"z_srf",z_srf,emb%z_srf,method="radius",fill=.TRUE.,missing_value=dble(mv))
        call map_field_conservative_map1(emb%map_toemb%map,"z_srf",real(z_srf,dp),tmp8,method="mean",fill=.TRUE.,missing_value=dble(mv))
        emb%z_srf = real(tmp8,wp)

        ! Get the 2D energy diffusion coefficient
        ! emb%kappa = par%en_D 
        emb%kappa = par%en_D_win + (par%en_D_sum-par%en_D_win)*(0.5-0.5*cos((day-15)*2.0*pi/day_year))
        
        ! Boundary sea-level temperature, tsl
        emb%tsl_bnd = mv 
        !call map_field(emb%map_toemb,"tsl_bnd",t2m_bnd+gamma*z_srf,emb%tsl_bnd,method="radius",fill=.TRUE.,missing_value=dble(mv))
        call map_field_conservative_map1(emb%map_toemb%map,"tsl_bnd",real(t2m_bnd+gamma*z_srf,dp),tmp8,method="mean",fill=.TRUE.,missing_value=dble(mv))
        emb%tsl_bnd = real(tmp8,wp)

        ! Ensure borders are populated  with nearest neighbors 
        where (emb%tsl_bnd(2,:) .eq. mv) emb%tsl_bnd(2,:) = emb%tsl_bnd(3,:) 
        where (emb%tsl_bnd(1,:) .eq. mv) emb%tsl_bnd(1,:) = emb%tsl_bnd(2,:) 
        
        where (emb%tsl_bnd(nx-1,:) .eq. mv) emb%tsl_bnd(nx-1,:) = emb%tsl_bnd(nx-2,:) 
        where (emb%tsl_bnd(nx,:)   .eq. mv) emb%tsl_bnd(nx,:)   = emb%tsl_bnd(nx-1,:) 
        
        where (emb%tsl_bnd(:,2) .eq. mv) emb%tsl_bnd(:,2) = emb%tsl_bnd(:,3) 
        where (emb%tsl_bnd(:,1) .eq. mv) emb%tsl_bnd(:,1) = emb%tsl_bnd(:,2) 
        
        where (emb%tsl_bnd(:,ny-1) .eq. mv) emb%tsl_bnd(:,ny-1) = emb%tsl_bnd(:,ny-2) 
        where (emb%tsl_bnd(:,ny)   .eq. mv) emb%tsl_bnd(:,ny)   = emb%tsl_bnd(:,ny-1) 
        
        ! Initialize temperature to boundary field
        emb%tsl = emb%tsl_bnd 

        return 

    end subroutine rembo_calc_iterinit

    subroutine rembo_calc_en(t2m,emb,par,day,z_srf,t2m_bnd,gamma,swn,lwn,shf,lhf,lhp,rco2)

        implicit none 

        real(wp),                intent(INOUT) :: t2m(:,:)
        type(diffusion_class),   intent(INOUT) :: emb  
        real(wp),                intent(IN)    :: z_srf(:,:)
        real(wp),                intent(IN)    :: t2m_bnd(:,:)
        real(wp),                intent(IN)    :: gamma(:,:)
        real(wp),                intent(IN)    :: swn(:,:)
        real(wp),                intent(IN)    :: lwn(:,:)
        real(wp),                intent(IN)    :: shf(:,:)
        real(wp),                intent(IN)    :: lhf(:,:)
        real(wp),                intent(IN)    :: lhp(:,:)
        real(wp),                intent(IN)    :: rco2(:,:)
        type(rembo_param_class), intent(IN)    :: par 
        integer,                 intent(IN)    :: day 
        
        ! Local variables  
        real(wp) :: tsl_fac 
        integer :: q, nx, ny 
        real(dp), allocatable :: tmp8(:,:) 
        real(dp), allocatable :: tmp8hi(:,:) 

        nx = emb%grid%G%nx 
        ny = emb%grid%G%ny 

        allocate(tmp8(nx,ny)) 
        allocate(tmp8hi(size(t2m,1),size(t2m,2))) 

        ! Get the tsl => column energy conversion
        ! tsl_fac = H_a[m] c_p[J kg-1 K-1] rho_a[kg m-3] = [J m-2 K-1]
        tsl_fac = par%en_Ha *1000.0 *1.225 !* 1.225 ! = 8.6e6

        ! ajr: calculated inside of `rembo_calc_iterinit`, not needed here
        !! Boundary sea-level temperature, tsl
        !call map_field(emb%map_toemb,"tsl_bnd",t2m_bnd+gamma*z_srf,emb%tsl_bnd,method="radius",fill=.TRUE.,missing_value=dble(mv))
        
        ! Sea-level temperature, tsl
!         call map_field(emb%map_toemb,"tsl",t2m+gamma*z_srf,emb%tsl,method="radius",fill=.TRUE.,missing_value=dble(mv))
        call map_field_conservative_map1(emb%map_toemb%map,"tsl",real(t2m+gamma*z_srf,dp),tmp8,method="mean",fill=.TRUE.,missing_value=dble(mv))
        emb%tsl = real(tmp8,wp)

        ! Radiative forcing, tsl_F [ J s-1 m-2] * [J-1 m2 K] == [K s-1]
        call map_field(emb%map_toemb,"tsl_F",(swn + lwn + (shf+lhf) + lhp + rco2) / tsl_fac , &
                       emb%tsl_F,method="radius",fill=.TRUE.,missing_value=dble(missing_value))
!         tmp8hi = (swn + lwn + (shf+lhf) + lhp + rco2) / tsl_fac
!         call map_field_conservative_map1(emb%map_toemb%map,"tsl_F",tmp8hi, &
!                         tmp8,method="mean",fill=.TRUE.,missing_value=dble(mv))
!         emb%tsl_F = real(tmp8,wp)

        ! Calculate radiative balance over the day
        do q = 1, par%en_nstep * 5
            !where (emb%mask .eq. 1) emb%en = emb%en_bnd
            call adv_diff_2D(emb%tsl,emb%tsl_bnd,emb%tsl_F,relax=emb%mask, &
                             dx=real(emb%grid%G%dx*emb%grid%xy_conv,wp), &
                             dy=real(emb%grid%G%dx*emb%grid%xy_conv,wp), &
                             dt=par%en_dt,kappa=emb%kappa,k_relax=par%en_kr) !, &
!                              v_x=emb%ug,v_y=emb%vg)
!             call solve_diff_2D_adi(emb%tsl,emb%tsl_bnd,emb%tsl_F,relax=emb%mask, &
!                                    dx=real(emb%grid%G%dx*emb%grid%xy_conv,wp), &
!                                    dy=real(emb%grid%G%dx*emb%grid%xy_conv,wp), &
!                                    dt=par%en_dt,kappa=emb%kappa,k_relax=par%en_kr)

        end do 

        call map_field(emb%map_fromemb,"tsl",emb%tsl,t2m,method="nng",fill=.TRUE.,missing_value=dble(mv),sigma=50.d0)
!         call map_field_conservative_map1(emb%map_fromemb%map,"tsl",real(emb%tsl,dp),tmp8hi,method="mean",fill=.TRUE.,missing_value=dble(mv))
!         t2m = real(tmp8hi,wp)
        t2m = t2m - gamma*z_srf

!         ! ajr: TO DO !!  (seems to go really slow, or smt is wrong)
!         t2m =  interp_bilinear(x=emb%grid%G%x,y=emb%grid%G%y,z=dble(emb%tsl), &
!                                xout=emb%map_fromemb%x,yout=emb%map_fromemb%y, &
!                                missing_value=dble(mv))
!         t2m = t2m - gamma*z_srf

        return 

    end subroutine rembo_calc_en

    subroutine rembo_calc_ccw(pr,ccw,c_w,emb,par,tcw,q_r,ww)

        implicit none 

        real(wp),                intent(INOUT) :: pr(:,:)
        real(wp),                intent(INOUT) :: ccw(:,:)
        real(wp),                intent(INOUT) :: c_w(:,:)
        type(diffusion_class),   intent(INOUT) :: emb 
        type(rembo_param_class), intent(IN)    :: par 
        real(wp),                intent(IN)    :: tcw(:,:)
        real(wp),                intent(IN)    :: q_r(:,:)
        real(wp),                intent(IN)    :: ww(:,:)

        ! Local variables
        real(wp)  :: dtsl_mean, lat0, lat1, sigma 
        integer   :: q  
        !real(wp)  :: dtsl_mean    ! Average temp anomaly away from boundary temps 

        !emb%tcw     = calc_qs(emb%tsl_bnd+dtsl_mean-par%gamma*emb%z_srf,emb%z_srf,par%e0,par%c1) *emb%rho_a *par%H_e
!         emb%tcw_bnd = calc_qs(emb%tsl_bnd-par%gamma*emb%z_srf,emb%z_srf,par%e0,par%c1) *emb%rho_a *par%H_e
!         emb%tcw     = emb%tcw_bnd
    
        !emb%q_r = emb%tcw / (calc_qsat(emb%tsl_bnd+dtsl_mean-par%gamma*emb%z_srf,emb%z_srf) *emb%rho_a *par%H_e)

        ! == Total water content ==
        call map_field(emb%map_toemb,"tcw",tcw,emb%tcw,method="radius",fill=.TRUE.,missing_value=dble(mv))
        
        ! == Relative humidity == 
        call map_field(emb%map_toemb,"q_r",q_r,emb%q_r,method="radius",fill=.TRUE.,missing_value=dble(mv))
        
        
        ! Get 2D diffusion constant
        emb%kappaw = par%ccw_D

        ! Note: removed any lat/elevation dependence of diffusion coefficient (ajr, 2015-08-31)
!         lat0 = minval(emb%grid%lat)
!         lat1 = maxval(emb%grid%lat)
!         emb%kappaw = emb%kappaw* & 
!             (1.d0 + par%en_kl*(2.d0*(emb%grid%lat-lat0)/(lat1-lat0)-1.d0))

!         emb%kappaw = emb%kappaw* (1.d0 + par%en_kz*emb%z_srf)
    
!         call map_field(emb%map_toemb,"en_F",swn + lwn + (shf+lhf) + lhp + rco2, &
!                        emb%en_F,method="radius",fill=.TRUE.,missing_value=dble(missing_value))

        ! Get vertical wind for precipitation
        sigma = emb%grid%G%dx*2.0  
        call map_field(emb%map_toemb,"ww",ww,emb%ww,method="nng",sigma=dble(sigma),fill=.TRUE.,missing_value=dble(mv))
!         call map_field(emb%map_toemb,"ww",ww,emb%ccw_pr,method="nng",sigma=dble(sigma),fill=.TRUE.,missing_value=dble(mv))
!         emb%ww = (emb%ww+emb%ccw_pr)/2.d0 
        ! ajr: not sure what the above average represented in rembo2-beta. For now just use ww directly...

        ! Get current moisture content [kg m-2]
        call map_field(emb%map_toemb,"ccw",ccw,emb%ccw,method="nng",sigma=dble(sigma),fill=.TRUE.,missing_value=dble(mv))
        
        ! ajr: to do: define ccw_bnd!!!
        emb%ccw_bnd = emb%ccw 

        ! Calculate the initial lo-res precipitation rate [kg m-2 s-1]
        call map_field(emb%map_toemb,"pr", pr, emb%ccw_pr,method="nng",sigma=dble(sigma),fill=.TRUE.,missing_value=dble(mv))
        ! ajr: TO DO : calculate a basic rate somehow avoiding history 


        ! Now calculate the condensation rate [kg m-2 s-1]
        emb%ccw_cw = calc_condensation(emb%tcw,emb%q_r,emb%ww,par%k_c,par%k_x)

        ! Get moisture balance forcing [kg m-2 s-1]
        emb%ccw_F = emb%ccw_cw - emb%ccw_pr 
        where(emb%mask==1) emb%ccw_F = 0.d0 

        ! Calculate moisture balance to equilibrium
        do q = 1, par%ccw_nstep
            !where (emb%mask .eq. 1) emb%ccw = emb%ccw_bnd
            call adv_diff_2D(emb%ccw,emb%ccw_bnd, & !*emb%tcw/emb%tcw_bnd, &
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

            ! Get moisture balance forcing [kg m-2 s-1]
            emb%ccw_F = emb%ccw_cw - emb%ccw_pr 
            where(emb%mask==1) emb%ccw_F = 0.d0 

        end do 

        ! Send cloud moisture content back to main domain pts
        sigma = emb%map_fromemb%G%dx*2.0 
        call map_field(emb%map_fromemb,"ccw",emb%ccw,ccw,method="nng",sigma=dble(sigma),fill=.TRUE.,missing_value=dble(mv))
        call map_field(emb%map_fromemb,"c_w",emb%ccw_cw,c_w,method="nng",sigma=dble(sigma),fill=.TRUE.,missing_value=dble(mv))
        
        ! ajr: TO DO 
!         ccw =  interp_bilinear(is_points=.TRUE.,x=emb%grid%G%x,y=emb%grid%G%y,z=emb%ccw, &
!                                xout=emb%map_fromemb%x,yout=emb%map_fromemb%y, &
!                                missing_value=missing_value)

        ! ajr: TO DO 
!         c_w = interp_bilinear(is_points=.TRUE.,x=emb%grid%G%x,y=emb%grid%G%y,z=emb%ccw_cw, &
!                                xout=emb%map_fromemb%x,yout=emb%map_fromemb%y, &
!                                missing_value=missing_value)
         
        ! ajr: TO DO 
!         pr = interp_bilinear(is_points=.TRUE.,x=emb%grid%G%x,y=emb%grid%G%y,z=emb%ccw_pr, &
!                                xout=emb%map_fromemb%x,yout=emb%map_fromemb%y, &
!                                missing_value=missing_value)

        return 

    end subroutine rembo_calc_ccw

end module rembo_atm
