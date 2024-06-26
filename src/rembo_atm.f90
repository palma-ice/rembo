module rembo_atm 
    
    use gaussian_filter 
    use solvers 

    use rembo_defs 
    use rembo_physics 
    use rembo1

    use coordinates_mapping_scrip, only : map_scrip_class, map_scrip_init, map_scrip_field, &
                                            gen_map_filename, nc_read_interp

    implicit none 

    real(wp), parameter :: cap  = 1000.0_wp  !< specific heat capacitiy of air [J/(kg K)]
    ! ajr: move to rembo_defs later...

    private
    public :: rembo1_calc_atmosphere
    public :: rembo_calc_atmosphere
    public :: rembo_gen_domain_mask

contains 

    subroutine rembo1_calc_atmosphere(now,emb,bnd,grid,par,day,year)
        ! Calculate the equilibrium rembo atmosphere solution 
        ! for a given day of a given year (can be representative of a month)

        implicit none 

        type(rembo_state_class),    intent(INOUT) :: now     ! rembo atmosphere variables
        type(diffusion_class),      intent(INOUT) :: emb     ! Diffusion variables
        type(rembo_boundary_class), intent(IN)    :: bnd     ! rembo boundary variables
        type(rgrid_class),          intent(IN)    :: grid    ! rembo grid information
        type(rembo_param_class),    intent(IN)    :: par     ! rembo parameters
        integer, intent(IN) :: day 
        integer, intent(IN) :: year 

        ! Local variables
        real(wp) :: als_max, als_min, afac, tmid    ! ajr: to do: move to parameters!!
        real(wp) :: beta                            ! ajr: to do: move to parameters!!
        integer  :: iter, n_iter 
        integer  :: iter_surf, n_iter_surf 
        integer  :: iter_rembo1, n_iter_rembo1

        real(wp), allocatable :: dh_snow(:,:)
        real(wp), allocatable :: dzs(:,:)

        allocate(dh_snow(grid%nx,grid%ny))
        allocate(dzs(grid%nx,grid%ny))

        ! Get the current lapse rate (gamma=winter, gamma2=summer)  <= CHECK !!
        now%gamma = par%gamma  &
           + (0.5*cos((day-15)*2.0*pi/par%c%day_year)-0.5)*(par%gamma-par%gamma2)
        !write(*,*) "gamma: ", day, minval(gamma), maxval(gamma)
        
        ! Get the surface pressure 
        now%sp = calc_sp(bnd%z_srf,par%c%g)
        
        ! Get the total column mass
        call calc_total_column_mass(now%tcm,bnd%z_srf,rho_0=1.3_wp,H_a=8e3_wp,H_toa=20e3_wp)
        
        ! Calc gradient of Z: dZdx, dZdy 
        call d_dx(now%dZdx,now%Z,dx=grid%dx)
        call d_dy(now%dZdy,now%Z,dx=grid%dy)

        ! Get geostrophic components and magnitude of velocity
        now%ug  = calc_u_geo(now%dZdy,bnd%f,par%c%g)
        now%vg  = calc_v_geo(now%dZdx,bnd%f,par%c%g)
        now%uvg = calc_magnitude(now%ug,now%vg)

        ! Initialize t2m solution with boundary solution 
        now%t2m = now%t2m_bnd 
        
        ! Initialize emb object for this day 
        call rembo_calc_iterinit(emb,bnd%mask,bnd%z_srf,now%t2m_bnd,now%ccw_bnd,now%gamma,par,day,par%c%g)

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

            ! Calculate the surface albedo (ice surfaces)
            als_min  =   0.24
            als_max  =   0.81
            afac     =  -0.21
            tmid     = 274.16
            where (bnd%f_ice .ge. 0.5) &
                now%al_s = calc_albedo_t2m(now%t2m,als_min,als_max,afac,tmid)

            ! Calculate the surface albedo (non-ice surfaces)
            als_min  =   0.13
            als_max  =   1.22
            afac     =  -0.10
            tmid     = 270.23
            where (bnd%f_ice .lt. 0.5) &
                now%al_s = calc_albedo_t2m(now%t2m,als_min,als_max,afac,tmid)
        
            ! als_min  =   0.10
            ! als_max  =   2.53
            ! afac     =  -0.021
            ! tmid     = 278.86
            ! where (bnd%f_ice .ge. 0.5) &
            !     now%al_s = calc_albedo_t2m(now%t2m,als_min,als_max,afac,tmid)

            ! als_min  =   0.10 
            ! where (bnd%f_ice .lt. 0.5) &
            !     now%al_s = calc_albedo_t2m(now%t2m,als_min,als_max,afac,tmid)
           
            ! Make sure no strange limits are broken
            where(now%al_s .gt. 0.81) now%al_s = 0.81 
            where(now%al_s .lt. 0.10) now%al_s = 0.10

        
            ! Calculate the planetary albedo 
            now%al_p = par%alp_a + par%alp_b*now%al_s - par%alp_c*now%tcw 
            where(now%al_p .lt. 0.0) now%al_p = 0.0
            where(now%al_p .gt. 1.0) now%al_p = 1.0  

            ! ajr: testing strong albedo reduction for ice-free points
            where (bnd%f_ice .lt. 0.5) now%al_p = 0.1

            ! Calculate the outgoing long-wave radiation at toa
            !now%lwu = par%lwu_a + par%lwu_b*(now%t2m-par%c%T0) + par%lwu_c*now%S
            !now%lwu = par%lwu_a + par%lwu_b*(now%t2m-par%c%T0)
            now%lwu = par%lwu_a + par%lwu_b*now%t2m
            
            ! Calculate the incoming short-wave radiation at toa
            !now%swd = (1.0-now%al_p)*now%S 
            now%swd = now%S 

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
            where(bnd%z_srf .gt. 0.0 .and. bnd%f_ice .lt. 0.5) now%lhf_s = 50.0 

            ! Iterate surface temperature and sensible heat flux calcs
            do iter_surf = 1, n_iter_surf
                
                ! Calculate surface temperature 

                ! To do: calculate from energy balance equation (roughly)

                now%tsurf = now%t2m 
                where (bnd%f_ice .gt. 0.0 .and. now%tsurf .gt. par%c%T0) now%tsurf = par%c%T0 

                ! Calculate sensible heat flux at the surface
!                 now%shf_s  = 0.0

!                 beta       = 0.0      ! 7-20 W m−2 K−1 (following Krebs-Kanzow etal, tc, 2018 and Braithwaite 1995 et al)
!                 now%shf_s  = beta*(T0-now%t2m)
        
                ! Following Krapp et al., 2017
                now%shf_s  = calc_sensible_heat_flux(now%tsurf,now%t2m,now%uv_s,now%rho_a, &
                                                        csh_pos=1.65e-3*2.4,csh_neg=1.65e-3,cap=cap)

                ! ajr: testing
                now%shf_s = 50.0 
                
            end do 

if (.TRUE.) then
            ! Calculate energy balance on low-resolution grid
            ! ajr: need to test whether it's better just to use slightly lower
            ! resolution for the whole atmosphere, but one grid. (High resolution
            ! can happen externally for smb model, where it is needed).
            ! call rembo_calc_en(now%t2m,emb,par,day,bnd%z_srf,now%t2m_bnd,now%gamma, &
            !                   swn=(now%swd - now%swd_s*(1.0-now%al_s)), &
            !                   lwn=(- now%lwu + now%lwu_s - now%lwd_s), &
            !                   shf=now%shf_s,lhf=now%lhf_s, &
            !                   lhp=(par%Lw*(now%pr-now%sf) + par%Ls*now%sf)*1.0, &
            !                   rco2=par%en_kdT + now%rco2_a, &
            !                   ug=now%ug,vg=now%vg,dx=grid%dx,g=par%c%g) 

            ! call rembo_calc_en(now%t2m,emb,par,day,bnd%z_srf,now%t2m_bnd,now%gamma, &
            !                   swn=now%swd*(1.0-now%al_p), &
            !                   lwn=-now%lwu, &
            !                   shf=now%shf_s*0.0_wp,lhf=now%lhf_s*0.0_wp, &
            !                   lhp=(par%Lw*(now%pr-now%sf) + par%Ls*now%sf)*1.0_wp, &
            !                   rco2=par%en_kdT + now%rco2_a, &
            !                   ug=now%ug,vg=now%vg,dx=grid%dx,g=par%c%g) 

            call rembo_calc_en_simple(now%t2m,now%tsl,par,day,bnd%z_srf, &
                              now%t2m_bnd,now%tsl_bnd,now%gamma, &
                              swn=now%swd*(1.0-now%al_p), &
                              lwn=-now%lwu, &
                              shf=now%shf_s*0.0_wp,lhf=now%lhf_s, &
                              lhp=(par%Lw*(now%pr-now%sf) + par%Ls*now%sf)*1.0_wp, &
                              !lhp=now%pr*0.0, &
                              rco2=par%en_kdT + now%rco2_a, &
                              ug=now%ug*0.0,vg=now%vg*0.0, &
                              mask=bnd%mask, &
                              dx=grid%dx,g=par%c%g)

            
!else 
!            ! Impose boundary temperature field for now
!            now%t2m = now%t2m_bnd 

end if 

            ! Calculate inversion correction for moisture balance
            now%ct2m = calc_ttcorr1(now%t2m,bnd%z_srf,-2.4e-3,-3.7e-1,106.0)
            
            ! Get saturated specific humidity and sat. total water content
            now%q_sat   = calc_qsat(now%t2m+now%ct2m,bnd%z_srf,par%c%g,par%c%T0) 
            now%tcw_sat = now%q_sat * now%rho_a * par%H_e

            ! Calculate the current total water content
            now%q_s = calc_qs(now%t2m+now%ct2m,bnd%z_srf,par%e0,par%c1,par%c%g,par%c%T0) 
            now%tcw = now%q_s * now%rho_a * par%H_e
            now%q_r = now%tcw / now%tcw_sat 

            ! Now calculate the condensation rate (kg m**-2 s**-1)
            now%c_w = calc_condensation(now%tcw,now%q_r,bnd%dzsdxy,now%ww, &
                                        par%k_c,par%k_x)

            ! ajr: disabled condensation for now, while testing energy
            !now%c_w = 0.0

if (.FALSE.) then
            ! Calculate the current cloud water content  (kg m**-2)
            call rembo_calc_ccw(now%pr,now%ccw,now%c_w,emb,par,now%tcw,now%q_r,now%ww,grid%dx) 

            ! Now calculate the high resolution precipitation rate (kg m**-2 s**-1)
            now%pr = calc_precip(now%ccw,now%ww,now%t2m,bnd%z_srf,par%k_w,par%k_z,par%k_t)  

            ! Calculate snowfall (kg m**-2 s**-1)
            now%sf = calc_snowfrac(now%t2m,par%sf_a,par%sf_b) * now%pr 
end if 

if (.FALSE.) then
    ! REMBO1 - old code!! Not working well yet....

            n_iter_rembo1 = 20 

            do iter_rembo1 = 1, n_iter_rembo1

                ! Get rate of snowfall melt (currently zero)
                dh_snow = 0.0

                ! Call REMBO1 atmosphere
                call rembo1_calc_atm(now%t2m,now%tcw,now%pr,now%pr-now%sf,now%sf, &
                                        now%t2m_bnd,now%q_r,now%S,now%al_p,dh_snow, &
                                        par%en_kdT+now%rco2_a,bnd%z_srf,bnd%dzsdxy, &
                                        bnd%mask,grid%dx,par%en_dt,par%c%sec_day)

            end do
end if 

        end do 

        return 

    end subroutine rembo1_calc_atmosphere

    subroutine rembo_calc_atmosphere(now,emb,bnd,grid,par,day,year)
        ! Calculate the equilibrium rembo atmosphere solution 
        ! for a given day of a given year (can be representative of a month)

        implicit none 

        type(rembo_state_class),    intent(INOUT) :: now     ! rembo atmosphere variables
        type(diffusion_class),      intent(INOUT) :: emb     ! Diffusion variables
        type(rembo_boundary_class), intent(IN)    :: bnd     ! rembo boundary variables
        type(rgrid_class),          intent(IN)    :: grid    ! rembo grid information
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
           + (0.5*cos((day-15)*2.0*pi/par%c%day_year)-0.5)*(par%gamma-par%gamma2)
        !write(*,*) "gamma: ", day, minval(gamma), maxval(gamma)
        
        ! Get the surface pressure 
        now%sp = calc_sp(bnd%z_srf,par%c%g)
        
        ! Calc gradient of Z: dZdx, dZdy 
        call d_dx(now%dZdx,now%Z,dx=grid%dx)
        call d_dy(now%dZdy,now%Z,dx=grid%dy)

        ! Get geostrophic components and magnitude of velocity
        now%ug  = calc_u_geo(now%dZdy,bnd%f,par%c%g)
        now%vg  = calc_v_geo(now%dZdx,bnd%f,par%c%g)
        now%uvg = calc_magnitude(now%ug,now%vg)

        ! Initialize t2m solution with boundary solution 
        now%t2m = now%t2m_bnd 
        
        ! Initialize emb object for this day 
        call rembo_calc_iterinit(emb,bnd%mask,bnd%z_srf,now%t2m_bnd,now%ccw_bnd,now%gamma,par,day,par%c%g)

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

            ! Calculate the surface albedo (ice surfaces)
            als_min  =   0.24
            als_max  =   0.81
            afac     =  -0.21
            tmid     = 274.16
            where (bnd%f_ice .ge. 0.5) &
                now%al_s = calc_albedo_t2m(now%t2m,als_min,als_max,afac,tmid)

            ! Calculate the surface albedo (non-ice surfaces)
            als_min  =   0.13
            als_max  =   1.22
            afac     =  -0.10
            tmid     = 270.23
            where (bnd%f_ice .lt. 0.5) &
                now%al_s = calc_albedo_t2m(now%t2m,als_min,als_max,afac,tmid)
        
            ! als_min  =   0.10
            ! als_max  =   2.53
            ! afac     =  -0.021
            ! tmid     = 278.86
            ! where (bnd%f_ice .ge. 0.5) &
            !     now%al_s = calc_albedo_t2m(now%t2m,als_min,als_max,afac,tmid)

            ! als_min  =   0.10 
            ! where (bnd%f_ice .lt. 0.5) &
            !     now%al_s = calc_albedo_t2m(now%t2m,als_min,als_max,afac,tmid)
           
            ! Make sure no strange limits are broken
            where(now%al_s .gt. 0.81) now%al_s = 0.81 
            where(now%al_s .lt. 0.10) now%al_s = 0.10

        
            ! Calculate the planetary albedo 
            now%al_p = par%alp_a + par%alp_b*now%al_s - par%alp_c*now%tcw 
            where(now%al_p .lt. 0.0) now%al_p = 0.0
            where(now%al_p .gt. 1.0) now%al_p = 1.0  

            ! Calculate the outgoing long-wave radiation at toa
            !now%lwu = par%lwu_a + par%lwu_b*(now%t2m-par%c%T0) + par%lwu_c*now%S
            !now%lwu = par%lwu_a + par%lwu_b*(now%t2m-par%c%T0)
            now%lwu = par%lwu_a + par%lwu_b*now%t2m
            
            ! Calculate the incoming short-wave radiation at toa
            !now%swd = (1.0-now%al_p)*now%S 
            now%swd = now%S 

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
            where(bnd%z_srf .gt. 0.0 .and. bnd%f_ice .lt. 0.5) now%lhf_s = 50.0 

            ! Iterate surface temperature and sensible heat flux calcs
            do iter_surf = 1, n_iter_surf
                
                ! Calculate surface temperature 

                ! To do: calculate from energy balance equation (roughly)

                now%tsurf = now%t2m 
                where (bnd%f_ice .gt. 0.0 .and. now%tsurf .gt. par%c%T0) now%tsurf = par%c%T0 

                ! Calculate sensible heat flux at the surface
!                 now%shf_s  = 0.0

!                 beta       = 0.0      ! 7-20 W m−2 K−1 (following Krebs-Kanzow etal, tc, 2018 and Braithwaite 1995 et al)
!                 now%shf_s  = beta*(T0-now%t2m)
        
                ! Following Krapp et al., 2017
                now%shf_s  = calc_sensible_heat_flux(now%tsurf,now%t2m,now%uv_s,now%rho_a, &
                                                        csh_pos=1.65e-3*2.4,csh_neg=1.65e-3,cap=cap)

            end do 

if (.TRUE.) then
            ! Calculate energy balance on low-resolution grid
            ! ajr: need to test whether it's better just to use slightly lower
            ! resolution for the whole atmosphere, but one grid. (High resolution
            ! can happen externally for smb model, where it is needed).
            ! call rembo_calc_en(now%t2m,emb,par,day,bnd%z_srf,now%t2m_bnd,now%gamma, &
            !                   swn=(now%swd - now%swd_s*(1.0-now%al_s)), &
            !                   lwn=(- now%lwu + now%lwu_s - now%lwd_s), &
            !                   shf=now%shf_s,lhf=now%lhf_s, &
            !                   lhp=(par%Lw*(now%pr-now%sf) + par%Ls*now%sf)*1.0, &
            !                   rco2=par%en_kdT + now%rco2_a, &
            !                   ug=now%ug,vg=now%vg,dx=grid%dx,g=par%c%g) 

            ! call rembo_calc_en(now%t2m,emb,par,day,bnd%z_srf,now%t2m_bnd,now%gamma, &
            !                   swn=now%swd*(1.0-now%al_p), &
            !                   lwn=-now%lwu, &
            !                   shf=now%shf_s*0.0_wp,lhf=now%lhf_s*0.0_wp, &
            !                   lhp=(par%Lw*(now%pr-now%sf) + par%Ls*now%sf)*1.0_wp, &
            !                   rco2=par%en_kdT + now%rco2_a, &
            !                   ug=now%ug,vg=now%vg,dx=grid%dx,g=par%c%g) 

            call rembo_calc_en_simple(now%t2m,now%tsl,par,day,bnd%z_srf, &
                              now%t2m_bnd,now%tsl_bnd,now%gamma, &
                              swn=now%swd*(1.0-now%al_p), &
                              lwn=-now%lwu, &
                              shf=now%shf_s*0.0_wp,lhf=now%lhf_s, &
                              lhp=(par%Lw*(now%pr-now%sf) + par%Ls*now%sf)*1.0_wp, &
                              rco2=par%en_kdT + now%rco2_a, &
                              ug=now%ug,vg=now%vg, &
                              mask=bnd%mask, &
                              dx=grid%dx,g=par%c%g)

            
else 
            ! Impose boundary temperature field for now
            now%t2m = now%t2m_bnd 

end if 

            ! Calculate inversion correction for moisture balance
            now%ct2m = calc_ttcorr1(now%t2m,bnd%z_srf,-2.4e-3,-3.7e-1,106.0)
            
            ! Get saturated specific humidity and sat. total water content
            now%q_sat   = calc_qsat(now%t2m+now%ct2m,bnd%z_srf,par%c%g,par%c%T0) 
            now%tcw_sat = now%q_sat * now%rho_a * par%H_e

            ! Calculate the current total water content
            now%q_s = calc_qs(now%t2m+now%ct2m,bnd%z_srf,par%e0,par%c1,par%c%g,par%c%T0) 
            now%tcw = now%q_s * now%rho_a * par%H_e
            now%q_r = now%tcw / now%tcw_sat 

!             ! Now calculate the condensation rate (kg m**-2 s**-1)
!             now%c_w = calc_condensation(now%tcw,now%q_r,bnd%dzsdxy,now%ww, &
!                                         par%k_c,par%k_x)

if (.TRUE.) then
            ! Calculate the current cloud water content  (kg m**-2)
            call rembo_calc_ccw(now%pr,now%ccw,now%c_w,emb,par,now%tcw,now%q_r,now%ww,grid%dx)   
end if 

            ! Now calculate the high resolution precipitation rate (kg m**-2 s**-1)
            now%pr = calc_precip(now%ccw,now%ww,now%t2m,bnd%z_srf,par%k_w,par%k_z,par%k_t)

            ! Calculate snowfall (kg m**-2 s**-1)
            now%sf = calc_snowfrac(now%t2m,par%sf_a,par%sf_b) * now%pr 

        end do 

        return 

    end subroutine rembo_calc_atmosphere

    subroutine rembo_calc_iterinit(emb,mask,z_srf,t2m_bnd,ccw_bnd,gamma,par,day,g)

        implicit none 

        type(diffusion_class),   intent(INOUT) :: emb  
        integer,                 intent(IN)    :: mask(:,:)
        real(wp),                intent(IN)    :: z_srf(:,:)
        real(wp),                intent(IN)    :: t2m_bnd(:,:)
        real(wp),                intent(IN)    :: ccw_bnd(:,:)
        real(wp),                intent(IN)    :: gamma(:,:)
        type(rembo_param_class), intent(IN)    :: par 
        integer,                 intent(IN)    :: day 
        real(wp),                intent(IN)    :: g 

        ! Local variables  
        real(wp) :: tsl_fac 
        integer :: q, nx, ny

        nx = emb%grid%nx 
        ny = emb%grid%ny 

        ! == emb relaxation mask (used for rembo_calc_en and rembo_calc_ccw)

        ! Relaxation mask (ensure borders are relaxing too)
        !call map_field(emb%map_toemb,"mask",mask,emb%mask,method="nn",fill=.TRUE.,missing_value=dble(mv))
        call map_scrip_field(emb%map_toemb,"mask",mask,emb%mask,method="count",missing_value=int(mv))

        emb%mask(1,:)           = -1 
        emb%mask(emb%grid%nx,:) = -1 
        emb%mask(:,1)           = -1 
        emb%mask(:,emb%grid%ny) = -1 
        call remove_islands(emb%mask)

        ! Surface elevation, z_srf
        ! Interpolate to low resolution emb domain
        call map_scrip_field(emb%map_toemb,"z_srf",z_srf,emb%z_srf,method="mean", &
                                        missing_value=real(mv,wp),fill_method="weighted")

        ! Get the tsl => column energy conversion
        ! tsl_fac = H_a[m] c_v[J kg-1 K-1] rho_a[kg m-3] = [J m-2 K-1]
        ! H_a = 8000 m
        !tsl_fac = 8000.0 *715.0 *1.225 !* 1.225 ! =~ 8e6
        
        ! climber-x formula to get tsl_fac
        ! tsl_fac = slp / g * cv_atm 
        ! slp = 101100.0 [kg m-1 s-2]; g = 9.81 [m s-2]; c_v_atm = 715.0 [J kg-1 K-1]
        ! tsl_fac => [kg m-1 s-2] * [m-1 s2] * [J kg-1 K-1] = [J m-2 K-1]
        tsl_fac = 101100.0 /g *715.0 

        ! Get the 2D energy diffusion coefficient
        ! emb%kappa_t = par%en_D 
        emb%kappa_t = par%en_D_win + (par%en_D_sum-par%en_D_win)*(0.5-0.5*cos((day-15)*2.0*pi/par%c%day_year))
        emb%kappa_t = emb%kappa_t / tsl_fac 

        ! Boundary sea-level temperature, tsl
        ! Interpolate to low resolution emb domain
        emb%tsl_bnd = mv 
        call map_scrip_field(emb%map_toemb,"tsl_bnd",t2m_bnd+gamma*z_srf,emb%tsl_bnd,method="mean", &
                                                        missing_value=real(mv,wp),fill_method="weighted")

        ! Boundary total column water vapor, ccw
        emb%ccw_bnd = mv 
        call map_scrip_field(emb%map_toemb,"ccw_bnd",ccw_bnd,emb%ccw_bnd,method="mean", &
                                                        missing_value=real(mv,wp),fill_method="weighted")

        ! ajr: testing disabling ccw_bnd
        emb%ccw_bnd = emb%ccw_bnd*0.01

        ! Initialize temperature and ccw to boundary field
        emb%tsl = emb%tsl_bnd 
        emb%ccw = emb%ccw_bnd 

        return 

    end subroutine rembo_calc_iterinit

    subroutine rembo_calc_en_simple(t2m,tsl,par,day,z_srf,t2m_bnd,tsl_bnd,gamma,swn,lwn,shf,lhf,lhp,rco2,ug,vg,mask,dx,g)

        implicit none 

        real(wp),                intent(INOUT) :: t2m(:,:)  
        real(wp),                intent(INOUT) :: tsl(:,:)  
        type(rembo_param_class), intent(IN)    :: par 
        integer,                 intent(IN)    :: day 
        real(wp),                intent(IN)    :: z_srf(:,:)
        real(wp),                intent(IN)    :: t2m_bnd(:,:)
        real(wp),                intent(IN)    :: tsl_bnd(:,:)
        real(wp),                intent(IN)    :: gamma(:,:)
        real(wp),                intent(IN)    :: swn(:,:)
        real(wp),                intent(IN)    :: lwn(:,:)
        real(wp),                intent(IN)    :: shf(:,:)
        real(wp),                intent(IN)    :: lhf(:,:)
        real(wp),                intent(IN)    :: lhp(:,:)
        real(wp),                intent(IN)    :: rco2(:,:)
        real(wp),                intent(IN)    :: ug(:,:)
        real(wp),                intent(IN)    :: vg(:,:)
        integer,                 intent(IN)    :: mask(:,:)
        real(wp),                intent(IN)    :: dx
        real(wp),                intent(IN)    :: g 

        ! Local variables  
        real(wp) :: tsl_fac 
        integer :: q, nx, ny

        real(wp), allocatable :: en(:,:)
        real(wp), allocatable :: en_F(:,:)
        real(wp), allocatable :: kappa(:,:)

        nx = size(t2m,1)
        ny = size(t2m,2)

        allocate(en(nx,ny))
        allocate(en_F(nx,ny))
        allocate(kappa(nx,ny))
        
        ! Get the tsl => column energy conversion
        ! tsl_fac = H_a[m] c_v[J kg-1 K-1] rho_a[kg m-3] = [J m-2 K-1]
        ! H_a = 8000 m
        !tsl_fac = 8000.0 *715.0 *1.225 !* 1.225 ! =~ 8e6
        
        ! climber-x formula to get tsl_fac
        ! tsl_fac = slp / g * cv_atm 
        ! slp = 101100.0 [kg m-1 s-2]; g = 9.81 [m s-2]; c_v_atm = 715.0 [J kg-1 K-1]
        ! tsl_fac => [kg m-1 s-2] * [m-1 s2] * [J kg-1 K-1] = [J m-2 K-1]
        tsl_fac = 101100.0 /g *715.0 

        ! Sea-level temperature, tsl
        tsl = t2m+gamma*z_srf

        ! Radiative forcing, en_F [ J s-1 m-2] * [J-1 m2 K] == [K s-1]
        en_F = (swn + lwn + (shf+lhf) + lhp + rco2) / tsl_fac

        ! Energy diffusion coefficient
        kappa = par%en_D_win + (par%en_D_sum-par%en_D_win)*(0.5-0.5*cos((day-15)*2.0*pi/par%c%day_year))
        kappa = kappa / tsl_fac 

        ! Calculate radiative balance over the day
        do q = 1, par%en_nstep * 5
            call solve_diffusion_advection_2D(tsl,ug,vg,en_F,kappa,tsl_bnd,mask,dx,dx,par%en_dt, &
                                    k_rel=par%en_kr,solver=par%solver,step=par%step,bc="infinite")
        end do 

        ! 2m temperature
        t2m = tsl - gamma*z_srf

        return 

    end subroutine rembo_calc_en_simple
    
    subroutine rembo_calc_en(t2m,emb,par,day,z_srf,t2m_bnd,gamma,swn,lwn,shf,lhf,lhp,rco2,ug,vg,dx,g)

        implicit none 

        real(wp),                intent(INOUT) :: t2m(:,:)
        type(diffusion_class),   intent(INOUT) :: emb  
        type(rembo_param_class), intent(IN)    :: par 
        integer,                 intent(IN)    :: day 
        real(wp),                intent(IN)    :: z_srf(:,:)
        real(wp),                intent(IN)    :: t2m_bnd(:,:)
        real(wp),                intent(IN)    :: gamma(:,:)
        real(wp),                intent(IN)    :: swn(:,:)
        real(wp),                intent(IN)    :: lwn(:,:)
        real(wp),                intent(IN)    :: shf(:,:)
        real(wp),                intent(IN)    :: lhf(:,:)
        real(wp),                intent(IN)    :: lhp(:,:)
        real(wp),                intent(IN)    :: rco2(:,:)
        real(wp),                intent(IN)    :: ug(:,:)
        real(wp),                intent(IN)    :: vg(:,:)
        real(wp),                intent(IN)    :: dx
        real(wp),                intent(IN)    :: g 

        ! Local variables  
        real(wp) :: tsl_fac 
        integer :: q, nx, ny 

        nx = emb%grid%nx 
        ny = emb%grid%ny 
        
        ! Get the tsl => column energy conversion
        ! tsl_fac = H_a[m] c_v[J kg-1 K-1] rho_a[kg m-3] = [J m-2 K-1]
        ! H_a = 8000 m
        !tsl_fac = 8000.0 *715.0 *1.225 !* 1.225 ! =~ 8e6
        
        ! climber-x formula to get tsl_fac
        ! tsl_fac = slp / g * cv_atm 
        ! slp = 101100.0 [kg m-1 s-2]; g = 9.81 [m s-2]; c_v_atm = 715.0 [J kg-1 K-1]
        ! tsl_fac => [kg m-1 s-2] * [m-1 s2] * [J kg-1 K-1] = [J m-2 K-1]
        tsl_fac = 101100.0 /g *715.0 

        ! Note: tsl_bnd has been updated in `rembo_calc_iterinit`

        ! Sea-level temperature, tsl
        call map_scrip_field(emb%map_toemb,"tsl",t2m+gamma*z_srf,emb%tsl,method="mean", &
                                                missing_value=real(mv,wp),fill_method="weighted")

        ! Radiative forcing, tsl_F [ J s-1 m-2] * [J-1 m2 K] == [K s-1]
        call map_scrip_field(emb%map_toemb,"tsl_F", (swn + lwn + (shf+lhf) + lhp + rco2) / tsl_fac, &
                        emb%tsl_F,method="mean",missing_value=real(mv,wp),fill_method="weighted")

        ! Wind [m s-1] 
        call map_scrip_field(emb%map_toemb,"ug",ug,emb%ug,method="mean", &
                                    missing_value=real(mv,wp),fill_method="weighted")
        call map_scrip_field(emb%map_toemb,"vg",vg,emb%vg,method="mean", &
                                    missing_value=real(mv,wp),fill_method="weighted")

        ! Calculate radiative balance over the day
        do q = 1, par%en_nstep * 5
            call solve_diffusion_advection_2D(emb%tsl,ug,vg,emb%tsl_F,emb%kappa_t,emb%tsl_bnd,emb%mask,emb%grid%dx, &
                                    emb%grid%dy,par%en_dt,k_rel=par%en_kr,solver=par%solver,step=par%step,bc="infinite")
        end do 

        ! Sea-level temperature, tsl
        call map_scrip_field(emb%map_fromemb,"tsl",emb%tsl,t2m,method="mean",missing_value=real(mv,wp), &
                                                    filt_method="gaussian",filt_par=[emb%grid%dx/2_wp,dx])
        t2m = t2m - gamma*z_srf

        return 

    end subroutine rembo_calc_en

    subroutine rembo_calc_ccw(pr,ccw,c_w,emb,par,tcw,q_r,ww,dx)

        implicit none 

        real(wp),                intent(INOUT) :: pr(:,:)
        real(wp),                intent(INOUT) :: ccw(:,:)
        real(wp),                intent(INOUT) :: c_w(:,:)
        type(diffusion_class),   intent(INOUT) :: emb 
        type(rembo_param_class), intent(IN)    :: par 
        real(wp),                intent(IN)    :: tcw(:,:)
        real(wp),                intent(IN)    :: q_r(:,:)
        real(wp),                intent(IN)    :: ww(:,:)
        real(wp),                intent(IN)    :: dx
        
        ! Local variables
        real(wp)  :: dtsl_mean, lat0, lat1, sigma 
        integer   :: q  
        real(wp), parameter  :: ccw_fac = 1.0

        !real(wp)  :: dtsl_mean    ! Average temp anomaly away from boundary temps 

        ! == Total water content ==
        call map_scrip_field(emb%map_toemb,"tcw",tcw,emb%tcw,method="mean", &
                                    missing_value=real(mv,wp),fill_method="weighted")

        ! == Relative humidity == 
        call map_scrip_field(emb%map_toemb,"q_r",q_r,emb%q_r,method="mean", &
                                    missing_value=real(mv,wp),fill_method="weighted")
        
        ! Get 2D diffusion constant
        emb%kappa_w = par%ccw_D
        emb%kappa_w = emb%kappa_w / ccw_fac 

        ! Note: removed any lat/elevation dependence of diffusion coefficient (ajr, 2015-08-31)
!         lat0 = minval(emb%grid%lat)
!         lat1 = maxval(emb%grid%lat)
!         emb%kappaw = emb%kappaw* & 
!             (1.d0 + par%en_kl*(2.d0*(emb%grid%lat-lat0)/(lat1-lat0)-1.d0))

!         emb%kappaw = emb%kappaw* (1.d0 + par%en_kz*emb%z_srf)
    
        
        ! Get vertical wind for precipitation
        call map_scrip_field(emb%map_toemb,"ww",ww,emb%ww,method="mean", &
                                    missing_value=real(mv,wp),fill_method="weighted")

        ! Get current moisture content [kg m-2]
        call map_scrip_field(emb%map_toemb,"ccw",ccw,emb%ccw,method="mean", &
                                    missing_value=real(mv,wp),fill_method="weighted")

        ! Get boundary current moisture content [kg m-2]
        ! call map_scrip_field(emb%map_toemb,"ccw_bnd",ccw_bnd,emb%ccw_bnd,method="mean", &
        !                             missing_value=real(mv,wp),fill_method="weighted")
        ! Note: ccw_bnd has been updated in `rembo_calc_iterinit`

        ! Calculate the initial lo-res precipitation rate [kg m-2 s-1]
        call map_scrip_field(emb%map_toemb,"pr",pr,emb%ccw_pr,method="mean", &
                                    missing_value=real(mv,wp),fill_method="weighted")

        ! ajr: TO DO : calculate a basic rate somehow avoiding history 


        ! Now calculate the condensation rate [kg m-2 s-1]
        emb%ccw_cw = calc_condensation(emb%tcw,emb%q_r,emb%ww,par%k_c,par%k_x,par%c%sec_day)

        ! Get moisture balance forcing [kg m-2 s-1]
        emb%ccw_F = emb%ccw_cw - emb%ccw_pr 
        !where(emb%mask==-1) emb%ccw_F = 0.d0 

        ! Calculate moisture balance to equilibrium
        do q = 1, par%ccw_nstep
            call solve_diffusion_advection_2D(emb%ccw,emb%ug,emb%vg,emb%ccw_F,emb%kappa_w,emb%ccw_bnd,emb%mask, &
                        emb%grid%dx,emb%grid%dy,par%ccw_dt,k_rel=par%ccw_kr,solver=par%solver,step=par%step,bc="infinite")

            ! Now calculate the precipitation rate on emb grid (kg m**-2 s**-1)
            ! *Cheaper than reinterpolating, but necessary to update moisture balance
            emb%ccw_pr = calc_precip(emb%ccw,emb%ww,emb%tsl-par%gamma*emb%z_srf,emb%z_srf,par%k_w,par%k_z,par%k_t) 

            ! Get moisture balance forcing [kg m-2 s-1]
            emb%ccw_F = emb%ccw_cw - emb%ccw_pr 
            where(emb%mask==1) emb%ccw_F = 0.d0 

        end do 

        ! Send cloud moisture content back to main domain pts
        call map_scrip_field(emb%map_fromemb,"ccw",emb%ccw,ccw,method="mean",missing_value=real(mv,wp), &
                                                    filt_method="gaussian",filt_par=[emb%grid%dx/2_wp,dx])
        call map_scrip_field(emb%map_fromemb,"c_w",emb%ccw_cw,c_w,method="mean",missing_value=real(mv,wp), &
                                                    filt_method="gaussian",filt_par=[emb%grid%dx/2_wp,dx])

        return 

    end subroutine rembo_calc_ccw

    subroutine rembo_gen_domain_mask(mask,zs,xx,yy,mask_domain,zs_min,radius)
        ! mask==-1: impose boundary variable
        ! mask==0: set variable to zero
        ! mask==1: solve model
        ! mask==2: solve model + relaxation term

        implicit none 
        
        integer,  intent(OUT) :: mask(:,:)
        real(wp), intent(IN) :: zs(:,:) 
        real(wp), intent(IN) :: xx(:,:)
        real(wp), intent(IN) :: yy(:,:) 
        logical,  intent(IN) :: mask_domain(:,:)
        real(wp), intent(IN) :: zs_min
        real(wp), intent(IN) :: radius

        ! Local variables
        real(wp) :: dist, mindist 
        integer :: nx, ny, i, j, i1, j1 
        !real(wp), parameter :: zs_min = 10.0   ! [m]
        
        real(wp) :: tmp(size(zs,1),size(zs,2))

        nx = size(zs,1)
        ny = size(zs,2)

        ! First assume model should be resolved everywhere
        mask = 1
        
        ! Next determine which points are ocean or low-elevation points 
        ! outside of the radius of interest.
        do j = 1, ny 
        do i = 1, nx 

            if (zs(i,j) .le. zs_min) then 
                ! Ocean point, or too-low elevation point, see if it should be modeled

                ! Loop over all land points to find minimum distance to high-elevation (land) points
                mindist = 1e8
                do i1 = 1, nx
                    do j1 = 1, ny 
                        if (zs(i1,j1) .gt. zs_min) then 
                            dist = sqrt( (xx(i1,j1)-xx(i,j))**2 + (yy(i1,j1)-yy(i,j))**2 )
                            if (dist .lt. mindist) mindist = dist  
                        end if 
                    end do 
                end do
            
                ! Ensure ocean points above minimum distance to land are solved, but
                ! with relaxation applied
                if (mindist .gt. radius) mask(i,j) = 2

            end if 

        end do 
        end do 

        ! Next ensure we do not solve model outside of the region mask 
        where(.not. mask_domain) mask = -1 

        ! Next remove any problematic 'island' points
        call remove_islands(mask)

        ! Finally, ensure the borders are set to boundary conditions
        mask(1:2,:)     = -1
        mask(nx-1:nx,:) = -1
        mask(:,1:2)     = -1
        mask(:,ny-1:ny) = -1

        return

    end subroutine rembo_gen_domain_mask 

    subroutine remove_islands(mask)

        implicit none 

        integer :: mask(:,:) 

        ! Local variables
        integer :: nx, ny, i, j
        real(wp) :: tmp(size(mask,1),size(mask,2))

        nx = size(mask,1)
        ny = size(mask,2)

        ! Remove islands
        tmp = mask

        do j = 2, ny-1
        do i = 2, nx-1 
            
            if (tmp(i,j)   .eq. -1 .and.  &
                tmp(i-1,j) .eq. 1 .and. tmp(i+1,j) .eq. 1 .and. &
                tmp(i,j-1) .eq. 1 .and. tmp(i,j+1) .eq. 1) then 

                mask(i,j) = 1

            else if (tmp(i,j)   .eq. 1 .and.  &
                     tmp(i-1,j) .eq. -1 .and. tmp(i+1,j) .eq. -1 .and. &
                     tmp(i,j-1) .eq. -1 .and. tmp(i,j+1) .eq. -1) then 
 
                mask(i,j) = -1

            end if 

        end do 
        end do 

        return 

    end subroutine remove_islands

end module rembo_atm
