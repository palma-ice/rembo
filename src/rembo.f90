module rembo
    ! Wrapper to hold all modules needed for librembo 

    use nml 
    use ncio 
    
    use rembo_defs 
    use rembo_atm 

    implicit none 

    public :: rembo_update
    public :: rembo_init 
    public :: rembo_end 

contains 


    subroutine rembo_update()
        ! Calculate atmosphere for each month of the year 

        implicit none 

        

        return 

    end subroutine rembo_update 
    
    subroutine rembo_init(dom,par_path,domain,grid,year)

        use solvers

        implicit none
        
        type(rembo_class), intent(INOUT) :: dom
        character(len=*),  intent(IN)    :: par_path  
        character(len=*),  intent(IN)    :: domain
        integer,           intent(IN)    :: year 
        type(grid_class),  intent(IN)    :: grid 

        ! Local variables         
        character(len=256) :: filename, file_boundary 
        integer            :: i
        real(wp)           :: dt_check(2) 
        character(len=64)  :: fmt1 

        ! Load the rembo parameters
        call rembo_par_load(dom%par,trim(par_path),domain)
        
        ! Initialize rembo domain and grid
        dom%grid       = grid 
        dom%par%npts   = grid%npts 
        dom%par%nx     = grid%G%nx 
        dom%par%ny     = grid%G%ny 
        dom%par%domain = trim(domain) 
        
        ! Allocate state variables
        call rembo_alloc(dom%now,dom%par%nx,dom%par%ny)
        call rembo_alloc(dom%mon,dom%par%nx,dom%par%ny)
        call rembo_alloc(dom%ann,dom%par%nx,dom%par%ny)

        ! Diffusion grid and variables
        ! ajr: TO DO 
!         call rembo_emb_init(dom%emb,grid,dx=dom%par%emb_dx)
        
        write(*,*) "rembo_init :: allocated rembo variables."

        ! Check diffusion time step consistency
        dt_check(1) = diff2D_timestep(real(dom%emb%grid%G%dx*dom%emb%grid%xy_conv,wp), &
                                      real(dom%emb%grid%G%dy*dom%emb%grid%xy_conv,wp), &
                                      dom%par%en_D)
        dt_check(2) = diff2D_timestep(real(dom%emb%grid%G%dx*dom%emb%grid%xy_conv,wp), &
                                      real(dom%emb%grid%G%dy*dom%emb%grid%xy_conv,wp), &
                                      dom%par%ccw_D)
!         dt_check(1) = diff2Dadi_timestep(dom%emb%grid%G%dx*dom%emb%grid%xy_conv, &
!                                          dom%emb%grid%G%dy*dom%emb%grid%xy_conv, &
!                                          dom%par%en_D)
!         dt_check(2) = diff2Dadi_timestep(dom%emb%grid%G%dx*dom%emb%grid%xy_conv, &
!                                          dom%emb%grid%G%dy*dom%emb%grid%xy_conv, &
!                                          dom%par%ccw_D)
        fmt1="(a,f10.1,a,f10.1,a)"
        write(*,fmt1) "Diffusion time step,   energy: ", dom%par%en_dt, " s ( max = ",dt_check(1)," s )"
        write(*,fmt1) "Diffusion time step, moisture: ", dom%par%ccw_dt," s ( max = ",dt_check(2)," s )"

        if (dom%par%en_dt .gt. dt_check(1) .or. dom%par%ccw_dt .gt. dt_check(2)) then 
            write(*,*) "Time step too big!"
            stop 
        end if 

        write(*,*) "rembo_init:: Initialization complete for domain: "// &
                   trim(dom%par%domain)

        return 

    end subroutine rembo_init

    subroutine rembo_init_state(dom)
        ! This subroutine intializes the current state variables
        ! depending on whether loading from a restart file or using default values
            
        implicit none 

        type(rembo_class) :: dom
        
        ! Initialize variables 
        if (.not. trim(dom%par%restart) .eq. "no") then 
            ! Add code to load previously stopped run
!             call rembo_restart(dom,trim(dom%par%restart))  ! ## TO DO ##

        else 
            ! Initialize variables with default values

            ! Annual values
            dom%now%co2_a = 280.d0

            ! Sub-annual values
            dom%now%t2m   = 273.15d0  
            dom%now%tcw   = 0.d0
            dom%now%ccw   = 0.d0
            dom%now%c_w   = 0.d0
            dom%now%ccw_prev = 0.d0 
            dom%now%pp    = 0.d0
            dom%now%al_p  = 0.5d0 
            dom%now%al_s  = 0.5d0
            dom%now%lwu   = 0.d0
            dom%now%swd   = 100.d0  

            dom%now%tsurf = dom%now%t2m 
            where( dom%now%tsurf .gt. 273.15d0 ) dom%now%tsurf = 273.15d0
            dom%now%swd_s = 0.d0 
            dom%now%lwd_s = 0.d0  
            dom%now%u_s   = 0.d0 
            dom%now%v_s   = 0.d0 

            ! Diffusion variables
            dom%emb%tsl    = 273.15d0  ! 0 degC
            dom%emb%tcw    = 0.d0 
            dom%emb%ug     = 0.d0 
            dom%emb%vg     = 0.d0 

            dom%emb%en     = 0.d0
            dom%emb%en_bnd = 0.d0  
            dom%emb%en_F   = 0.d0 
            dom%emb%kappa  = dom%par%en_D 

            dom%emb%ccw     = 0.03d0    ! Average winter value 
            dom%emb%ccw_bnd = 0.d0
            dom%emb%ccw_F   = 0.d0 
            dom%emb%kappaw  = dom%par%ccw_D
            
        end if 

        return 

    end subroutine rembo_init_state 

    subroutine rembo_end(dom)

        implicit none 

        type(rembo_class), intent(INOUT) :: dom

        ! Atmospheric variables
        call rembo_dealloc(dom%now)
        call rembo_dealloc(dom%mon)
        call rembo_dealloc(dom%ann)

        write(*,*) "rembo_end :: deallocated rembo variables: "//trim(dom%par%domain)
        write(*,*)

        return 

    end subroutine rembo_end 
    
    subroutine rembo_print(dom,m,d,year)

        implicit none 

        type(rembo_class) :: dom 
        integer  :: m, d, year
        real(wp) :: npts 
        character(len=256) :: fmt_head, fmt_table 

        npts = real(count(dom%now%mask == 2),wp)

        fmt_head  = "(a6,a5,  5a8,2x,  a8,  3a8,2x,  4a8)"
        fmt_table = "(i6,i5,5f8.2,2x,f8.2,3f8.2,2x,4f8.2)"

        if (m .eq. 1) &
        write(*,fmt_head) "year","day","swd","tas", &
                           "tcw","ccw","c_w","pr"
                           
        write(*,fmt_table) year, d, &
            sum(dom%now%swd,        mask=dom%now%mask==2)/npts, &     ! [W/m^2]
            sum(dom%now%t2m,        mask=dom%now%mask==2)/npts, &     ! [K]
            sum(dom%now%tcw,        mask=dom%now%mask==2)/npts, &     ! [mm]
            sum(dom%now%ccw,        mask=dom%now%mask==2)/npts, &     ! [mm]
            sum(dom%now%c_w*sec_day,mask=dom%now%mask==2)/npts, &     ! [mm/d]
            sum(dom%now%pp*sec_day, mask=dom%now%mask==2)/npts        ! [mm/d]

        return 

    end subroutine rembo_print 
    
    subroutine rembo_par_load(par,filename,domain)

        type(rembo_param_class), intent(INOUT) :: par 
        character(len=*),        intent(IN)    :: filename 
        character(len=*),        intent(IN)    :: domain 

        par%domain = trim(domain)
        call nml_read(filename,"rembo_params","restart",    par%restart    )

        call nml_read(filename,"rembo_params","H_e",  par%H_e   )
        call nml_read(filename,"rembo_params","emb_dx",  par%emb_dx  )
        call nml_read(filename,"rembo_params","rembo1",  par%rembo1  )
        call nml_read(filename,"rembo_params","dist_rel",  par%dist_rel  )
        call nml_read(filename,"rembo_params","en_dt",  par%en_dt  )
        call nml_read(filename,"rembo_params","en_D",  par%en_D  )
        call nml_read(filename,"rembo_params","en_D_win",  par%en_D_win  )
        call nml_read(filename,"rembo_params","en_D_sum",  par%en_D_sum  )
        call nml_read(filename,"rembo_params","en_kr",  par%en_kr  )
        call nml_read(filename,"rembo_params","en_kz",  par%en_kz  )
        call nml_read(filename,"rembo_params","en_kl",  par%en_kl  )
        call nml_read(filename,"rembo_params","en_kdT",  par%en_kdT  )
        call nml_read(filename,"rembo_params","en_Ha",  par%en_Ha  )
        call nml_read(filename,"rembo_params","ccw_dt",  par%ccw_dt  )
        call nml_read(filename,"rembo_params","ccw_D",  par%ccw_D  )
        call nml_read(filename,"rembo_params","ccw_kr",  par%ccw_kr  )
        call nml_read(filename,"rembo_params","k_c",  par%k_c  )
        call nml_read(filename,"rembo_params","k_p",  par%k_p  )
        call nml_read(filename,"rembo_params","k_z",  par%k_z  )
        call nml_read(filename,"rembo_params","k_x",  par%k_x  )
        call nml_read(filename,"rembo_params","e0",  par%e0  )
        call nml_read(filename,"rembo_params","c1",  par%c1  )
        call nml_read(filename,"rembo_params","k_p_now",  par%k_p_now  )
        call nml_read(filename,"rembo_params","k_t",  par%k_t  )
        call nml_read(filename,"rembo_params","k_w",  par%k_w  )
        call nml_read(filename,"rembo_params","ps_a",  par%ps_a  )
        call nml_read(filename,"rembo_params","ps_b",  par%ps_b  )
        call nml_read(filename,"rembo_params","f_k",  par%f_k  )
        call nml_read(filename,"rembo_params","nk1",  par%nk1  )
        call nml_read(filename,"rembo_params","nk2",  par%nk2  )
        call nml_read(filename,"rembo_params","nk3",  par%nk3  )
        call nml_read(filename,"rembo_params","teff_sigma",  par%teff_sigma  )
        call nml_read(filename,"rembo_params","gamma",  par%gamma  )
        call nml_read(filename,"rembo_params","gamma2",  par%gamma2  )
        call nml_read(filename,"rembo_params","S0",  par%S0  )
        call nml_read(filename,"rembo_params","Lw",  par%Lw  )
        call nml_read(filename,"rembo_params","Lm",  par%Lm  )
        call nml_read(filename,"rembo_params","Ls",  par%Ls  )
        call nml_read(filename,"rembo_params","rho_w",  par%rho_w  )
        call nml_read(filename,"rembo_params","cp",  par%cp  )
        call nml_read(filename,"rembo_params","ci",  par%ci  )
        call nml_read(filename,"rembo_params","alp_a",  par%alp_a  )
        call nml_read(filename,"rembo_params","alp_b",  par%alp_b  )
        
        call nml_read(filename,"rembo_params","alp_c",  par%alp_c  )
        call nml_read(filename,"rembo_params","lwu_a",  par%lwu_a  )
        call nml_read(filename,"rembo_params","lwu_b",  par%lwu_b  )
        call nml_read(filename,"rembo_params","lwu_c",  par%lwu_c  )
        call nml_read(filename,"rembo_params","swds_a",  par%swds_a  )
        call nml_read(filename,"rembo_params","swds_b",  par%swds_b  )
        call nml_read(filename,"rembo_params","swds_c",  par%swds_c  )
        call nml_read(filename,"rembo_params","lwds_a",  par%lwds_a  )
        call nml_read(filename,"rembo_params","lwds_b",  par%lwds_b  )
        call nml_read(filename,"rembo_params","lwds_c",  par%lwds_c  )
        call nml_read(filename,"rembo_params","lwds_d",  par%lwds_d  )
        
        call nml_read(filename,"rembo_params","shfs_Cm",  par%shfs_Cm)
        call nml_read(filename,"rembo_params","shfs_p",  par%shfs_p  )

        call nml_read(filename,"rembo1_params","ta",  par%r1_ta  )
        call nml_read(filename,"rembo1_params","tb",  par%r1_tb  )

        ! How many time steps in 1 day, aprx?
        par%en_nstep  = floor(sec_day / par%en_dt)
        par%ccw_nstep = floor(sec_day / par%ccw_dt)
        
        return

    end subroutine rembo_par_load


    subroutine rembo_alloc(now,nx,ny)

        implicit none 

        type(rembo_state_class), intent(INOUT) :: now 
        integer, intent(IN) :: nx, ny 

        ! Make sure this object is fully deallocated first
        call rembo_dealloc(now)

        ! Allocate variables 

        return 

    end subroutine rembo_alloc 
    
    subroutine rembo_dealloc(now)

        implicit none 

        type(rembo_state_class), intent(INOUT) :: now

        ! Deallocate variables 

        return 

    end subroutine rembo_dealloc

end module rembo 