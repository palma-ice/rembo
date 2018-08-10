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
        call rembo_physics_par_load(dom%par,trim(par_path))

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

    subroutine rembo_end()

        implicit none 


        return 

    end subroutine rembo_end 
    
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