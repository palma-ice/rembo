module solvers

    use, intrinsic :: iso_fortran_env, only : input_unit, output_unit, error_unit
    use precision, only : wp
    use solver_linear 

    implicit none
    
    private
    public :: solve_diffusion_2D
    public :: solve_advection_2D
    public :: solve_diffusion_advection_2D
    public :: adv2D_timestep
    public :: diff2D_timestep
    public :: diff2Dadi_timestep
    public :: solve2D
    public :: solve_diff_2D_adi

contains

    subroutine solve_diffusion_2D(uu,F,kappa,ubnd,mask,dx,dy,dt,k_rel,solver,step, &
                                                                    bc,bc1,bc2,bc3,bc4)

        implicit none 

        real(wp), intent(INOUT) :: uu(:,:)
        real(wp), intent(IN)    :: F(:,:)
        real(wp), intent(IN)    :: kappa(:,:)
        real(wp), intent(IN)    :: ubnd(:,:)
        integer,  intent(IN)    :: mask(:,:)
        real(wp), intent(IN)    :: dx
        real(wp), intent(IN)    :: dy
        real(wp), intent(IN)    :: dt
        real(wp), intent(IN)    :: k_rel
        character(len=*), intent(IN) :: solver
        character(len=*), intent(IN) :: step
        character(len=*), intent(IN), optional :: bc
        character(len=*), intent(IN), optional :: bc1
        character(len=*), intent(IN), optional :: bc2
        character(len=*), intent(IN), optional :: bc3
        character(len=*), intent(IN), optional :: bc4

        ! Local variables
        integer :: nx, ny 
        real(wp), allocatable :: zeros(:,:)

        nx = size(uu,1)
        ny = size(uu,2)
        
        allocate(zeros(nx,ny))
        zeros = 0.0 

        ! Call more general diffusion-advection routine with no velocity
        call solve_diffusion_advection_2D(uu,zeros,zeros,F,kappa,ubnd,mask,dx,dy,dt,k_rel,solver,step, &
                                                                                    bc,bc1,bc2,bc3,bc4)

        return

    end subroutine solve_diffusion_2D

    subroutine solve_advection_2D(uu,v_x,v_y,F,ubnd,mask,dx,dy,dt,k_rel,solver,step, &
                                                                        bc,bc1,bc2,bc3,bc4)

        implicit none 

        real(wp), intent(INOUT) :: uu(:,:)
        real(wp), intent(IN)    :: v_x(:,:)
        real(wp), intent(IN)    :: v_y(:,:)
        real(wp), intent(IN)    :: F(:,:)
        real(wp), intent(IN)    :: ubnd(:,:)
        integer,  intent(IN)    :: mask(:,:)
        real(wp), intent(IN)    :: dx
        real(wp), intent(IN)    :: dy
        real(wp), intent(IN)    :: dt
        real(wp), intent(IN)    :: k_rel
        character(len=*), intent(IN) :: solver
        character(len=*), intent(IN) :: step
        character(len=*), intent(IN), optional :: bc
        character(len=*), intent(IN), optional :: bc1
        character(len=*), intent(IN), optional :: bc2
        character(len=*), intent(IN), optional :: bc3
        character(len=*), intent(IN), optional :: bc4

        ! Local variables
        integer :: nx, ny 
        real(wp), allocatable :: zeros(:,:)

        nx = size(uu,1)
        ny = size(uu,2)
        
        allocate(zeros(nx,ny))
        zeros = 0.0 

        ! Call more general diffusion-advection routine with no diffusion constant
        call solve_diffusion_advection_2D(uu,v_x,v_y,F,zeros,ubnd,mask,dx,dy,dt,k_rel,solver,step, &
                                                                                    bc,bc1,bc2,bc3,bc4)

        return

    end subroutine solve_advection_2D

    subroutine solve_diffusion_advection_2D(uu,v_x,v_y,F,kappa,ubnd,mask,dx,dy,dt,k_rel,solver,step, &
                                                                                    bc,bc1,bc2,bc3,bc4)

        implicit none 

        real(wp), intent(INOUT) :: uu(:,:)
        real(wp), intent(IN)    :: v_x(:,:)
        real(wp), intent(IN)    :: v_y(:,:)
        real(wp), intent(IN)    :: F(:,:)
        real(wp), intent(IN)    :: kappa(:,:)
        real(wp), intent(IN)    :: ubnd(:,:)
        integer,  intent(IN)    :: mask(:,:)
        real(wp), intent(IN)    :: dx
        real(wp), intent(IN)    :: dy
        real(wp), intent(IN)    :: dt
        real(wp), intent(IN)    :: k_rel
        character(len=*), intent(IN) :: solver
        character(len=*), intent(IN) :: step
        character(len=*), intent(IN), optional :: bc
        character(len=*), intent(IN), optional :: bc1
        character(len=*), intent(IN), optional :: bc2
        character(len=*), intent(IN), optional :: bc3
        character(len=*), intent(IN), optional :: bc4

        ! Local variables
        integer :: nx, ny

        real(wp), allocatable :: f0(:,:)
        real(wp), allocatable :: f1(:,:)
        real(wp), allocatable :: f2(:,:)
        real(wp), allocatable :: f3(:,:)
        real(wp), allocatable :: u0(:,:)
        real(wp), allocatable :: u1(:,:)
        real(wp), allocatable :: u2(:,:)
        real(wp), allocatable :: u3(:,:)

        character(len=56) :: bcs(4)

        nx = size(uu,1)
        ny = size(uu,2)

        ! Set boundary conditions
        call set_boundary_conditions(bcs,bc,bc1,bc2,bc3,bc4)
        
        ! Perform timestepping

        select case(trim(step))

            case("fe")
                ! Forward-euler timestepping

                allocate(f0(nx,ny))

                ! Calculate derivative
                call calc_tendency(f0,uu,v_x,v_y,F,kappa,ubnd,mask,dx,dy,dt,k_rel,solver,bcs)

                ! Update uu
                uu = uu + dt*f0

            case("rk4")
                ! Runge-kutta fourth order timestepping
                
                ! Get four values of the derivative between t0 and t0+dt
                
                ! call f ( t0, u0, f0 )

                ! t1 = t0 + dt / 2.0
                ! u1 = u0 + dt * f0 / 2.0
                ! call f ( t1, u1, f1 )

                ! t2 = t0 + dt / 2.0
                ! u2 = u0 + dt * f1 / 2.0
                ! call f ( t2, u2, f2 )

                ! t3 = t0 + dt
                ! u3 = u0 + dt * f2
                ! call f ( t3, u3, f3 )

                ! Combine them to estimate the solution U at time T1
                ! u = u0 + dt * ( f0 + 2.0*f1 + 2.0*f2 + f3 ) / 6.0

                allocate(f0(nx,ny))
                allocate(f1(nx,ny))
                allocate(f2(nx,ny))
                allocate(f3(nx,ny))
                
                allocate(u0(nx,ny))
                allocate(u1(nx,ny))
                allocate(u2(nx,ny))
                allocate(u3(nx,ny))

                u0 = uu 

                call calc_tendency(f0,u0,v_x,v_y,F,kappa,ubnd,mask,dx,dy,dt,k_rel,solver,bcs)

                u1 = u0 + (dt/2.0)*f0

                call calc_tendency(f1,u1,v_x,v_y,F,kappa,ubnd,mask,dx,dy,dt,k_rel,solver,bcs)

                u2 = u0 + (dt/2.0)*f1
                
                call calc_tendency(f2,u2,v_x,v_y,F,kappa,ubnd,mask,dx,dy,dt,k_rel,solver,bcs)

                u3 = u0 + dt*f2

                call calc_tendency(f3,u3,v_x,v_y,F,kappa,ubnd,mask,dx,dy,dt,k_rel,solver,bcs)

                ! Combine them to estimate the solution U at time t+dt
                uu = u0 + dt * ( f0 + 2.0*f1 + 2.0*f2 + f3 ) / 6.0

            case DEFAULT

                write(*,*) "solve_diffusion_advection_2D:: timestepping method not recongized."
                write(*,*) "step = ", trim(step)
                stop 
                
        end select
        
        ! Impose hard boundary conditions
        !where (mask .eq. -1) uu = ubnd 

        return

    end subroutine solve_diffusion_advection_2D

    subroutine calc_tendency(dudt,uu,v_x,v_y,F,kappa,ubnd,mask,dx,dy,dt,k_rel,solver,bcs)
        ! Calculate advection diffusion equation tendency
        ! du/dt = kappa * grad^2 T - u<dot>T + F + F_relax

        implicit none 

        real(wp), intent(INOUT) :: dudt(:,:)
        real(wp), intent(IN)    :: uu(:,:)
        real(wp), intent(IN)    :: v_x(:,:)
        real(wp), intent(IN)    :: v_y(:,:)
        real(wp), intent(IN)    :: F(:,:)
        real(wp), intent(IN)    :: kappa(:,:)
        real(wp), intent(IN)    :: ubnd(:,:)
        integer,  intent(IN)    :: mask(:,:)
        real(wp), intent(IN)    :: dx
        real(wp), intent(IN)    :: dy
        real(wp), intent(IN)    :: dt
        real(wp), intent(IN)    :: k_rel
        character(len=*), intent(IN) :: solver
        character(len=*), intent(IN) :: bcs(4)

        ! Local variables
        integer  :: i, j, nx, ny 
        real(wp) :: diff 
        real(wp), allocatable :: dudt_adv(:,:)
        real(wp), allocatable :: dudt_diff(:,:)
        real(wp), allocatable :: dudt_relax(:,:)

        nx = size(dudt,1)
        ny = size(dudt,2) 
        
        allocate(dudt_diff(nx,ny))
        allocate(dudt_adv(nx,ny))
        allocate(dudt_relax(nx,ny))

        ! Get relaxation rate at boundaries 
        dudt_relax = 0.0
        where (mask .eq. -2)
            ! Imposed relaxation rate (k_rel is restoring rate fraction per second, positive value)
            dudt_relax = -k_rel*(uu-ubnd)
        end where 
        
        select case(trim(solver))

            case("none")
                ! No solver right now, set dudt to zero

                dudt = 0.0 

            case("expl-diff")
                ! Explicit diffusion

                ! Get diffusive tendency, i.e. kappa*(Laplacian operator) (du/dt)
                call calc_tendency_diffuse2D(dudt_diff,uu,kappa,mask,dx,dy,bcs)

                ! No advection considered
                dudt_adv = 0.0
        
                ! Get total tendency for this timestep
                dudt = dudt_diff + dudt_adv + F + dudt_relax

            case("expl-adv")
                ! Explicit advection 

                ! No diffusion considered
                dudt_diff = 0.0 

                ! Get advective tendency (du/dt)
                call calc_tendency_advec2D_upwind(dudt_adv,uu,v_x,v_y,mask,dx,dy,bcs)

                ! Get total tendency for this timestep
                dudt = dudt_diff + dudt_adv + F + dudt_relax
                
            case("expl")
                ! Explicit diffusion and advection 

                ! Get diffusive tendency, i.e. kappa*(Laplacian operator) (du/dt)
                call calc_tendency_diffuse2D(dudt_diff,uu,kappa,mask,dx,dy,bcs)

                ! Get advective tendency (du/dt)
                call calc_tendency_advec2D_upwind(dudt_adv,uu,v_x,v_y,mask,dx,dy,bcs)

                ! Get total tendency for this timestep
                dudt = dudt_diff + dudt_adv + F + dudt_relax
                
            case("impl")
                ! Implicit diffusion and advection 

                ! Get total tendency (du/dt): diffusive-advective + F + relaxation
                call calc_tendency_diffadvec2D_upwind_impl(dudt,uu,v_x,v_y,kappa, &
                                                           F+dudt_relax,ubnd,mask,dx,dy,dt,bcs)
                
            case ("adi-impl")
                ! ADI diffusion, implicit advection

                ! TO DO 
                write(*,*) "calc_tendency:: method not yet implemented."
                write(*,*) "solver = ", trim(solver)
                stop 

            case DEFAULT

                write(*,*) "calc_tendency:: method not recongized."
                write(*,*) "solver = ", trim(solver)
                stop 
                
        end select

        return

    end subroutine calc_tendency 

    subroutine calc_tendency_diffuse2D(dudt,uu,kappa,mask,dx,dy,bcs)
        ! Explicit calculation of kappa * 2D Laplace: kappa*d2u/dxdy
        ! If input units are [u], returns [u/m2]
    
        implicit none 

        real(wp), intent(OUT) :: dudt(:,:)
        real(wp), intent(IN)  :: uu(:,:)
        real(wp), intent(IN)  :: kappa(:,:)
        integer,  intent(IN)  :: mask(:,:)
        real(wp), intent(IN)  :: dx
        real(wp), intent(IN)  :: dy
        character(len=*), intent(IN) :: bcs(4)

        ! Local variables 
        integer :: i, j, nx, ny
        integer :: im1, ip1, jm1, jp1 
        real(wp) :: inv_dx2, inv_dy2
        real(wp) :: du

        nx = size(dudt,1)
        ny = size(dudt,2) 

        inv_dx2 = 1.d0 / (dx*dx)
        inv_dy2 = 1.d0 / (dy*dy)

        dudt = 0.d0 

        do j = 1,ny
        do i = 1,nx
            
            if (mask(i,j) .ne. -1) then

                ! Get neighbor indices
                call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,bcs)

                ! Laplacian, five-point stencil finite difference method (du/(dx*dy))
                du =  inv_dx2*(uu(im1,j)-2.d0*uu(i,j)+uu(ip1,j)) &
                    + inv_dy2*(uu(i,jm1)-2.d0*uu(i,j)+uu(i,jp1))
                
                ! Get tendency (du/dt)
                dudt(i,j) = kappa(i,j)*du

            else
                ! Impose no change

                dudt(i,j) = 0.0

            end if 

        end do 
        end do

        return

    end subroutine calc_tendency_diffuse2D 

    subroutine calc_tendency_advec2D_upwind(dudt,uu,vx,vy,mask,dx,dy,bcs)
        ! Second-order upwind advection scheme
        ! If input units are [u], returns [u/s]
        
        implicit none 

        real(wp), intent(OUT) :: dudt(:,:)
        real(wp), intent(IN)  :: uu(:,:)
        real(wp), intent(IN)  :: vx(:,:)
        real(wp), intent(IN)  :: vy(:,:)
        integer,  intent(IN)  :: mask(:,:)
        real(wp), intent(IN)  :: dx
        real(wp), intent(IN)  :: dy 
        character(len=*), intent(IN) :: bcs(4)

        ! Local variables
        integer :: i, j, nx, ny 

        nx = size(dudt,1)
        ny = size(dudt,2)

        dudt = 0.0 

        do j = 1, ny
        do i = 1, nx 

            if (mask(i,j) .ne. -1) then
                call calc_advec_horizontal_point(dudt(i,j),uu,vx,vy,dx,i,j,bcs)
            else
                dudt(i,j) = 0.0
            end if

        end do 
        end do


        return

    end subroutine calc_tendency_advec2D_upwind

    subroutine calc_advec_horizontal_point(dvardt,var,ux,uy,dx,i,j,bcs)
        ! Newly implemented advection algorithms (ajr)
        ! Output: [[var] [time]-1]

        ! e.g., [var] = [K]; [time] = [a-1]; [m-1] * [m a-1] * [K] = [K a-1]

        implicit none

        real(wp), intent(OUT) :: dvardt
        real(wp), intent(IN)  :: var(:,:)       ! nx,ny  Enth, T, age, etc...
        real(wp), intent(IN)  :: ux(:,:)        ! nx,ny
        real(wp), intent(IN)  :: uy(:,:)        ! nx,ny 
        real(wp), intent(IN)  :: dx  
        integer,    intent(IN)  :: i, j 
        character(len=*), intent(IN) :: bcs(4) 

        ! Local variables 
        integer :: k, nx, ny 
        integer :: im1, ip1, jm1, jp1
        real(wp) :: ux_aa, uy_aa 
        real(wp) :: dx_inv, dx_inv2
        real(wp) :: advecx, advecy, advec_rev 

        ! Define some constants 
        dx_inv  = 1.0_wp / dx 
        dx_inv2 = 1.0_wp / (2.0_wp*dx)

        nx  = size(var,1)
        ny  = size(var,2)

        advecx  = 0.0 
        advecy  = 0.0 
        dvardt  = 0.0 

        ! Get neighbor indices
        call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,bcs)

        ! Estimate direction of current flow into cell (x and y), centered in grid point
        ux_aa = 0.5_wp*(ux(i,j)+ux(im1,j))
        uy_aa = 0.5_wp*(uy(i,j)+uy(i,jm1))
        
        ! Explicit form (to test different order approximations)
        if (ux(im1,j) .gt. 0.0_wp .and. ux(i,j) .lt. 0.0_wp .and. i .ge. 3 .and. i .le. nx-2) then 
            ! Convergent flow - take the mean 

            ! 2nd order
            !advecx    = dx_inv2 * ux(i-1,j)*(-(4.0*var(i-1,j)-var(i-2,j)-3.0*var(i,j)))
            !advec_rev = dx_inv2 * ux(i,j)*((4.0*var(i+1,j)-var(i+2,j)-3.0*var(i,j)))

            ! 1st order
            advecx    = dx_inv * ux(im1,j)*(-(var(im1,j)-var(i,j)))
            advec_rev = dx_inv * ux(i,j)*((var(ip1,j)-var(i,j)))
            
            advecx    = 0.5_wp * (advecx + advec_rev) 

        else if (ux_aa .gt. 0.0 .and. i .ge. 3) then  
            ! Flow to the right - inner points

            ! 2nd order
            !advecx = dx_inv2 * ux(i-1,j)*(-(4.0*var(i-1,j)-var(i-2,j)-3.0*var(i,j)))

            ! 1st order
            advecx = dx_inv * ux(im1,j)*(-(var(im1,j)-var(i,j)))
            
        else if (ux_aa .gt. 0.0 .and. i .eq. 2) then  
            ! Flow to the right - border points

            ! 1st order
            advecx = dx_inv * ux(im1,j)*(-(var(im1,j)-var(i,j)))
            
        else if (ux_aa .lt. 0.0 .and. i .le. nx-2) then 
            ! Flow to the left

            ! 2nd order
            !advecx = dx_inv2 * ux(i,j)*((4.0*var(i+1,j)-var(i+2,j)-3.0*var(i,j)))

            ! 1st order 
            advecx = dx_inv * ux(i,j)*((var(ip1,j)-var(i,j)))
            
        else if (ux_aa .lt. 0.0 .and. i .eq. nx-1) then 
            ! Flow to the left

            ! 1st order 
            advecx = dx_inv * ux(i,j)*((var(ip1,j)-var(i,j)))
            
        else 
            ! No flow or divergent 

            advecx = 0.0

        end if 

        if (uy(i,j-1) .gt. 0.0_wp .and. uy(i,j) .lt. 0.0_wp .and. j .ge. 3 .and. j .le. ny-2) then 
            ! Convergent flow - take the mean 

            ! 2nd order
            !advecy    = dx_inv2 * uy(i,j-1)*(-(4.0*var(i,j-1)-var(i,j-2)-3.0*var(i,j)))
            !advec_rev = dx_inv2 * uy(i,j)*((4.0*var(i,j+1)-var(i,j+2)-3.0*var(i,j)))
            
            ! 1st order
            advecy    = dx_inv * uy(i,jm1)*(-(var(i,jm1)-var(i,j)))
            advec_rev = dx_inv * uy(i,j)*((var(i,jp1)-var(i,j)))
            
            advecy    = 0.5_wp * (advecy + advec_rev) 

        else if (uy_aa .gt. 0.0 .and. j .ge. 3) then   
            ! Flow to the right  - inner points

            ! 2nd order
            !advecy = dx_inv2 * uy(i,j-1)*(-(4.0*var(i,j-1)-var(i,j-2)-3.0*var(i,j)))

            ! 1st order
            advecy = dx_inv * uy(i,jm1)*(-(var(i,jm1)-var(i,j)))
            
        else if (uy_aa .gt. 0.0 .and. j .eq. 2) then   
            ! Flow to the right - border points

            ! 1st order
            advecy = dx_inv * uy(i,jm1)*(-(var(i,jm1)-var(i,j)))
            
        else if (uy_aa .lt. 0.0 .and. j .le. ny-2) then 
            ! Flow to the left

            ! 2nd order
            !advecy = dx_inv2 * uy(i,j)*((4.0*var(i,j+1)-var(i,j+2)-3.0*var(i,j)))
            
            ! 1st order
            advecy = dx_inv * uy(i,j)*((var(i,jp1)-var(i,j)))
            
        else if (uy_aa .lt. 0.0 .and. j .eq. ny-1) then 
            ! Flow to the left

            ! 1st order
            advecy = dx_inv * uy(i,j)*((var(i,jp1)-var(i,j)))
                
        else
            ! No flow 
            advecy = 0.0 

        end if 
        
        ! Combine advection terms for total contribution 
        dvardt = -(advecx+advecy)

        return 

    end subroutine calc_advec_horizontal_point
    
    ! ==== LIS-based solving routines ====

    subroutine calc_tendency_diffadvec2D_upwind_impl(dHdt,H,ux,uy,kappa,F,Hbnd,mask,dx,dy,dt,bcs)
        ! General routine to apply 2D advection equation to variable `var` 
        ! with source term `F`. Various solvers are possible

        real(wp),       intent(OUT)   :: dHdt(:,:)              ! [dHdt] Variable rate of change
        real(wp),       intent(IN)    :: H(:,:)                 ! [H]    Variable to be updated
        real(wp),       intent(IN)    :: ux(:,:)                ! [m/[time]] 2D velocity, x-direction (ac-nodes)
        real(wp),       intent(IN)    :: uy(:,:)                ! [m/[time]] 2D velocity, y-direction (ac-nodes)
        real(wp),       intent(IN)    :: kappa(:,:)             ! [m2/[time]] Diffusion constant
        real(wp),       intent(IN)    :: F(:,:)                 ! [d[H]/d[time]] Source term for variable
        real(wp),       intent(IN)    :: Hbnd(:,:)              ! [H]    Variable boundary values
        integer,        intent(IN)    :: mask(:,:)              ! Mask
        real(wp),       intent(IN)    :: dx                     ! [m]   Horizontal resolution, x-direction
        real(wp),       intent(IN)    :: dy                     ! [m]   Horizontal resolution, y-direction
        real(wp),       intent(IN)    :: dt                     ! [[time]]   Timestep 
        character(len=*), intent(IN)  :: bcs(4)

        ! Local variables
        integer :: nx, ny, i, j 
        type(linear_solver_class) :: lgs
        character(len=256)        :: lis_opt 
        real(wp), allocatable :: H_new(:,:) 

        nx = size(H,1)
        ny = size(H,2)

        allocate(H_new(nx,ny))

        ! Initialize linear solver variables
        call linear_solver_init(lgs,nx,ny,nvar=1,n_terms=5)

        ! Populate advection matrices Ax=b
        call linear_solver_matrix_diffusion_advection_csr_2D(lgs,H,ux,uy,kappa,F,Hbnd,mask,dx,dy,dt,bcs)
        
        ! Solve linear equation
        lis_opt = "-i bicg -p ilu -maxiter 1000 -tol 1.0e-12 -initx_zeros false"
        !lis_opt = "-i minres -p jacobi -maxiter 1000 -tol 1.0e-12 -initx_zeros false"
        call linear_solver_matrix_solve(lgs,lis_opt)
        
        !call linear_solver_print_summary(lgs,io_unit_err)

        ! Store new solution
        call linear_solver_save_variable(H_new,lgs)

        ! Determine rate of change 
        dHdt = (H_new - H) / dt 

        return 

    end subroutine calc_tendency_diffadvec2D_upwind_impl

    subroutine linear_solver_save_variable(H,lgs)
        ! Extract ice thickness solution from lgs object. 

        implicit none 

        real(wp), intent(OUT) :: H(:,:)                 ! [X]
        type(linear_solver_class), intent(IN) :: lgs

        ! Local variables 
        integer :: i, j, nr 

        do nr = 1, lgs%nmax

            i = lgs%n2i(nr)
            j = lgs%n2j(nr)

            H(i,j) = lgs%x_value(nr)

        end do

        return

    end subroutine linear_solver_save_variable

    subroutine linear_solver_matrix_diffusion_advection_csr_2D(lgs,H,ux,uy,kappa,F,Hbnd,mask,dx,dy,dt,bcs)
        ! Define sparse matrices A*x=b in format 'compressed sparse row' (csr)
        ! for 2D advection equations with velocity components
        ! ux and uy defined on ac-nodes (right and top borders of i,j grid cell)
        ! and variable to be advected H defined on aa nodes.
        ! Store sparse matrices in linear_solver_class object 'lgs' for later use.
        
        ! Boundary conditions (bcs) counterclockwise unit circle 
        ! 1: x, right-border
        ! 2: y, upper-border 
        ! 3: x, left--border 
        ! 4: y, lower-border 
        
        implicit none 

        type(linear_solver_class), intent(INOUT) :: lgs
        real(wp), intent(IN)      :: H(:,:)         ! [X] Variable of interest (aa nodes)
        real(wp), intent(IN)      :: ux(:,:)        ! [m a-1] Horizontal velocity x-direction (ac nodes)
        real(wp), intent(IN)      :: uy(:,:)        ! [m a-1] Horizontal velocity y-direction (ac nodes)
        real(wp), intent(IN)      :: kappa(:,:)     ! [m2 a-1] Diffusion constant
        real(wp), intent(IN)      :: F(:,:)         ! [m a-1] Net source/sink terms (aa nodes)
        real(wp), intent(IN)      :: Hbnd(:,:)      ! [X] Boundary variable (aa nodes)
        integer,  intent(IN)      :: mask(:,:)      ! Advection mask
        real(wp), intent(IN)      :: dx             ! [m] Horizontal step x-direction
        real(wp), intent(IN)      :: dy             ! [m] Horizontal step y-direction 
        real(wp), intent(IN)      :: dt             ! [a] Time step 
        character(len=56), intent(IN) :: bcs(4)     ! Boundary conditions

        ! Local variables  
        integer  :: i, j, k, nx, ny
        integer  :: im1, ip1, jm1, jp1 
        integer  :: n, nr, nc
        real(wp) :: dt_darea
        real(wp) :: dt_dx2
        real(wp) :: dt_dy2
        real(wp) :: alpha
        real(wp) :: beta
        
        real(wp), allocatable  :: ux_1(:,:), ux_2(:,:)
        real(wp), allocatable  :: uy_1(:,:), uy_2(:,:)
        real(wp), allocatable  :: Hx_1(:,:), Hx_2(:,:)
        real(wp), allocatable  :: Hy_1(:,:), Hy_2(:,:)

        real(wp), parameter :: WOVI = 1.0       ! Weighting parameter for the over-implicit scheme 
                                                ! WOVI=0: Explicit scheme
                                                ! WOVI=1: Implicit scheme
                                                ! WOVI=1.5: Over-implicit scheme

        nx = size(H,1)
        ny = size(H,2) 

        dt_darea = dt/(dx*dy)
        dt_dx2   = dt/(dx*dx)
        dt_dy2   = dt/(dy*dy)

                ! Safety check for initialization
        if (.not. allocated(lgs%x_value)) then 
            ! Object 'lgs' has not been initialized yet.

            write(error_unit,*) "Error: linear_solver_matrix_boundaries_csr_2D:: lgs object not &
            &initialized yet. First call linear_solver_init(lgs,nx,ny,nvar,n_terms)."
            stop "Error: linear_solver_matrix_boundaries_csr_2D."

        end if

        ! Allocate local variables
        allocate(ux_1(nx,ny))
        allocate(ux_2(nx,ny))
        allocate(uy_1(nx,ny))
        allocate(uy_2(nx,ny))
        
        allocate(Hx_1(nx,ny))
        allocate(Hx_2(nx,ny))
        allocate(Hy_1(nx,ny))
        allocate(Hy_2(nx,ny))
        
        ! ====================================================================================
        ! Step 1: populate variables representing velocity components (ux,uy) at each
        ! cell face (left,right,bottom,top) and upstream variable (Hx,Hy) at each cell face. 

        ux_1 = 0.0
        ux_2 = 0.0
        uy_1 = 0.0
        uy_2 = 0.0
        
        Hx_1 = 0.0
        Hx_2 = 0.0
        Hy_1 = 0.0
        Hy_2 = 0.0

        do j = 1, ny
        do i = 1, nx 

            ! Get neighbor indices
            call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,bcs)

            ! Get velocity components on each cell face 
            ux_1(i,j) = ux(im1,j)
            ux_2(i,j) = ux(i,j)
            uy_1(i,j) = uy(i,jm1)
            uy_2(i,j) = uy(i,j)

            if (ux_1(i,j) >= 0.0) then
                Hx_1(i,j) = H(im1,j)
            else
                Hx_1(i,j) = H(i,j)
            end if

            if (ux_2(i,j) >= 0.0) then
                Hx_2(i,j) = H(i,j)
            else
                Hx_2(i,j) = H(ip1,j)
            end if

            if (uy_1(i,j) >= 0.0) then
                Hy_1(i,j) = H(i,jm1)
            else
                Hy_1(i,j) = H(i,j)
            end if

            if (uy_2(i,j) >= 0.0) then
                Hy_2(i,j) = H(i,j)
            else
                Hy_2(i,j) = H(i,jp1)
            end if

        end do
        end do

        ! ====================================================================================
        ! Step 2: Assembly of the system of linear equations
        !             (matrix storage: compressed sparse row CSR)

        ! Initialize values to zero 
        lgs%a_index = 0.0
        lgs%a_value = 0.0 
        lgs%b_value = 0.0 
        lgs%x_value = 0.0 

        lgs%a_ptr(1) = 1

        k = 0

        do n = 1, lgs%nmax

            ! Get i,j indices of current point
            i = lgs%n2i(n)
            j = lgs%n2j(n)

            nr = n   ! row counter

            ! Get neighbor indices assuming periodic domain
            ! (all other boundary conditions are treated with special cases below)
            im1 = i-1
            if (im1 .eq. 0)    im1 = nx 
            ip1 = i+1
            if (ip1 .eq. nx+1) ip1 = 1 

            jm1 = j-1
            if (jm1 .eq. 0)    jm1 = ny
            jp1 = j+1
            if (jp1 .eq. ny+1) jp1 = 1 
            

            ! Handle special cases first, otherwise populate with normal inner discretization

            if (mask(i,j) .eq. 0) then 
                ! Zero value imposed 

                k = k+1
                lgs%a_index(k) = nr
                lgs%a_value(k) = 1.0_wp   ! diagonal element only
                
                lgs%b_value(nr) = 0.0_wp
                lgs%x_value(nr) = 0.0_wp

            else if (mask(i,j) .eq. -1) then 
                ! Prescribed value imposed 

                k = k+1
                lgs%a_index(k) = nr
                lgs%a_value(k) = 1.0_wp   ! diagonal element only
                
                lgs%b_value(nr) = Hbnd(i,j)
                lgs%x_value(nr) = H(i,j)

            else if ( (.not. trim(bcs(1)) .eq. "periodic") .and. i .eq. nx) then
                ! Right border

                if (bcs(1) .eq. "infinite") then

                    k = k+1
                    lgs%a_index(k) = lgs%ij2n(i,j)                  ! column counter for H(i,j)
                    lgs%a_value(k) =  1.0_wp
                    
                    k = k+1
                    lgs%a_index(k) = lgs%ij2n(im1,j)                ! column counter for H(im1,j)
                    lgs%a_value(k) = -1.0_wp

                    lgs%b_value(nr) = 0.0_wp
                    lgs%x_value(nr) = H(i,j)

                else
                    ! Assume zero for now

                    k = k+1
                    lgs%a_index(k) = nr
                    lgs%a_value(k) = 1.0_wp   ! diagonal element only
                    
                    lgs%b_value(nr) = 0.0_wp
                    lgs%x_value(nr) = 0.0_wp

                end if

            else if ( (.not. trim(bcs(2)) .eq. "periodic") .and. j .eq. ny) then
                ! Top border

                if (bcs(2) .eq. "infinite") then

                    k = k+1
                    lgs%a_index(k) = lgs%ij2n(i,j)                  ! column counter for H(i,j)
                    lgs%a_value(k) =  1.0_wp
                    
                    k = k+1
                    lgs%a_index(k) = lgs%ij2n(i,jm1)                ! column counter for H(i,jm1)
                    lgs%a_value(k) = -1.0_wp

                    lgs%b_value(nr) = 0.0_wp
                    lgs%x_value(nr) = H(i,j)

                else
                    ! Assume zero for now

                    k = k+1
                    lgs%a_index(k) = nr
                    lgs%a_value(k) = 1.0_wp   ! diagonal element only
                    
                    lgs%b_value(nr) = 0.0_wp
                    lgs%x_value(nr) = 0.0_wp

                end if

            else if ( (.not. trim(bcs(3)) .eq. "periodic") .and. i .eq. 1) then
                ! Left border

                if (bcs(3) .eq. "infinite") then

                    k = k+1
                    lgs%a_index(k) = lgs%ij2n(i,j)                  ! column counter for H(i,j)
                    lgs%a_value(k) =  1.0_wp
                    
                    k = k+1
                    lgs%a_index(k) = lgs%ij2n(ip1,j)                ! column counter for H(ip1,j)
                    lgs%a_value(k) = -1.0_wp

                    lgs%b_value(nr) = 0.0_wp
                    lgs%x_value(nr) = H(i,j)

                    ! Initial guess == previous H

                    lgs%x_value(nr) = H(i,j)                            
                

                else
                    ! Assume zero for now

                    k = k+1
                    lgs%a_index(k) = nr
                    lgs%a_value(k) = 1.0_wp   ! diagonal element only
                    
                    lgs%b_value(nr) = 0.0_wp
                    lgs%x_value(nr) = 0.0_wp

                end if

            else if ( (.not. trim(bcs(4)) .eq. "periodic") .and. j .eq. 1) then
                ! Bottom border

                if (bcs(4) .eq. "infinite") then

                    k = k+1
                    lgs%a_index(k) = lgs%ij2n(i,j)                  ! column counter for H(i,j)
                    lgs%a_value(k) =  1.0_wp
                    
                    k = k+1
                    lgs%a_index(k) = lgs%ij2n(i,jp1)                ! column counter for H(i,jp1)
                    lgs%a_value(k) = -1.0_wp

                    lgs%b_value(nr) = 0.0_wp
                    lgs%x_value(nr) = H(i,j)

                else
                    ! Assume zero for now

                    k = k+1
                    lgs%a_index(k) = nr
                    lgs%a_value(k) = 1.0_wp   ! diagonal element only
                    
                    lgs%b_value(nr) = 0.0_wp
                    lgs%x_value(nr) = 0.0_wp

                end if

            else
                ! Inner point
                ! Populate with advection and diffusion terms

                alpha = kappa(i,j) * dt_dx2
                beta  = kappa(i,j) * dt_dy2

                k = k+1
                lgs%a_index(k) = lgs%ij2n(i,jm1)                    ! for H(i,jm1)
                if (uy_1(i,j) > 0.0) &                              !   -Advection
                    lgs%a_value(k) = lgs%a_value(k) - dt_darea*uy_1(i,j)*dx*WOVI
                lgs%a_value(k) = lgs%a_value(k) - beta*WOVI         !   -Diffusion

                k = k+1
                lgs%a_index(k) = lgs%ij2n(im1,j)                    ! for H(im1,j)
                if (ux_1(i,j) > 0.0) &                              !   -Advection
                    lgs%a_value(k) = lgs%a_value(k) - dt_darea*ux_1(i,j)*dy*WOVI
                lgs%a_value(k) = lgs%a_value(k) - alpha*WOVI        !   -Diffusion

                k = k+1
                lgs%a_index(k) = nr                                 ! for H(i,j)
                lgs%a_value(k) = 1.0                                ! (diagonal element)
                if (uy_1(i,j) < 0.0) &                              !   -Advection
                    lgs%a_value(k) = lgs%a_value(k) &
                                    - dt_darea*uy_1(i,j)*dx*WOVI
                if (ux_1(i,j) < 0.0) &                              !   -Advection
                    lgs%a_value(k) = lgs%a_value(k) &
                                    - dt_darea*ux_1(i,j)*dy*WOVI
                if (ux_2(i,j) > 0.0) &
                    lgs%a_value(k) = lgs%a_value(k) &               !   -Advection
                                    + dt_darea*ux_2(i,j)*dy*WOVI
                if (uy_2(i,j) > 0.0) &
                    lgs%a_value(k) = lgs%a_value(k) &               !   -Advection
                                    + dt_darea*uy_2(i,j)*dx*WOVI
                lgs%a_value(k) = lgs%a_value(k) + (2.0*alpha+2.0*beta)*WOVI  !   -Diffusion

                k = k+1
                lgs%a_index(k) = lgs%ij2n(ip1,j)                    ! for H(ip1,j)
                if (ux_2(i,j) < 0.0) &                              !   -Advection
                    lgs%a_value(k) = lgs%a_value(k) + dt_darea*ux_2(i,j)*dy*WOVI
                lgs%a_value(k) = lgs%a_value(k) - alpha*WOVI        !   -Diffusion

                k = k+1
                lgs%a_index(k) = lgs%ij2n(i,jp1)                    ! for H(i,jp1)
                if (uy_2(i,j) < 0.0) &                              !   -Advection
                    lgs%a_value(k) = lgs%a_value(k) + dt_darea*uy_2(i,j)*dx*WOVI
                lgs%a_value(k) = lgs%a_value(k) - beta*WOVI         !   -Diffusion

                ! Right-hand side 

                lgs%b_value(nr) = H(i,j) + dt*F(i,j)
                
                lgs%b_value(nr) = lgs%b_value(nr) + &               !   -Advection (explicit)
                                -(1.0-WOVI) * dt_darea &            
                                     * (  ( ux_2(i,j)*Hx_2(i,j)*dy      &
                                           -ux_1(i,j)*Hx_1(i,j)*dy )    &
                                        + ( uy_2(i,j)*Hy_2(i,j)*dx      &
                                           -uy_1(i,j)*Hy_1(i,j)*dx ) )

                lgs%b_value(nr) = lgs%b_value(nr) + &               !   -Diffusion (explicit)
                                -(1.0-WOVI) * &
                                    ( -alpha*(H(ip1,j)-2.0*H(i,j)+H(im1,j)) &
                                      -beta *(H(i,jp1)-2.0*H(i,j)+H(i,jm1)) )

                ! Initial guess == previous H

                lgs%x_value(nr) = H(i,j)                            
                
            end if

            lgs%a_ptr(nr+1) = k+1   ! row is completed, store index to next row
            
        end do

        ! Done: A, x and b matrices in Ax=b have been populated 
        ! and stored in lgs object. 

        return

    end subroutine linear_solver_matrix_diffusion_advection_csr_2D

    subroutine get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,bcs)
        ! Boundary conditions (bcs) counterclockwise unit circle 
        ! 1: x, right-border
        ! 2: y, upper-border 
        ! 3: x, left--border 
        ! 4: y, lower-border 
        
        implicit none

        integer, intent(OUT) :: im1 
        integer, intent(OUT) :: ip1 
        integer, intent(OUT) :: jm1 
        integer, intent(OUT) :: jp1 
        integer, intent(IN)  :: i 
        integer, intent(IN)  :: j
        integer, intent(IN)  :: nx 
        integer, intent(IN)  :: ny
        character(len=*), intent(IN) :: bcs(4)

        if (trim(bcs(1)) .eq. "infinite") then
            ip1 = min(i+1,nx) 
        else    ! periodic, zero
            ip1 = i+1
            if (ip1 .eq. nx+1) ip1 = 1 
        end if 

        if (trim(bcs(2)) .eq. "infinite") then
            jp1 = min(j+1,ny) 
        else    ! periodic, zero
            jp1 = j+1
            if (jp1 .eq. ny+1) jp1 = 1 
        end if 

        if (trim(bcs(3)) .eq. "infinite") then
            im1 = max(i-1,1)
        else    ! periodic, zero
            im1 = i-1
            if (im1 .eq. 0) im1 = nx
        end if 

        if (trim(bcs(4)) .eq. "infinite") then
            jm1 = max(j-1,1)
        else    ! periodic, zero
            jm1 = j-1
            if (jm1 .eq. 0)    jm1 = ny
        end if 

        return

    end subroutine get_neighbor_indices

    subroutine set_boundary_conditions(bcs,bc,bc1,bc2,bc3,bc4)

        implicit none

        character(len=*), intent(OUT) :: bcs(4)
        character(len=*), intent(IN), optional :: bc
        character(len=*), intent(IN), optional :: bc1
        character(len=*), intent(IN), optional :: bc2
        character(len=*), intent(IN), optional :: bc3
        character(len=*), intent(IN), optional :: bc4

        ! Safety check, make sure BCs were provided
        ! Must have at least either the general bc, or when not given,
        ! then all four specific bcs.
        if (.not. present(bc) .and. &
                ( .not. present(bc1) .and. .not. present(bc2) &
            .and. .not. present(bc3) .and. .not. present(bc4) )    ) then
        
            write(error_unit,*) "set_boundary_conditions:: Error: Boundary condition &
            &arguments not given properly. Must either be provided with general bc, or, when &
            &not provided, then all four specific bcs."
            stop "set_boundary_conditions:: Error."

        end if

        ! If general bc provided, assign it to all four boundaries
        if (present(bc)) then            
            bcs(1:4) = trim(bc)
        end if

        ! Overwrite individual bcs if needed 
        if (present(bc1)) bcs(1) = trim(bc1)
        if (present(bc2)) bcs(2) = trim(bc2)
        if (present(bc3)) bcs(3) = trim(bc3)
        if (present(bc4)) bcs(4) = trim(bc4)

        return

    end subroutine set_boundary_conditions

    function adv2D_timestep(dx,dy,vx_max, vy_max) result(dt)
        implicit none 
        real(wp) :: dx, dy, vx_max, vy_max 
        real(wp) :: dt 

        ! Condition v*dt/dx <= 1 
        ! dt = dx / v
        dt = min(dx/max(vx_max,1d-6),dy/max(vy_max,1d-6))

        return 

    end function adv2D_timestep 

    function diff2D_timestep(dx,dy,kappa) result(dt)
        implicit none 
        real(wp) :: dx, dy, kappa 
        real(wp) :: dt 

        ! 1D condition kappa*dt/dx^2 <= 1
        ! dt = dx^2/kappa 

        dt = (1.d0/(2.d0*kappa)) *(dx*dy)**2/(dx**2+dy**2)
        return 

    end function diff2D_timestep 
    
    function diff2Dadi_timestep(dx,dy,kappa) result(dt)
        implicit none 
        real(wp) :: dx, dy, kappa 
        real(wp) :: dt 

        ! Condition: 4*kappa*dt / (dx*dx) <= 1

        dt = (dx*dy) / (4.d0*kappa)

        return 

    end function diff2Dadi_timestep 

!  === ADI related routines that need refactoring, probably with LIS ===




    subroutine solve2D(uu,u0,F,relax,alpha,gamma,beta,solver)

        implicit none 

        real(wp), intent(INOUT), dimension(:,:) :: uu
        real(wp), intent(IN),    dimension(:,:) :: u0, F, relax, alpha 
        real(wp), Intent(IN)                    :: gamma, beta
        character(len=*) :: solver 

        select case(solver)
            case("explicit")
                call solve_explicit2D(uu,u0,F,relax,alpha,gamma,beta)

            case("adi")
                call solve_adi2D(uu,u0,F,relax,alpha,gamma,beta)

            case DEFAULT

                write(*,*) "solvers:: solve2D:: ", "error solver type not defined: "//trim(solver)
                stop 

        end select

        return

    end subroutine solve2D 

    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Subroutine :  e x p l i c i t
    ! Purpose    :  Perform 1 timestep of explicit routine on 2d grid
    ! Author     :  Alex Robinson (24. June 2008)
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine solve_explicit2D(uu,u0,F,relax,alpha,gamma,beta)

        implicit none

        real(wp), intent(INOUT), dimension(:,:) :: uu
        real(wp), intent(IN),    dimension(:,:) :: u0, F, relax, alpha 
        real(wp), Intent(IN)                    :: gamma, beta

        real(wp), allocatable,   dimension(:,:) :: utmp, alpha1, beta1, relax1
        real(wp) :: gamma1 

        integer :: nx, ny 
        integer :: i, j, k    
        
        nx = size(uu,1)
        ny = size(uu,2)

        allocate(utmp(nx,ny),alpha1(nx,ny),beta1(nx,ny),relax1(nx,ny))

        ! Constants
        alpha1 = alpha 
        gamma1  = gamma
        beta1  = 1.d0-4.d0*alpha1 - gamma1*beta + gamma1*relax
        relax1 = gamma1*relax*u0 ! Assign 'boundary' values to relaxation matrix

        ! Update boundary conditions
        uu(1:nx,1)  = u0(1:nx,1)
        uu(1:nx,ny) = u0(1:nx,ny)
        uu(1,1:ny)  = u0(1,1:ny)
        uu(nx,1:ny) = u0(nx,1:ny)

        ! Define matrix with current values for calculations (to fill in boundary regions)
        utmp = uu

        ! explicit scheme
        do i= 2, nx-1
            do j= 2, ny-1
            
                uu(i,j)  = beta1(i,j) * utmp(i,j)   &
                           + alpha1(i,j) * (utmp(i,j+1)+utmp(i,j-1)+utmp(i+1,j)+utmp(i-1,j)) &
                           + gamma1*F(i,j) - relax1(i,j)
            end do
        end do

        return

    end subroutine solve_explicit2D

    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Subroutine :  a d i
    ! Purpose    :  Perform 1 timestep of adi routine on 2d grid
    !               kappa varying in space (prescribed)
    ! Author     :  Alex Robinson (01. Aug 2008)
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine solve_adi2D(uu,u0,F,relax,alpha,gamma,beta)

        implicit none

        real(wp), intent(INOUT), dimension(:,:) :: uu
        real(wp), intent(IN),    dimension(:,:) :: u0, F, relax, alpha 
        real(wp), Intent(IN)                    :: gamma, beta

        real(wp), allocatable,   dimension(:,:) :: alpha1, coeffs, relax1
        real(wp)                                :: gamma1, beta_ex, beta_im 

        real(wp), allocatable,   dimension(:)   :: ax, bx, cx, uprimex, cprimex, mx, rx
        real(wp), allocatable,   dimension(:)   :: ay, by, cy, uprimey, cprimey, my, ry
        
        integer :: nx, ny 
        integer :: i, j, k    
        
        nx = size(uu,1)
        ny = size(uu,2)

        allocate(alpha1(nx,ny),coeffs(nx,ny),relax1(nx,ny))
        allocate(ax(nx),bx(nx),cx(nx),uprimex(nx),cprimex(nx),mx(nx),rx(nx))
        allocate(ay(ny),by(ny),cy(ny),uprimey(ny),cprimey(ny),my(ny),ry(ny))
        
        ! Constants, divide by two because it is 1/2 timestep
        alpha1  = alpha / 2.d0
        gamma1  = gamma / 2.d0
        beta_ex = 0.d0            ! Explicit coefficient, 0 if using implicit
        beta_im = beta            ! Implicit coefficient, 0 if using explicit

        relax1 = gamma1*relax*u0     ! Assign 'boundary' values to relaxation matrix
        coeffs = 1.d0 - 2.d0*alpha1 - gamma1*beta_ex + gamma1*relax

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
                         + alpha1(1,j) * (uu(1,j+1) + uu(1,j-1))    &
                         + gamma1*F(1,j) - relax1(1,j) ) / bx(1)     
          
            ! Define the three diagonal rows of the matrix: x
            ax(2:nx-1) = -alpha1(2:nx-1,j)
            bx(2:nx-1) = 1 + 2.d0 * alpha1(2:nx-1,j) + gamma1*beta_im
            cx(2:nx-1) = -alpha1(2:nx-1,j) 
                    
            do i = 2, nx

                cprimex(i) = cx(i) / (bx(i) - ax(i) * cprimex(i-1))

                ! RHS
                uprimex(i) = ( ( coeffs(i,j) * uu(i,j)      &
                             + alpha1(i,j) * (uu(i,j+1) + uu(i,j-1))      & 
                             + gamma1*F(i,j) - relax1(i,j) ) - ax(i)    &
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
                         + alpha1(i,1) * (uu(i+1,1) + uu(i-1,1))   &  
                         + gamma1*F(i,1) - relax1(i,1) ) / by(1)
          
            ! Define the three diagonal rows of the matrix: y
            ay(2:ny-1) = -alpha1(i,2:ny-1)
            by(2:ny-1) = 1 + 2.d0 * alpha1(i,2:ny-1) + gamma1*beta_im
            cy(2:ny-1) = -alpha1(i,2:ny-1)

            do j = 2, ny

                cprimey(j) = cy(j) / ( by(j) - ay(j) * cprimey(j-1) )

                ! RHS               
                uprimey(j) = ( ( coeffs(i,j) * uu(i,j)      &
                             + alpha1(i,j) * (uu(i+1,j) + uu(i-1,j))     &  
                             + gamma1*F(i,j) - relax1(i,j) ) - ay(j)    &
                             * uprimey(j-1)) / (by(j) - ay(j) * cprimey(j-1))                  

            end do
          
            do j = ny-1, 2, -1
                uu(i,j) = uprimey(j) - cprimey(j) * uu(i,j+1)
            end do
          
        end do  ! End of i-loop

        return

    end subroutine solve_adi2D
    
    subroutine solve_diff_2D_adi(uu,u0,F,relax,dx,dy,dt,kappa,k_relax)

        implicit none 

        real(wp), dimension(:,:) :: uu, u0, F, kappa
        integer,          dimension(:,:) :: relax 
        real(wp) :: dx, dy, dt, k_relax 
        real(wp), dimension(:,:), allocatable :: utmp, udiff
        integer :: nx, ny 
        real(wp), dimension(:,:), allocatable :: alpha 
        real(wp) :: gamma, beta 

        nx = size(uu,1)
        ny = size(uu,2) 
        
        allocate(utmp(nx,ny),udiff(nx,ny))
        allocate(alpha(nx,ny)) 

        ! Constants, in time units of sec
        alpha = kappa * dt / (dx*dy)
        gamma = dt 
        beta  = 0

        utmp  = uu 
        udiff = uu 
!         call diffuse2D(udiff,dx,dy)
!         uu = utmp + dt*(kappa*udiff + F) - k_relax*relax*(utmp-u0)

!         call solve_adi2D(udiff,u0,F,k_relax*relax,alpha,gamma,beta)
!         uu = udiff 
        
        call solve_adi2D(udiff,u0,F,0.0*relax,alpha,gamma,beta)
        uu = udiff - dt*k_relax*relax*(udiff-u0)

        return

    end subroutine solve_diff_2D_adi  

    subroutine set_border(uu,u0,u00,nbx,nby)
        
        implicit none 
        
        real(wp), dimension(:,:) :: uu
        real(wp), dimension(:,:), optional :: u0
        real(wp), optional :: u00 
        integer, optional :: nbx, nby 
        integer :: nx, ny, nbx1, nby1 

        nbx1 = 1
        if (present(nbx)) nbx1 = nbx 
        nby1 = 1 
        if (present(nby)) nby1 = nby

        nx = size(uu,1)
        ny = size(uu,2) 

        if (present(u0)) then 
            uu(1:nbx1,:)       = u0(1:nbx1,:)
            uu(nx-nbx1+1:nx,:) = u0(nx-nbx1+1:nx,:)
            uu(:,1:nby1)       = u0(:,1:nby1)
            uu(:,ny-nby1+1:ny) = u0(:,ny-nby1+1:ny)
        else if (present(u00)) then 
            uu(1:nbx1,:)       = u00
            uu(nx-nbx1+1:nx,:) = u00
            uu(:,1:nby1)       = u00
            uu(:,ny-nby1+1:ny) = u00
        else
            write(*,*) "set_border:: either u0 or u00 should be defined as an argument."
            stop 
        end if 

        return

    end subroutine set_border 

end module solvers
