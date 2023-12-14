module solvers

    use precision, only : wp
    use solver_linear 

    implicit none
    
    private
    public :: solve_diffusion_advection_2D
    public :: solve2D
    public :: solve_diff_2D_adi, diff2Dadi_timestep
    public :: adv2D_timestep, diff2D_timestep
    public :: set_border

contains

    subroutine solve_diffusion_advection_2D(uu,v_x,v_y,F,kappa,ubnd,mask,dx,dy,dt,k_rel,solver,step)

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

        nx = size(uu,1)
        ny = size(uu,2)

        select case(trim(step))

            case("fe")
                ! Forward-euler timestepping

                allocate(f0(nx,ny))

                ! Calculate derivative
                call calc_tendency(f0,uu,v_x,v_y,F,kappa,ubnd,mask,dx,dy,dt,k_rel,solver)

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

                call calc_tendency(f0,u0,v_x,v_y,F,kappa,ubnd,mask,dx,dy,dt,k_rel,solver)

                u1 = u0 + (dt/2.0)*f0

                call calc_tendency(f1,u1,v_x,v_y,F,kappa,ubnd,mask,dx,dy,dt,k_rel,solver)

                u2 = u0 + (dt/2.0)*f1
                
                call calc_tendency(f2,u2,v_x,v_y,F,kappa,ubnd,mask,dx,dy,dt,k_rel,solver)

                u3 = u0 + dt*f2

                call calc_tendency(f3,u3,v_x,v_y,F,kappa,ubnd,mask,dx,dy,dt,k_rel,solver)

                ! Combine them to estimate the solution U at time t+dt
                uu = u0 + dt * ( f0 + 2.0*f1 + 2.0*f2 + f3 ) / 6.0

            case DEFAULT

                write(*,*) "solve_diffusion_advection_2D:: timestepping method not recongized."
                write(*,*) "step = ", trim(step)
                stop 
                
        end select
        
        return

    end subroutine solve_diffusion_advection_2D

    subroutine calc_tendency(dudt,uu,v_x,v_y,F,kappa,ubnd,mask,dx,dy,dt,k_rel,solver)
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
        

        ! Local variables
        integer  :: i, j, nx, ny 
        real(wp) :: diff 
        real(wp), allocatable :: dudt_adv(:,:)
        real(wp), allocatable :: dudt_diff(:,:)
        real(wp), allocatable :: dudt_relax(:,:)
        integer,  allocatable :: mask_adv(:,:)
        
        nx = size(dudt,1)
        ny = size(dudt,2) 
        
        allocate(dudt_diff(nx,ny))
        allocate(dudt_adv(nx,ny))
        allocate(dudt_relax(nx,ny))
        allocate(mask_adv(nx,ny))

        ! Define advection mask (used for implicit solver only for now)
        mask_adv = 1
        where (mask .eq. 1) mask_adv = -1
        mask_adv(1,:)   = -1
        mask_adv(nx,:)  = -1
        mask_adv(:,1)   = -1
        mask_adv(:,ny)  = -1
        
        select case(trim(solver))

            case("expl-diff")
                ! Explicit diffusion

                ! Get diffusive tendency, i.e. kappa*(Laplacian operator) (du/dt)
                call calc_tendency_diffuse2D(dudt_diff,uu,kappa,dx,dy)

                ! No advection considered
                dudt_adv = 0.0

                ! Get relaxation rate at boundaries (k_rel is restoring rate fraction per second, positive value)
                dudt_relax = -k_rel*mask*(uu-ubnd)

                ! Get total tendency for this timestep
                dudt = dudt_diff - dudt_adv + F + dudt_relax

            case("expl-adv")
                ! Explicit advection 

                ! No diffusion considered
                dudt_diff = 0.0 

                ! Get advective tendency (du/dt)
                call calc_tendency_advec2D_upwind(dudt_adv,uu,v_x,v_y,dx,dy)

                ! Get relaxation rate at boundaries (k_rel is restoring rate fraction per second, positive value)
                dudt_relax = -k_rel*mask*(uu-ubnd)

                ! Get total tendency for this timestep
                dudt = dudt_diff - dudt_adv + F + dudt_relax
                
            case("expl")
                ! Explicit diffusion and advection 

                ! Get diffusive tendency, i.e. kappa*(Laplacian operator) (du/dt)
                call calc_tendency_diffuse2D(dudt_diff,uu,kappa,dx,dy)

                ! Get advective tendency (du/dt)
                call calc_tendency_advec2D_upwind(dudt_adv,uu,v_x,v_y,dx,dy)

                ! Get relaxation rate at boundaries (k_rel is restoring rate fraction per second, positive value)
                dudt_relax = -k_rel*mask*(uu-ubnd)

                ! Get total tendency for this timestep
                dudt = dudt_diff - dudt_adv + F + dudt_relax
                
            case("expl-impl")
                ! Explicit diffusion, implicit advection

                ! Get diffusive tendency, i.e. kappa*(Laplacian operator) (du/dt)
                call calc_tendency_diffuse2D(dudt_diff,uu,kappa,dx,dy)
                !dudt_diff = 0.0 

                ! Get advective tendency (du/dt)
                call calc_tendency_advec2D_upwind_impl(dudt_adv,uu,v_x,v_y,F*0.0,mask_adv, &
                                                    dx,dy,dt,solver="impl-lis",boundaries="infinite")
                !dudt_adv = 0.0

                ! Get relaxation rate at boundaries (k_rel is restoring rate fraction per second, positive value)
                dudt_relax = -k_rel*mask*(uu-ubnd)

                ! Get total tendency for this timestep
                dudt = dudt_diff - dudt_adv + F + dudt_relax
                
            case("impl")
                ! Implicit diffusion and advection 

                ! TO DO 
                write(*,*) "calc_tendency:: method not yet implemented."
                write(*,*) "solver = ", trim(solver)
                stop 

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

    subroutine calc_tendency_diffuse2D(dudt,uu,kappa,dx,dy)
        ! Explicit calculation of kappa * 2D Laplace: kappa*d2u/dxdy
        ! If input units are [u], returns [u/m2]
    
        implicit none 

        real(wp), intent(OUT) :: dudt(:,:)
        real(wp), intent(IN)  :: uu(:,:)
        real(wp), intent(IN)  :: kappa(:,:)
        real(wp), intent(IN)  :: dx
        real(wp), intent(IN)  :: dy

        ! Local variables 
        integer :: i, j, nx, ny 
        real(wp) :: inv_dx2, inv_dy2
        real(wp) :: du

        nx = size(dudt,1)
        ny = size(dudt,2) 

        inv_dx2 = 1.d0 / (dx*dx)
        inv_dy2 = 1.d0 / (dy*dy)

        dudt = 0.d0 

        do j = 2,ny-1
        do i = 2,nx-1
        
            ! Laplacian, five-point stencil finite difference method (du/(dx*dy))
            du =  inv_dx2*(uu(i-1,j)-2.d0*uu(i,j)+uu(i+1,j)) &
                + inv_dy2*(uu(i,j-1)-2.d0*uu(i,j)+uu(i,j+1))
            
            ! Get tendency (du/dt)
            dudt(i,j) = kappa(i,j)*du

        end do 
        end do

        return

    end subroutine calc_tendency_diffuse2D 

    
    subroutine calc_tendency_advec2D_upwind(dudt,uu,vx,vy,dx,dy)
        ! Second-order upwind advection scheme
        ! If input units are [u], returns [u/s]
        
        implicit none 

        real(wp), intent(OUT) :: dudt(:,:)
        real(wp), intent(IN)  :: uu(:,:)
        real(wp), intent(IN)  :: vx(:,:)
        real(wp), intent(IN)  :: vy(:,:)
        real(wp), intent(IN)  :: dx
        real(wp), intent(IN)  :: dy 

        ! Local variables
        integer :: i, j, nx, ny 
        real(wp) :: inv_2dx, inv_2dy
        real(wp), allocatable :: utmp(:,:)
        real(wp), allocatable :: du_x_neg(:,:)
        real(wp), allocatable :: du_x_pos(:,:)
        real(wp), allocatable :: du_y_neg(:,:)
        real(wp), allocatable :: du_y_pos(:,:) 
        
        nx = size(dudt,1)
        ny = size(dudt,2) 

        allocate(du_x_neg(nx,ny))
        allocate(du_x_pos(nx,ny))
        allocate(du_y_neg(nx,ny))
        allocate(du_y_pos(nx,ny))

        inv_2dx = 1.d0 / (2.d0*dx)
        inv_2dy = 1.d0 / (2.d0*dy)

        du_x_neg = 0.d0  
        du_x_pos = 0.d0 
        du_y_neg = 0.d0  
        du_y_pos = 0.d0 
        
        do j = 1, ny
            do i = 3,nx
                du_x_neg(i,j) = inv_2dx* ( 3.d0*uu(i,j)  -4.d0*uu(i-1,j)+1.d0*uu(i-2,j) )
            end do
            do i = 1,nx-2
                du_x_pos(i,j) = inv_2dx* (-1.d0*uu(i+2,j)+4.d0*uu(i+1,j)-3.d0*uu(i,j) )
            end do
        end do

        do i = 1,nx
            do j = 3,ny
                du_y_neg(i,j) = inv_2dy* ( 3.d0*uu(i,j)  -4.d0*uu(i,j-1)+1.d0*uu(i,j-2) )
            end do
            do j = 1,ny-2
                du_y_pos(i,j) = inv_2dy* (-1.d0*uu(i,j+2)+4.d0*uu(i,j+1)-3.d0*uu(i,j) )
            end do
        end do

        where(vx .le. 0.d0) du_x_neg = 0.d0 
        where(vx .gt. 0.d0) du_x_pos = 0.d0 
        where(vy .le. 0.d0) du_y_neg = 0.d0 
        where(vy .gt. 0.d0) du_y_pos = 0.d0 
            
        dudt =  (vx*du_x_neg + vx*du_x_pos) &
              + (vy*du_y_neg + vy*du_y_pos)

        return

    end subroutine calc_tendency_advec2D_upwind






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
    


    ! ==== LIS-based solving routines ====

    subroutine calc_tendency_advec2D_upwind_impl(dvdt,var,ux,uy,var_dot,mask_adv,dx,dy,dt,solver,boundaries)
        ! General routine to apply 2D advection equation to variable `var` 
        ! with source term `var_dot`. Various solvers are possible

        real(wp),       intent(OUT)   :: dvdt(:,:)              ! [dvdt] Variable rate of change
        real(wp),       intent(IN)    :: var(:,:)               ! [var]  Variable to be advected
        real(wp),       intent(IN)    :: ux(:,:)                ! [m/a] 2D velocity, x-direction (ac-nodes)
        real(wp),       intent(IN)    :: uy(:,:)                ! [m/a] 2D velocity, y-direction (ac-nodes)
        real(wp),       intent(IN)    :: var_dot(:,:)           ! [dvar/dt] Source term for variable
        integer,        intent(IN)    :: mask_adv(:,:)          ! Advection mask
        real(wp),       intent(IN)    :: dx                     ! [m]   Horizontal resolution, x-direction
        real(wp),       intent(IN)    :: dy                     ! [m]   Horizontal resolution, y-direction
        real(wp),       intent(IN)    :: dt                     ! [a]   Timestep 
        character(len=*), intent(IN)    :: solver               ! Solver to use for the ice thickness advection equation
        character(len=*), intent(IN)    :: boundaries           ! Boundary conditions to impose

        ! Local variables
        integer :: nx, ny  
        type(linear_solver_class) :: lgs
        character(len=256)        :: adv_lis_opt 

        real(wp), allocatable :: var_now(:,:) 

        nx = size(var,1)
        ny = size(var,2)

        allocate(var_now(nx,ny))

        ! Assign local variable to be modified 
        var_now = var 

        select case(trim(solver))
            ! Choose solver to use 

            case("none") 

                ! Do nothing: no advection considered. 

            case("impl-lis")

                ! Initialize linear solver variables
                call linear_solver_init(lgs,nx,ny,nvar=1,n_terms=5)

                ! Populate advection matrices Ax=b
                call linear_solver_matrix_advection_csr_2D(lgs,var_now,ux,uy,var_dot,mask_adv,dx,dy,dt,boundaries)
                
                ! Solve linear equation
                adv_lis_opt = "-i bicg -p ilu -maxiter 1000 -tol 1.0e-12 -initx_zeros false"
                !adv_lis_opt = "-i minres -p jacobi -maxiter 1000 -tol 1.0e-12 -initx_zeros false"
                call linear_solver_matrix_solve(lgs,adv_lis_opt)
                
                !call linear_solver_print_summary(lgs,io_unit_err)

                ! Store advection solution
                call linear_solver_save_advection(var_now,lgs)

            case DEFAULT 

                write(*,*) "calc_tendency_advec2D_upwind_impl:: Error: solver not recognized."
                write(*,*) "solver = ", trim(solver)
                stop 

        end select 
        
        ! Determine rate of change 
        dvdt = (var_now - var) / dt 

        return 

    end subroutine calc_tendency_advec2D_upwind_impl

    subroutine linear_solver_save_advection(H,lgs)
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

    end subroutine linear_solver_save_advection

    subroutine linear_solver_matrix_advection_csr_2D(lgs,H,ux,uy,F,mask,dx,dy,dt,boundaries)
        ! Define sparse matrices A*x=b in format 'compressed sparse row' (csr)
        ! for 2D advection equations with velocity components
        ! ux and uy defined on ac-nodes (right and top borders of i,j grid cell)
        ! and variable to be advected H defined on aa nodes.
        ! Store sparse matrices in linear_solver_class object 'lgs' for later use.

        implicit none 

        type(linear_solver_class), intent(INOUT) :: lgs
        real(wp), intent(INOUT)   :: H(:,:)         ! [X] Variable of interest (aa nodes)
        real(wp), intent(IN)      :: ux(:,:)        ! [m a-1] Horizontal velocity x-direction (ac nodes)
        real(wp), intent(IN)      :: uy(:,:)        ! [m a-1] Horizontal velocity y-direction (ac nodes)
        real(wp), intent(IN)      :: F(:,:)         ! [m a-1] Net source/sink terms (aa nodes)
        integer,  intent(IN)      :: mask(:,:)      ! Advection mask
        real(wp), intent(IN)      :: dx             ! [m] Horizontal step x-direction
        real(wp), intent(IN)      :: dy             ! [m] Horizontal step y-direction 
        real(wp), intent(IN)      :: dt             ! [a] Time step 
        character(len=*), intent(IN) :: boundaries 

        ! Local variables  
        integer  :: i, j, k, nx, ny
        integer  :: im1, ip1, jm1, jp1 
        integer  :: n, nr, nc
        real(wp) :: dt_darea
        character(len=56) :: bcs(4)

        real(wp), allocatable  :: ux_1(:,:), ux_2(:,:)
        real(wp), allocatable  :: uy_1(:,:), uy_2(:,:)
        real(wp), allocatable  :: Hx_1(:,:), Hx_2(:,:)
        real(wp), allocatable  :: Hy_1(:,:), Hy_2(:,:)

        real(wp), parameter :: WOVI = 1.0     ! Weighting parameter for the over-implicit scheme 

        nx = size(H,1)
        ny = size(H,2) 

        dt_darea = dt/(dx*dy)

        ! Boundary conditions (bcs) counterclockwise unit circle 
        ! 1: x, right-border
        ! 2: y, upper-border 
        ! 3: x, left--border 
        ! 4: y, lower-border 
        
        ! Define border conditions (only choices are: no-slip, free-slip, periodic)
        select case(trim(boundaries)) 

            case("infinite")

                bcs(1:4) = "infinite" 
            
            case("zeros")

                bcs(1:4) = "zero" 

            case("periodic","periodic-xy")

                bcs(1:4) = "periodic" 

            case DEFAULT 

                bcs(1:4) = "zero"

        end select 


        ! Safety check for initialization
        if (.not. allocated(lgs%x_value)) then 
            ! Object 'lgs' has not been initialized yet, do so now.

            call linear_solver_init(lgs,nx,ny,nvar=1,n_terms=5)

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
            call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)

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

        !-------- Assembly of the system of linear equations
        !             (matrix storage: compressed sparse row CSR) --------

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
                
                lgs%b_value(nr) = 0.0_wp
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

                k = k+1
                lgs%a_index(k) = lgs%ij2n(i,jm1)                    ! for H(i,jm1)
                if (uy_1(i,j) > 0.0) &
                    lgs%a_value(k) = -dt_darea*uy_1(i,j)*dx*WOVI

                k = k+1
                lgs%a_index(k) = lgs%ij2n(im1,j)                    ! for H(im1,j)
                if (ux_1(i,j) > 0.0) &
                    lgs%a_value(k) = -dt_darea*ux_1(i,j)*dy*WOVI

                k = k+1
                lgs%a_index(k) = nr                                 ! for H(i,j)
                lgs%a_value(k) = 1.0                                ! (diagonal element)
                if (uy_1(i,j) < 0.0) &
                    lgs%a_value(k) = lgs%a_value(k) &
                                    - dt_darea*uy_1(i,j)*dx*WOVI
                if (ux_1(i,j) < 0.0) &
                    lgs%a_value(k) = lgs%a_value(k) &
                                    - dt_darea*ux_1(i,j)*dy*WOVI
                if (ux_2(i,j) > 0.0) &
                    lgs%a_value(k) = lgs%a_value(k) &
                                    + dt_darea*ux_2(i,j)*dy*WOVI
                if (uy_2(i,j) > 0.0) &
                    lgs%a_value(k) = lgs%a_value(k) &
                                    + dt_darea*uy_2(i,j)*dx*WOVI

                k = k+1
                lgs%a_index(k) = lgs%ij2n(ip1,j)                    ! for H(ip1,j)
                if (ux_2(i,j) < 0.0) &
                    lgs%a_value(k) = dt_darea*ux_2(i,j)*dy*WOVI

                k = k+1
                lgs%a_index(k) = lgs%ij2n(i,jp1)                    ! for H(i,jp1)
                if (uy_2(i,j) < 0.0) &
                    lgs%a_value(k) = dt_darea*uy_2(i,j)*dx*WOVI


                ! Right-hand side 

                lgs%b_value(nr) = H(i,j) + dt*F(i,j) &
                                -(1.0-WOVI) * dt_darea &
                                     * (  ( ux_2(i,j)*Hx_2(i,j)*dy      &
                                           -ux_1(i,j)*Hx_1(i,j)*dy )    &
                                        + ( uy_2(i,j)*Hy_2(i,j)*dx      &
                                           -uy_1(i,j)*Hy_1(i,j)*dx  ) )

                ! Initial guess == previous H

                lgs%x_value(nr) = H(i,j)                            
                
            end if

            lgs%a_ptr(nr+1) = k+1   ! row is completed, store index to next row
            
        end do

        ! Done: A, x and b matrices in Ax=b have been populated 
        ! and stored in lgs object. 

        return

    end subroutine linear_solver_matrix_advection_csr_2D

    subroutine get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)

        implicit none

        integer, intent(OUT) :: im1 
        integer, intent(OUT) :: ip1 
        integer, intent(OUT) :: jm1 
        integer, intent(OUT) :: jp1 
        integer, intent(IN)  :: i 
        integer, intent(IN)  :: j
        integer, intent(IN)  :: nx 
        integer, intent(IN)  :: ny
        
        character(len=*), intent(IN) :: boundaries

        select case(trim(boundaries))

            case("infinite")
                im1 = max(i-1,1)
                ip1 = min(i+1,nx) 
                jm1 = max(j-1,1)
                jp1 = min(j+1,ny) 

            case("MISMIP3D","TROUGH")
                im1 = max(i-1,1)
                ip1 = min(i+1,nx) 
                jm1 = j-1
                if (jm1 .eq. 0)    jm1 = ny
                jp1 = j+1
                if (jp1 .eq. ny+1) jp1 = 1 
                
            case DEFAULT 
                ! Periodic

                im1 = i-1
                if (im1 .eq. 0)    im1 = nx 
                ip1 = i+1
                if (ip1 .eq. nx+1) ip1 = 1 

                jm1 = j-1
                if (jm1 .eq. 0)    jm1 = ny
                jp1 = j+1
                if (jp1 .eq. ny+1) jp1 = 1 

        end select 

        return

    end subroutine get_neighbor_indices

end module solvers
