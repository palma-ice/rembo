module solvers
    
    implicit none
    
    integer,  parameter :: wp  = kind(1.0)

    private
    public :: solve2D
    public :: solve_diff_2D_adi, diff2Dadi_timestep
    public :: adv_diff_2D, advec2D_upwind, diffuse2D
    public :: adv2D_timestep, diff2D_timestep
    public :: set_border

contains

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
        uu = udiff - k_relax*relax*(udiff-u0)

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
    
    subroutine adv_diff_2D(uu,u0,F,relax,dx,dy,dt,kappa,k_relax,v_x,v_y)

        implicit none 

        real(wp), dimension(:,:) :: uu, u0, F, kappa
        real(wp), dimension(:,:), optional :: v_x, v_y 
        integer,          dimension(:,:) :: relax 
        real(wp) :: dx, dy, dt, k_relax 
        real(wp), dimension(:,:), allocatable :: utmp, uadv, udiff
        integer :: nx, ny 

        nx = size(uu,1)
        ny = size(uu,2) 
        
        allocate(utmp(nx,ny),udiff(nx,ny),uadv(nx,ny))

        utmp  = uu 

        udiff = uu 
        call diffuse2D(udiff,dx,dy)

        if (present(v_x) .and. present(v_y)) then
            uadv = uu 
            call advec2D_upwind(uadv,v_x,v_y,dx,dy)
        else
            uadv = 0.d0 
        end if 

        ! Update uu
        uu = utmp + dt*(kappa*udiff - uadv + F) - k_relax*relax*(utmp-u0)

        
        return

    end subroutine adv_diff_2D 

    ! Explicit calculation of 2D Laplace: d2u/dxdy
    ! If input units are [u], returns [u/m2]
    subroutine diffuse2D(uu,dx,dy)

        implicit none 

        real(wp), dimension(:,:) :: uu
        real(wp) :: dx, dy 
        real(wp) :: inv_dx2, inv_dy2
        integer :: i, j, nx, ny 
        real(wp), dimension(:,:), allocatable :: utmp 

        nx = size(uu,1)
        ny = size(uu,2) 

        allocate(utmp(nx,ny))
        utmp = uu 

        inv_dx2 = 1.d0 / (dx*dx)
        inv_dy2 = 1.d0 / (dy*dy)

        uu = 0.d0 

        do i = 2,nx-1
            do j = 2,ny-1
                           ! five-point stencil finite difference method
                uu(i,j) =  inv_dx2*(utmp(i-1,j)-2.d0*utmp(i,j)+utmp(i+1,j)) &
                         + inv_dy2*(utmp(i,j-1)-2.d0*utmp(i,j)+utmp(i,j+1))
            end do 
        end do

        return

    end subroutine diffuse2D 

    ! Second-order upwind advection scheme
    ! If input units are [u], returns [u/s]
    subroutine advec2D_upwind(uu,vx,vy,dx,dy)

        implicit none 

        real(wp), dimension(:,:) :: uu, vx, vy
        real(wp) :: dx, dy 
        real(wp) :: inv_2dx, inv_2dy
        integer :: i, j, nx, ny 
        real(wp), dimension(:,:), allocatable :: utmp
        real(wp), dimension(:,:), allocatable :: uu_x_neg, uu_x_pos 
        real(wp), dimension(:,:), allocatable :: uu_y_neg, uu_y_pos 
        
        nx = size(uu,1)
        ny = size(uu,2) 

        allocate(utmp(nx,ny))
        allocate(uu_x_neg(nx,ny),uu_x_pos(nx,ny))
        allocate(uu_y_neg(nx,ny),uu_y_pos(nx,ny))

        utmp     = uu 
        
        inv_2dx = 1.d0 / (2.d0*dx)
        inv_2dy = 1.d0 / (2.d0*dy)

        uu = 0.d0 
        uu_x_neg = 0.d0  
        uu_x_pos = 0.d0 
        uu_y_neg = 0.d0  
        uu_y_pos = 0.d0 
        
        do i = 3,nx-2
            do j = 3,ny-2
                uu_x_neg(i,j) = inv_2dx* ( 3.d0*utmp(i,j)  -4.d0*utmp(i-1,j)+1.d0*utmp(i-2,j) )
                uu_x_pos(i,j) = inv_2dx* (-1.d0*utmp(i+2,j)+4.d0*utmp(i+1,j)-3.d0*utmp(i,j) )

                uu_y_neg(i,j) = inv_2dy* ( 3.d0*utmp(i,j)  -4.d0*utmp(i,j-1)+1.d0*utmp(i,j-2) )
                uu_y_pos(i,j) = inv_2dy* (-1.d0*utmp(i,j+2)+4.d0*utmp(i,j+1)-3.d0*utmp(i,j) )
            end do 
        end do

        where(vx .le. 0.d0) uu_x_neg = 0.d0 
        where(vx .gt. 0.d0) uu_x_pos = 0.d0 
        where(vy .le. 0.d0) uu_y_neg = 0.d0 
        where(vy .gt. 0.d0) uu_y_pos = 0.d0 
            
        uu =  (vx*uu_x_neg + vx*uu_x_pos) &
            + (vy*uu_y_neg + vy*uu_y_pos)

        return

    end subroutine advec2D_upwind 

end module solvers
