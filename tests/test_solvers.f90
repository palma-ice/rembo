program test_solvers

    use precision
    use solvers
    use ncio 

    implicit none

    character(len=512) :: filename = "test.nc"

    integer,  parameter :: nx = 101
    integer,  parameter :: ny = 101
    real(wp), parameter :: dx = 10e3
    real(wp), parameter :: dy = 10e3
    real(wp), parameter :: k_rel = 1e-1

    real(wp), parameter :: D = 1e12
    real(wp), parameter :: tsl_fac = 101100.0 /9.81*715.0 

    real(wp) :: uu(nx,ny)
    real(wp) :: F(nx,ny)
    real(wp) :: kappa(nx,ny)
    real(wp) :: ubnd(nx,ny)
    integer  :: mask(nx,ny)
    real(wp) :: v_x(nx,ny)
    real(wp) :: v_y(nx,ny)

    integer  :: iter, n_iter
    real(wp) :: time_init, time_end, dt 
    real(wp) :: time 
    real(wp) :: dt_out 

    integer :: nx_mid, ny_mid

    ! Get middle of domain
    nx_mid = (nx-1)/2
    ny_mid = (ny-1)/2

    ! Set relaxation mask
    mask = 0
    mask(1:5,:)     = 1
    mask(nx-5:nx,:) = 1
    mask(:,1:5)     = 1
    mask(:,ny-5:ny) = 1

    uu = 250.0
    uu(1:5,:)       = 273.15
    uu(nx-5:nx,:)   = 273.15
    uu(:,1:5)       = 273.15
    uu(:,ny-5:ny)   = 273.15

    ubnd = 273.15

    ! Define forcing field
    F = 0.0
    F(nx_mid-5:nx_mid+5,ny_mid-5:ny_mid+5) = 50.0 / tsl_fac

    ! Define velocities
    v_x = 0.0
    v_y = 10.0 

    ! Define kappa = D/ce
    kappa = D / tsl_fac

    time_init = 0.0
    time_end  = 10000.0 
    dt        = 10.0
    n_iter    = ceiling( (time_end-time_init)/dt ) + 1 

    dt_out    = 1000.0 

    time = time_init 

    ! Initialize the output file
    call write_init(filename,nx,ny,dx,dy,time_init)

    do iter = 1, n_iter 

        time = time_init + (iter-1)*dt

        call solve_diffusion_advection_2D(uu,v_x,v_y,F,kappa,ubnd,mask,dx,dy,dt,k_rel, &
                                                                    solver="expl",step="fe")

        call write_step(filename,uu,F,time)

        write(*,*) "time = ", time

    end do 

contains

    subroutine write_init(filename,nx,ny,dx,dy,time_init)

        implicit none

        character(len=*),  intent(IN) :: filename 
        integer,           intent(IN) :: nx
        integer,           intent(IN) :: ny
        real(wp),          intent(IN) :: dx
        real(wp),          intent(IN) :: dy
        real(wp),          intent(IN) :: time_init

        ! Local variables 

        ! Initialize netcdf file and dimensions
        call nc_create(filename)
        call nc_write_dim(filename,"xc",    x=0.0_wp,dx=dx*1e-3,nx=nx, units="kilometers")
        call nc_write_dim(filename,"yc",    x=0.0_wp,dx=dy*1e-3,nx=ny, units="kilometers")
        call nc_write_dim(filename,"time",  x=time_init,dx=1.0_wp,nx=1,units="seconds",unlimited=.TRUE.)
        
        return

    end subroutine write_init

    subroutine write_step(filename,uu,F,time)

        implicit none

        character(len=*),  intent(IN) :: filename 
        real(wp),          intent(IN) :: uu(:,:)
        real(wp),          intent(IN) :: F(:,:)
        real(wp),          intent(IN) :: time

        ! Local variables
        integer  :: ncid, n
        integer  :: nx, ny
        real(wp) :: time_prev 
        
        nx = size(uu,1)
        ny = size(uu,2)

        ! Open the file for writing
        call nc_open(filename,ncid,writable=.TRUE.)

        ! Determine current writing time step 
        n = nc_size(filename,"time",ncid)
        call nc_read(filename,"time",time_prev,start=[n],count=[1],ncid=ncid) 
        if (abs(time-time_prev).gt.1e-5) n = n+1 

        ! Update the time step
        call nc_write(filename,"time",time,dim1="time",start=[n],count=[1],ncid=ncid)

        ! Update the variables
        call nc_write(filename,"uu",uu,units="m",long_name="Field variable", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],count=[nx,ny,1],ncid=ncid)
        call nc_write(filename,"F",F,units="m",long_name="Field forcing", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],count=[nx,ny,1],ncid=ncid)

        ! Close the netcdf file
        call nc_close(ncid)

        return

    end subroutine write_step

end program test_solvers