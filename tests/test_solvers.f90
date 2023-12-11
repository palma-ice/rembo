program test_solvers

    use precision
    use solvers
    use ncio 

    implicit none

    character(len=512) :: filename = "test.nc"

    integer,  parameter :: nx = 101
    integer,  parameter :: ny = 101
    real(wp), parameter :: dx = 1.0
    real(wp), parameter :: dy = 1.0

    real(wp) :: uu(nx,ny)
    integer  :: mask(nx,ny)

    integer  :: iter, n_iter
    real(wp) :: time_init, time_end, dt 
    real(wp) :: time 

    ! Set relaxation mask
    mask = 0
    mask(1,:)  = 1
    mask(nx,:) = 1
    mask(:,1)  = 1
    mask(:,ny) = 1

    time_init = 0.0
    time_end  = 10.0 
    dt        = 1.0
    n_iter    = ceiling( (time_end-time_init)/dt ) + 1 

    time = time_init 

    uu = 0.0
    uu(1,:)  = 100.0
    uu(nx,:) = 100.0
    uu(:,1)  = 100.0
    uu(:,ny) = 100.0

    ! Initialize the output file
    call write_init(filename,nx,ny,dx,dy,time_init)

    do iter = 1, n_iter 

        time = time_init + (iter-1)*dt

        uu = real(iter,wp)

        call write_step(filename,uu,time)

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
        call nc_write_dim(filename,"xc",    x=0.0_wp,dx=dx,nx=nx, units="meters")
        call nc_write_dim(filename,"yc",    x=0.0_wp,dx=dy,nx=ny, units="meters")
        call nc_write_dim(filename,"time",  x=time_init,dx=1.0_wp,nx=1,units="seconds",unlimited=.TRUE.)
        
        return

    end subroutine write_init

    subroutine write_step(filename,uu,time)

        implicit none

        character(len=*),  intent(IN) :: filename 
        real(wp),          intent(IN) :: uu(:,:)
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

        ! Update the variable
        call nc_write(filename,"uu",uu,units="m",long_name="Field variable", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],count=[nx,ny,1],ncid=ncid)

        ! Close the netcdf file
        call nc_close(ncid)
        
        return

    end subroutine write_step

end program test_solvers