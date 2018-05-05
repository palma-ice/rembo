
module sim_sampler 

    use ncio 
    use latinhypercube 

    implicit none 

    double precision, parameter :: MISSING = -9999.d0 

    type ensemble_summary_class 
        double precision, dimension(:), allocatable :: obs, var, err, err_std 
    end type 

    type ensemble_class 
        character(len=256)   :: name, dimname1 
        integer              :: nsim
        double precision, dimension(:),   allocatable :: dim1, obs, obs_err 
        double precision, dimension(:,:), allocatable :: var, err, err_std, wgt

        type(ensemble_summary_class) :: stat

    end type 

    interface sampler_latin 
        module procedure sampler_latin_1D, sampler_latin_ND 
    end interface 

    private
    public :: sampler, sampler_latin
    public :: stats_rmse, stats_stdev
    public :: ensemble_class, ensemble_init, ensemble_write 


contains

    subroutine ensemble_init(ens,name,nsim,dimname1,dim1)

        implicit none 

        type(ensemble_class) :: ens
        integer :: nsim, ndim1 
        character(len=*) :: name, dimname1
        double precision :: dim1(:) 

        ens%name     = trim(name) 
        ens%dimname1 = trim(dimname1)
        ens%nsim     = nsim 

        if (allocated(ens%dim1)) deallocate(ens%dim1)
        allocate(ens%dim1(size(dim1)))
        ens%dim1  = dim1 

        ! Allocate the metric variables
        ndim1 = size(ens%dim1)
        allocate(ens%var(nsim,ndim1))
        allocate(ens%obs(ndim1),ens%obs_err(ndim1))
        allocate(ens%err(nsim,ndim1))
        allocate(ens%err_std(nsim,ndim1))
        allocate(ens%wgt(nsim,ndim1))
        
        ! Preinitialize weights and errors
        ens%var     = 0.d0 
        ens%obs     = 0.d0 
        ens%obs_err = 1.d0 
        ens%wgt     = 1.d0 
        ens%err     = MISSING
        ens%err_std = MISSING 

        ! Allocate summary variables 
        allocate(ens%stat%obs(nsim))
        allocate(ens%stat%var(nsim))
        allocate(ens%stat%err(nsim))
        allocate(ens%stat%err_std(nsim))
        
        write(*,*) "Initialized ensemble: ",trim(ens%name),nsim,ndim1
        return 

    end subroutine ensemble_init

    subroutine ensemble_write(ens,filename,create)

        implicit none 

        type(ensemble_class) :: ens 
        character(len=*)     :: filename 
        logical, optional    :: create 
        logical              :: create_file 

        create_file = .TRUE. 
        if (present(create)) create_file = create 

        if (create_file) then 
            call nc_create(filename)
            call nc_write_dim(filename,"sim",x=1,dx=1,nx=ens%nsim)
            call nc_write_dim(filename,ens%dimname1,x=ens%dim1)
        end if 

        call nc_write(filename,"var",    ens%var,    dim1="sim",dim2=ens%dimname1,missing_value=MISSING)
        call nc_write(filename,"obs",    ens%obs,    dim1=ens%dimname1,missing_value=MISSING)
        call nc_write(filename,"obs_err",ens%obs_err,dim1=ens%dimname1,missing_value=MISSING)
        call nc_write(filename,"err",    ens%err,    dim1="sim",dim2=ens%dimname1,missing_value=MISSING)
        call nc_write(filename,"err_std",ens%err_std,dim1="sim",dim2=ens%dimname1,missing_value=MISSING)
        call nc_write(filename,"wgt",    ens%wgt,    dim1="sim",dim2=ens%dimname1,missing_value=MISSING)

        return 

    end subroutine ensemble_write


    subroutine sampler(samples,nsamples,pmin,pmax)

        implicit none 

        double precision  :: pmin, pmax 
        integer           :: nsamples, i  
        double precision, dimension(:), allocatable :: samples

        write(*,*) "Generating parameter samples: ", pmin, pmax, nsamples 

        if (allocated(samples)) deallocate(samples)
        allocate(samples(nsamples))

        samples(1) = pmin 
        do i = 1, nsamples
            samples(i) = pmin + (i-1)*(pmax-pmin)/dble(nsamples-1)
        end do 

        return 

    end subroutine sampler 

    subroutine sampler_latin_1D(samples,nsamples,pmin,pmax,seed,print)

        implicit none 

        integer           :: nsamples 
        double precision  :: pmin, pmax 
        double precision, dimension(:), allocatable :: samples
        double precision, dimension(:,:), allocatable :: samples_tmp
        integer, optional :: seed  
        logical, optional :: print 

        if (allocated(samples)) deallocate(samples)
        allocate(samples(nsamples))
        allocate(samples_tmp(1,nsamples))

        call sampler_latin_ND(samples_tmp,nsamples,[pmin],[pmax],seed)
        samples = samples_tmp(1,:)

        return 

    end subroutine sampler_latin_1D 

    subroutine sampler_latin_ND(samples,nsamples,pmin,pmax,seed,print)

        implicit none 

        integer            :: nsamples, npar, i  
        double precision   :: pmin(:), pmax(:) 
        double precision, dimension(:,:), allocatable :: samples
        integer, optional  :: seed 
        integer            :: seed_value  
        logical, optional  :: print
        logical            :: print_to_screen 
        integer, parameter :: print_max = 20 
        character(len=12)  :: pnames(20)
        npar = size(pmin)

        do i = 1, 20 
            if (i .lt. 10) then 
                write(pnames(i),"(a,i1)") "par0",i 
            else 
                write(pnames(i),"(a,i2)") "par",i 
            end if 
        end do 

        print_to_screen = .TRUE. 
        if (present(print)) print_to_screen = print 

        if (present(seed)) then 
            seed_value = seed 
        else 
            call get_seed ( seed_value )
            call random_initialize ( seed_value )
        end if 

        if (allocated(samples)) deallocate(samples)
        allocate(samples(npar,nsamples))

        call latin_random ( npar, nsamples, seed_value, samples )

        ! Scale samples up to parameter range of interest
        do i = 1, npar 
            samples(i,:) = samples(i,:)*(pmax(i)-pmin(i)) + pmin(i)
        end do  

        if (print_to_screen) then 
            write(*,*) "Latin-Hypercube samples (seed=",seed_value,")"
            write(*,"(5x,50a12)") pnames(1:npar)
            write(*,"(5x,50g12.4)") pmin 
            write(*,"(5x,50g12.4)") pmax
            write(*,"(a)") "-----"

            do i = 1, min(nsamples,print_max)
                write(*,"(i5,50g12.4)") i, samples(:,i)
            end do

            if (nsamples .gt. print_max) write(*,*) "..." 
        end if

        return 

    end subroutine sampler_latin_ND



    
    function stats_rmse(x,y,wgt,mask) result(err)

        implicit none 
        double precision, intent(IN) :: x(:), y(:)
        double precision, intent(IN), optional :: wgt 
        logical, intent(IN), optional :: mask(:) 
        double precision :: err, wgts(size(x))

        wgts = 1.d0 
        if (present(wgt)) wgts = wgt 
        if (present(mask)) then 
            where(.not. mask) wgts = 0.d0 
        end if 
        wgts = wgts/sum(wgts)

        if (present(mask)) then 
            err = dsqrt(sum(wgts*(x-y)**2,mask=mask))
        else 
            err = dsqrt(sum(wgts*(x-y)**2))
        end if 

        return 

    end function stats_rmse 


    function stats_stdev(x,y,wgt,mask) result(err)

        implicit none 
        double precision, intent(IN) :: x(:), y(:)
        double precision, intent(IN), optional :: wgt 
        logical, intent(IN), optional :: mask(:) 
        double precision :: err, wgts(size(x))

        wgts = 1.d0 
        if (present(wgt)) wgts = wgt 
        if (present(mask)) then 
            where(.not. mask) wgts = 0.d0 
        end if 
        wgts = wgts/sum(wgts)

        if (present(mask)) then 
            err = dsqrt(sum(wgts*(x-y)**2,mask=mask))
        else 
            err = dsqrt(sum(wgts*(x-y)**2))
        end if 

!         mean = SUM(X)/REAL(n_elts)
!         Std_Dev = SQRT(SUM((X-mean)**2)/REAL(n_elts))

        return 

    end function stats_stdev 

end module sim_sampler 


