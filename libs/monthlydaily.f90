module monthlydaily

    use precision, only : wp

    implicit none

    type monthlydaily_class
        integer,  allocatable :: day(:)
        integer,  allocatable :: month(:)
        integer,  allocatable :: day_of_month(:)
        integer,  allocatable :: m0(:)
        integer,  allocatable :: m1(:)
        real(wp), allocatable :: wt0(:)
        real(wp), allocatable :: wt1(:)

        integer :: nday_month
        integer :: nmon_year
        integer :: nday_year

    end type

    private
    public :: monthlydaily_class
    public :: monthlydaily_init 
    
contains

    subroutine monthlydaily_init(mm,nmon_year,nday_month)

        implicit none 

        type(monthlydaily_class), intent(OUT) :: mm
        integer, intent(IN) :: nmon_year
        integer, intent(IN) :: nday_month 

        mm%nmon_year  = nmon_year
        mm%nday_month = nday_month
        mm%nday_year  = mm%nmon_year*mm%nday_month

        ! Allocate object arrays to the right size
        call monthlydaily_alloc(mm)

        ! Calculate interpolation weights
        call interp_weights_monthly_to_daily(mm%m0,mm%m1,mm%wt0,mm%wt1,mm%nday_month,mm%nmon_year,verbose=.TRUE.)

        return

    end subroutine monthlydaily_init

    subroutine interp_weights_monthly_to_daily(m0,m1,wt0,wt1,nday_month,nmon_year,verbose)

        implicit none

        integer,  intent(out) :: m0(:)
        integer,  intent(out) :: m1(:)
        real(wp), intent(out) :: wt0(:)
        real(wp), intent(out) :: wt1(:)
        integer,  intent(IN)  :: nday_month
        integer,  intent(IN)  :: nmon_year 
        logical,  intent(IN)  :: verbose 

        ! Local variables
        integer :: dnow, q, k
        integer :: q0, q1, q2 
        integer :: mid
        
        ! ####################################################################
        ! Interpolate data in time: monthly => daily
        ! ####################################################################
        dnow = 0

        do q = 1, nmon_year
            do k = 1, nday_month

                ! Get current day of year
                dnow = dnow+1

                ! Determine indices of current, previous and next month

                q1 = q                     ! q1 is current month

                q0 = q1-1                  ! q0 is previous month
                if (q0 .eq. 0) q0 = 12     ! Loop to december

                if (q1 .eq. 13) q1 = 1     ! Loop to january for current month

                q2 = q1+1                  ! q2 is next month
                if (q2 .eq. 13) q2 = 1     ! Loop to january for next month

                ! Get Halfway point of current month (+/- 1 day)
                mid = nday_month / 2   

                ! Determine weights and month indices for interpolation

                if (k .lt. mid) then 
                    ! Weights between previous and current month

                    m0(dnow)  = q0
                    m1(dnow)  = q1
                    wt0(dnow) = dble(mid-k) / real(nday_month,wp)
                    wt1(dnow) = 1.0 - wt0(dnow)
                
                else if (k .gt. mid) then 
                    
                    m0(dnow)  = q1
                    m1(dnow)  = q2
                    wt0(dnow) = 1.0 - dble(k-mid) / real(nday_month,wp)
                    wt1(dnow) = 1.0 - wt0(dnow)
                
                else
                    ! Middle of the month
                    ! All weight given to current month

                    m0(dnow) = q1 
                    m1(dnow) = q1
                    wt0(dnow) = 1.0
                    wt1(dnow) = 0.0 

                end if

                if (verbose) then
                    write(*,*) q, k, m0(dnow), m1(dnow), wt0(dnow), wt1(dnow)
                end if

            end do
        end do
    
        return

    end subroutine interp_weights_monthly_to_daily

    subroutine monthlydaily_alloc(mm)

        implicit none

        type(monthlydaily_class), intent(INOUT) :: mm

        ! Local variables
        integer :: n 

        n = mm%nday_year

        ! First make sure everything is deallocated
        call monthlydaily_dealloc(mm)

        ! Now allocate to proper size
        allocate(mm%day(n))
        allocate(mm%month(n))
        allocate(mm%day_of_month(n))
        allocate(mm%m0(n))
        allocate(mm%m1(n))
        allocate(mm%wt0(n))
        allocate(mm%wt1(n))
        
        return

    end subroutine monthlydaily_alloc

    subroutine monthlydaily_dealloc(mm)

        implicit none

        type(monthlydaily_class), intent(INOUT) :: mm

        if (allocated(mm%day))              deallocate(mm%day)
        if (allocated(mm%month))            deallocate(mm%month)
        if (allocated(mm%day_of_month))     deallocate(mm%day_of_month)
        if (allocated(mm%m0))               deallocate(mm%m0)
        if (allocated(mm%m1))               deallocate(mm%m1)
        if (allocated(mm%wt0))              deallocate(mm%wt0)
        if (allocated(mm%wt1))              deallocate(mm%wt1)
        
        return

    end subroutine monthlydaily_dealloc
    
end module monthlydaily