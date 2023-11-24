program test_monthlydaily
    
    use precision
    use monthlydaily 

    implicit none

    type(monthlydaily_class) :: mm

    call monthlydaily_init(mm,nmon_year=12,nday_month=30)

end program test_monthlydaily