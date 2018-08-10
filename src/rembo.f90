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
    
    subroutine rembo_init()

        implicit none 


        return 

    end subroutine rembo_init 

    subroutine rembo_end()

        implicit none 


        return 

    end subroutine rembo_end 
    
    subroutine rembo_alloc()

        implicit none 


        return 

    end subroutine rembo_alloc 
    
    subroutine rembo_dealloc()

        implicit none 


        return 

    end subroutine rembo_dealloc

end module rembo 