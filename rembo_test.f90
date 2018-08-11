program rembo_test

    use coord 
    use rembo 

    implicit none 

    type(grid_class)  :: grid 
    type(rembo_class) :: rembo1 

    call grid_init(grid,name="GRL-20KM",mtype="stereographic",units="kilometers", &
                                   lon180=.TRUE.,dx=20.d0,nx=91,dy=20.d0,ny=151, &
                                   lambda=-40.d0,phi=72.d0,alpha=8.4d0)

    call rembo_init(rembo1,par_path="par/Greenland.nml",domain="Greenland",grid=grid,year=1980)


    call rembo_print(rembo1,m=1,d=1,year=1980)

    write(*,*) "rembo_test.x finished succesfully."

end program rembo_test 

