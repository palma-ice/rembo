library(myr)

load_mar = function(filename)
{
    dat = my.read.nc(filename)
    dat$t2m   = dat$t3m + 273.15 
    dat$tsurf = dat$ts + 273.15 

    return(dat)
}

# Load data 
if (TRUE) {

    fldr = "../ice_data/Greenland"
    grid_name = "GRL-20KM"

    outfldr   = file.path("../output",grid_name)

    # Load MAR climatology 
    #mar = load_mar(file.path(fldr,grid_name,paste0(grid_name,"_MARv3.5-ERA-30km-monthly_1981-2010.nc")))
    mar = load_mar(file.path(fldr,grid_name,paste0(grid_name,"_MARv3.9-monthly-ERA_1981-2010.nc")))
    
    # Load REMBO output 
    rem = my.read.nc(file.path(outfldr,"rembo.nc"))

}


# Plot comparison
if (TRUE) {

    
    
    par(mfrow=c(2,3))
    par(plt=c(0.05,0.9,0.05,0.95))
    
    m = 1 
    zlim = range(rem$t2m[,,m],mar$t2m[,,m])
    image.plot(rem$xc,rem$yc,rem$t2m[,,m],axes=FALSE,ann=FALSE,zlim=zlim)
    image.plot(mar$xc,mar$yc,mar$t2m[,,m],axes=FALSE,ann=FALSE,zlim=zlim)
    image.plot(mar$xc,mar$yc,rem$t2m[,,m]-mar$t2m[,,m],axes=FALSE,ann=FALSE,zlim=c(-10,10))

    m = 7 
    zlim = range(rem$t2m[,,m],mar$t2m[,,m])
    image.plot(rem$xc,rem$yc,rem$t2m[,,m],axes=FALSE,ann=FALSE,zlim=zlim)
    image.plot(mar$xc,mar$yc,mar$t2m[,,m],axes=FALSE,ann=FALSE,zlim=zlim)
    image.plot(mar$xc,mar$yc,rem$t2m[,,m]-mar$t2m[,,m],axes=FALSE,ann=FALSE,zlim=c(-10,10))


}