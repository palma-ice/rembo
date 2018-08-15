library(myr)

source("functions_plotting.r")



load_mar = function(filename)
{
    dat = my.read.nc(filename)
    dat$t2m   = dat$t3m + 273.15 
    dat$tsurf = dat$ts + 273.15 

    return(dat)
}

# Load data 
if (TRUE) {

    domain   = "Greenland"
    grid_name = "GRL-20KM" 

    infldr  = file.path("../ice_data",domain,grid_name)
    outfldr = file.path("../output",grid_name)

    # Load MAR climatology 
    #mar = load_mar(file.path(infldr,paste0(grid_name,"_MARv3.5-ERA-30km-monthly_1981-2010.nc")))
    mar = load_mar(file.path(infldr,paste0(grid_name,"_MARv3.9-monthly-ERA_1981-2010.nc")))
    
    topo   = my.read.nc(file.path(infldr,paste0(grid_name,"_TOPO-RTOPO-2.0.1.nc")))
    
}

# Load REMBO output
if (TRUE) {

    # Load REMBO output 
    rem = my.read.nc(file.path(outfldr,"rembo.nc"))

}

########################################
#
# PLOTTING and TESTING
#
########################################

mask = topo$mask %in% c(2,3)
# mask = marc$mask == 2 & topo$mask != 3 
# mask = marc$mask  > 0 & remboc$mask  > 0   # Ice sheet and land only  
# mask = marc$mask == 2 & remboc$mask == 2   # Ice sheet only 
mask.alpha = 0 
months = c(1,7,13)
ptype  = "png" 

col.axis = "grey40"
col.precip = c("white",jet.colors[2:length(jet.colors)])
col.smb    = rev(jet.colors)
col.melt   = c("white","yellow","red","darkred")
col.alb    = c("grey50","darkblue","skyblue","white")

# Plot comparison
if (TRUE) {

    plot_year(rem,mar,vnm="t2m",onm="t2m",long_name="2m-temp. (K)",top=topo,mask=mask,type=ptype)  
    plot_comparison(rem,mar,vnm="t2m",onm="t2m",long_name="2m-temp. (K)",top=topo,mask=mask,alpha=0,type=ptype,months=c(1,7),zlim=c(-6,6))
    
    
    # par(mfrow=c(2,3))
    # par(plt=c(0.05,0.9,0.05,0.95))

    # m = 1 
    # zlim = range(rem$t2m[,,m],mar$t2m[,,m])
    # image.plot(rem$xc,rem$yc,rem$t2m[,,m],axes=FALSE,ann=FALSE,zlim=zlim)
    # image.plot(mar$xc,mar$yc,mar$t2m[,,m],axes=FALSE,ann=FALSE,zlim=zlim)
    # image.plot(mar$xc,mar$yc,rem$t2m[,,m]-mar$t2m[,,m],axes=FALSE,ann=FALSE,zlim=c(-10,10))

    # m = 7 
    # zlim = range(rem$t2m[,,m],mar$t2m[,,m])
    # image.plot(rem$xc,rem$yc,rem$t2m[,,m],axes=FALSE,ann=FALSE,zlim=zlim)
    # image.plot(mar$xc,mar$yc,mar$t2m[,,m],axes=FALSE,ann=FALSE,zlim=zlim)
    # image.plot(mar$xc,mar$yc,rem$t2m[,,m]-mar$t2m[,,m],axes=FALSE,ann=FALSE,zlim=c(-10,10))


}