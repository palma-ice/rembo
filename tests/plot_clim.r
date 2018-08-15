library(myr)
source("functions_plotting.r")



load_mar = function(filename)
{
    dat = my.read.nc(filename)
    dat$t2m   = dat$t3m + 273.15 
    dat$tsurf = dat$ts + 273.15 

    return(dat)
}

load_racmo = function(filename)
{
    dat = my.read.nc(filename)
    
    return(dat)
}

gen_masks = function(mask,z_srf,months=c(1:12))
{   # Input a 2D boolean mask, return 2D mask and 2D mask for each month

    mask[is.na(mask)] = FALSE 
    dim(mask) = dim(z_srf)
    mask12 = array(mask,dim=c(dim(mask),12))

    for (m in 1:12) {
        if (! m %in% months) mask12[,,m] = FALSE
    }

    return(list(mask=mask,mask12=mask12,ij=which(mask),ij12=which(mask12)))
}

# Load data 
if (TRUE) {

    domain   = "Greenland"
    grid_name = "GRL-20KM" 

    # domain   = "Antarctica"
    # grid_name = "ANT-40KM" 
    
    infldr  = file.path("../ice_data",domain,grid_name)
    outfldr = file.path("../output",grid_name)

    if (domain == "Greenland") {
        # Load MAR climatology 
        rcm = load_mar(file.path(infldr,paste0(grid_name,"_MARv3.5-ERA-30km-monthly_1981-2010.nc")))
        #rcm = load_mar(file.path(infldr,paste0(grid_name,"_MARv3.9-monthly-ERA_1981-2010.nc")))

        region_number = 3.2 

    } else if (domain == "Antarctica") {
        # Load RACMO climatology
        rcm = load_racmo(file.path(infldr,paste0(grid_name,"_RACMO23-ERAINT-HYBRID_1981-2010.nc")))

        region_number = 2.2
    }

    topo   = my.read.nc(file.path(infldr,paste0(grid_name,"_TOPO-RTOPO-2.0.1.nc")))
    
    basins  = my.read.nc(file.path(infldr,paste0(grid_name,"_BASINS-nasa.nc")))
    regions = my.read.nc(file.path(infldr,paste0(grid_name,"_REGIONS.nc")))
    
    mask  = topo$mask %in% c(2,3)   & abs(regions$mask-region_number) < 0.05 
    mask2 = topo$mask %in% c(1,2,3) & abs(regions$mask-region_number) < 0.05 
    masks = gen_masks(mask,topo$z_srf)

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

    #plot_year(rem,rcm,vnm="t2m",onm="t2m",long_name="2m-temp. (K)",top=topo,mask=mask2,type=ptype)  
    plot_comparison(rem,rcm,vnm="t2m",onm="t2m",long_name="2m-temp. (K)",top=topo,mask=mask2,alpha=0,type=ptype,months=c(1,7),zlim=c(-6,6))
    
    #plot_year(rem,rcm,vnm="al_s",onm="al",long_name="Surface albedo",top=topo,mask=mask2,type=ptype)  
    plot_comparison(rem,rcm,vnm="al_s",onm="al",long_name="Surface albedo",top=topo,mask=mask2,alpha=0,type=ptype,months=c(1,7),col=col.alb)
    
    # par(mfrow=c(2,3))
    # par(plt=c(0.05,0.9,0.05,0.95))

}

calc_albedo <- function(p,data,subset=c(1:length(data$t2m)),opt=FALSE)
{

    amin = p[1]
    amax = p[2]
    afac = p[3]
    tmid = p[4]

    #alb = amin+(amax-amin)*(0.5*tanh(afac*(data$t2m-tmid))+0.5) 
    
    alb = amin+(amax-amin)*afac*(data$t2m-tmid)

    if (opt) {
        err = rmse(alb-data$al_s,ii=subset)
        return(err)
    } else {
        return(alb)
    }
}

# Surface albedo 
if (FALSE) {

    #mm = gen_masks(topo$mask %in% c(2,3),topo$z_srf)   # Ice mask
    #mm = gen_masks(topo$mask %in% c(1,2,3),topo$z_srf)   # Ice and land mask

    #mm = gen_masks(rcm$msk >= 50,topo$z_srf)   # Ice mask (MAR)
    #mm = gen_masks(rcm$mask == 4,topo$z_srf,months=c(6:8))   # Ice and land mask (MAR)
    mm = gen_masks(rcm$mask == 4 & rcm$msk < 50,topo$z_srf,months=c(6,7,8))   # Land mask (MAR)
    
    
    test     = data.frame(t2m=rcm$t2m[mm$ij12],al_s=rcm$al[mm$ij12])
    new      = data.frame(t2m=seq(230,290,by=0.5))
    new$al_s = calc_albedo(p=c(0.2,0.83,-0.18,275.35),data=new)


    
    ## Optimize fits
    fit = optim(par=c(0.2,0.8,-0.18,275),calc_albedo,data=test,opt=TRUE)
    out = calc_albedo(fit$par,data=test)
    new$al_s = calc_albedo(fit$par,data=new)

    plot(test$t2m,test$al_s,pch=20,cex=0.7,col="grey60")
    lines(new$t2m,new$al_s,col=2,lwd=2)

}

