

require(ncdf)

source("parameterizations.r")

lapse = 6.0e-3


load_icedata = function(domain,gridname)
{
    
    nm_topo = "TOPO-B13"
    if (domain=="Antarctica") {
        nm_topo = "TOPO-BEDMAP2"
    }

    dat = list()

    # Load topography data
    dat$topo = load_topo(paste0("EURICE/ice_data/",domain,"/",gridname,"/",gridname,"_",nm_topo,".nc"))
    # dat$topo1D = gen_1D(dat$topo)

    # Load basin information
    dat$topo$basin = load_basin(paste0("EURICE/ice_data/",domain,"/",gridname,"/",gridname,"_","BASINS-nasa",".nc"))

    # Load ceres data
    dat$ceres   = load_ceres_clim(paste0("EURICE/ice_data/",domain,"/",gridname,"/",gridname,"_CERES_2001-2013.nc"))
    # dat$ceres1D = gen_1D(dat$ceres)

    # Load ECMWF data of interest and get climatology
    dat$era   = load_ecmwf_clim(paste0("EURICE/ice_data/",domain,"/",gridname,"/ERA-INT/",gridname,"_ERA-INT_1981-2010.nc"))
    # dat$erac1D = gen_1D(dat$era)

    if (domain=="Antarctica") {
        dat$rcm   = load_racmo_clim(paste0("EURICE/ice_data/",domain,"/",gridname,"/",gridname,"_RACMO2-ANT3K55_HadCM3-A1B-monthly_2000-2010.nc"))
        # dat$rcm1D = gen_1D(rcm)
    } else {
        dat$rcm   =   load_mar_clim(paste0("EURICE/ice_data/",domain,"/",gridname,"/",gridname,"_MARv3.5-ERA-30km-monthly_1981-2010.nc"))
        # dat$rcm1D  = gen_1D(dat$rcm)
    }

    # Mask to define dataset (grl or ant)
    dat$mask = dat$topo$zs*0 
    
    dat$mask[] = 1
    if (domain=="Antarctica") dat$mask[] = 2

    return(dat)
}

load_fields <- function(filename,fields,sec_day=8.64e4)
{
    nc = open.ncdf(filename)

    # Load dimensions
    dat = list()
    for (q in 1:length(nc$dim)) {
        nm = names(nc$dim)[q]
        dat[[nm]] = get.var.ncdf(nc,nm)
    }

    if (! "time"  %in% names(dat)) dat[["time"]]  = 1 
    if (! "month" %in% names(dat)) dat[["month"]] = c(1:12)
    if (! "day"   %in% names(dat)) dat[["day"]]   = c(1:360)

    # Load static fields and fields of interest
    nms = c("lon2D","lat2D","x2D","y2D","zs","mask","basin",fields)
    for (q in 1:length(nms)) {
        nm = nms[q]
        if (nm %in% names(nc$var)) dat[[nm]] = get.var.ncdf(nc,nm)
    }

    close.ncdf(nc)

    if ("mask" %in% names(dat)) dat[["mask"]][is.na(dat[["mask"]])] = 0 
    if ("zs" %in% names(dat))   dat[["zs"]][is.na(dat[["zs"]])]     = 0 
    
    if ("sf" %in% fields & "rf" %in% fields) {
        dat[["pp"]] = dat[["sf"]] + dat[["rf"]]
    }

    if ("pp" %in% fields & "ps" %in% fields) {
        dat[["pr"]] = dat[["pp"]] - dat[["ps"]]
    }

    # Convert to mm/d as needed
    if (length(grep("MAR",filename))>0) sec_day = 1 
    if ("pp" %in% fields)      dat[["pp"]]      = dat[["pp"]]      *sec_day 
    if ("ps" %in% fields)      dat[["ps"]]      = dat[["ps"]]      *sec_day 
    if ("pr" %in% fields)      dat[["pr"]]      = dat[["pr"]]      *sec_day 
    if ("melt" %in% fields)    dat[["melt"]]    = dat[["melt"]]    *sec_day 
    if ("refr" %in% fields)    dat[["refr"]]    = dat[["refr"]]    *sec_day 
    if ("massbal" %in% fields) dat[["massbal"]] = dat[["massbal"]] *sec_day 
    if ("acc" %in% fields)     dat[["acc"]]     = dat[["acc"]]     *sec_day 
    if ("abl" %in% fields)     dat[["abl"]]     = dat[["abl"]]     *sec_day 

    # MAR only adjustments
    if (length(grep("MAR",filename))>0) {
        if ("t3m" %in% fields) dat[["t3m"]] = dat[["t3m"]] + 273.15    # C    => K
        if ("ts" %in% fields) dat[["ts"]] = dat[["ts"]]    + 273.15    # C    => K 
    }

    # Calculate air density
    dat[["rho_a"]] = calc_airdens(zs=dat[["zs"]])

    return(dat)
}


load_rembo <- function(filename,sec_day=8.64e4)
{

    dat = list(dataset=filename,scenario=NA,timestep=NA)

    # Open netcdf file
    nc = open.ncdf(filename)

    # Load dimensions
    dat$x = get.var.ncdf(nc,"xc")
    dat$y = get.var.ncdf(nc,"yc")
    dat$time  = get.var.ncdf(nc,"time")
    dat$month = get.var.ncdf(nc,"month")

    # Calculate dimension lengths
    nt = length(dat$time)
    nm = length(dat$month)
    nx = length(dat$x)
    ny = length(dat$y)

    # Load static fields
    dat$lon   = get.var.ncdf(nc,"lon2D")
    dat$lat   = get.var.ncdf(nc,"lat2D")
    dat$x2D   = get.var.ncdf(nc,"x2D")
    dat$y2D   = get.var.ncdf(nc,"y2D")
    dat$zs    = get.var.ncdf(nc,"zs")
    dat$mask  = get.var.ncdf(nc,"mask")
    dat$basin = get.var.ncdf(nc,"basin")

    vnms = names(nc$var)[! names(nc$var) %in% 
            c("xc","yc","time","month","lon2D","lat2D","x2D","y2D","zs","mask","basin")]
    # Climate
    # vnms = c("tt","ttcorr","pp","ps","q_s","q_sat","q_r","tcw","tcw_sat",
    #          "ccw","c_w","Z","ug","vg","uvg","ww","n","S","swd_toa","lwu_toa","ap","mnet",
    #          "t2m","swd","lwd","u","v","uv","sp")
    # vnms = c("t2m","pp","ps","q_s","q_r","tcw","tcw_sat",
    #          "ccw","Z","ug","vg","uvg","ww","cc","S","swd_toa","lwu_toa","alp",
    #          "swd","lwd","u","v","uv","sp")
    for (q in 1:length(vnms)) {
        vnm = vnms[q]
        dat[[vnm]] = get.var.ncdf(nc,vnm)
    }

    # # Surface 
    # vnms = c("tsurf","hsnow","alb","melt","refr","massbal","acc","lhf","shf")
    # for (q in 1:length(vnms)) {
    #     vnm = vnms[q]
    #     dat[[vnm]] = get.var.ncdf(nc,vnm)
    # }

    close.ncdf(nc) 

    ### Calculate additional fields
    
    # Rainfall
    dat$prrn  = dat$pr - dat$prsn

    # Get sea level temperature
    zs3D = array(dat$zs,dim=dim(dat$tas))
    dat$tsl = dat$tas + lapse*zs3D 

    return(dat)
}

load_basin <- function(filename)
{   # Return basin array 

    # Open netcdf file
    nc  = open.ncdf(filename)
    dat = get.var.ncdf(nc,"basin")
    close.ncdf(nc) 

    return(dat)
}

load_topo <- function(filename)
{

    dat = list(dataset=filename,scenario=NA,timestep=NA)

    # Open netcdf file
    nc = open.ncdf(filename)

    # Load dimensions
    dat$x = get.var.ncdf(nc,"xc")
    dat$y = get.var.ncdf(nc,"yc")

    # Load static fields
    dat$lon   = get.var.ncdf(nc,"lon2D")
    dat$lat   = get.var.ncdf(nc,"lat2D")
    dat$x2D   = get.var.ncdf(nc,"x2D")
    dat$y2D   = get.var.ncdf(nc,"y2D")
    dat$zs    = get.var.ncdf(nc,"zs")
    dat$zb    = get.var.ncdf(nc,"zb")
    dat$H     = get.var.ncdf(nc,"H")

    if ("mask" %in% names(nc$var)) {
        dat$mask  = get.var.ncdf(nc,"mask")
    } else {
        dat$mask  = get.var.ncdf(nc,"mask_ice")
    }

    # # For now, eliminate extra land (Baffin Island, Iceland, etc)
    # ii = which(dat$mask==3)
    # dat$zs[ii]   = 0 
    # dat$H[ii]    = 0 
    # dat$mask[ii] = 0 

    # Get gradient too 
    dzs = hgrad(dat$zs,dx=diff(dat$x[1:2])*1e3)
    dat$dzsdx  = dzs$x 
    dat$dzsdy  = dzs$y 
    dat$dzsdxy = dzs$xy 

    close.ncdf(nc) 

    return(dat)
}

load_ecmwf_clim <- function(filename)
{

    # Open netcdf file
    nc = open.ncdf(filename)

    dat = list(dataset=filename) 

    # Load dimensions
    dat$x = get.var.ncdf(nc,"xc")
    dat$y = get.var.ncdf(nc,"yc")
    dat$month = get.var.ncdf(nc,"month")
    dat$plev  = c(1000,950,850,750,700,650,600,550,500)

    # Load static fields
    dat$lon  = get.var.ncdf(nc,"lon2D")
    dat$lat  = get.var.ncdf(nc,"lat2D")
    dat$x2D  = get.var.ncdf(nc,"x2D")
    dat$y2D  = get.var.ncdf(nc,"y2D")

    dat$zs   = get.var.ncdf(nc,"zs") #/ 9.80665
    dat$mask = dat$zs*0 
    dat$mask[dat$zs>100] = 1 

    # Get gradient too 
    dzs = hgrad(dat$zs,dx=diff(dat$x[1:2])*1e3)
    dat$dzsdx  = dzs$x 
    dat$dzsdy  = dzs$y 
    dat$dzsdxy = dzs$xy 

    # Load the surface fields 
    vnms = c("sp","tcw","tcc","ciw","clw","u10","v10","t2m","al","sst")
    for (q in 1:length(vnms)) {
        vnm = vnms[q]
        dat[[vnm]] = get.var.ncdf(nc,vnm)
    }

    close.ncdf(nc) 

    vnms = c("p_t","p_z","p_u","p_v","p_w")
    for (q0 in 1:length(dat$plev)) {
        plev = dat$plev[q0]

        # Make new filename and load pressure level data 
        filename2 = gsub("ERA-INT_",paste0("ERA-INT-",plev,"Mb_"),filename)
        nc = open.ncdf(filename2)
        
        for (q in 1:length(vnms)) {
            vnm = vnms[q]
            if (q0==1) dat[[vnm]] = array(NA,dim=c(length(dat$plev),dim(dat$t2m)))
            dat[[vnm]][q0,,,] = get.var.ncdf(nc,vnm)
        }

        close.ncdf(nc) 

    }
        
    ### Calculate additional fields
    
    # Get geopotential height 
    dat$p_Z    = dat$p_z / 9.80655

    # Get the 750Mb level specifically
    l = which(dat$plev==750)
    dat$p750_t = dat$p_t[l,,,]
    dat$p750_z = dat$p_z[l,,,]
    dat$p750_u = dat$p_u[l,,,]
    dat$p750_v = dat$p_v[l,,,]
    dat$p750_w = dat$p_w[l,,,]
    dat$p750_Z = dat$p_Z[l,,,]

    # Get the 1000Mb level for 'sea level temp.'
    l = which(dat$plev==1000)
    dat$p1000_t = dat$p_t[l,,,]

    # Get the sea level temp 'T_a' from 500Mb level 
    l = which(dat$plev==500)
    dat$p500_t = dat$p_t[l,,,]
    zs = array(dat$zs,dim=dim(dat$p500_t))
    z500 = calc_elev(500)
    dat$ta = dat$p500_t + lapse*(z500-zs)

    # Magnitude of horizontal wind field
    dat$uv10     = sqrt(dat$u10^2 + dat$v10^2)
    dat$p750_uv  = sqrt(dat$p750_u^2 + dat$p750_v^2)

    # Convert omega to vertical velocity in observations for comparison 
    dat$p750_omega = dat$p750_w 
    dat$p750_w = conv_omega_w(dat$p750_omega,T=dat$p750_t,p=750*100)

    # Get cloud content water
    dat$ccw = dat$ciw + dat$clw 

    # Calculate saturated specific humidity quantities too 
    
    # Calculate air density
    dat$rho_a = calc_airdens(zs=dat$zs)

    # Calculate surface temperature without inversion based on interpolation of 
    # higher pressure layers and lapse rate
    zs = array(dat$zs,dim=dim(dat$t2m))
    dat$gamma = get_gamma(dat$plev,p_t=dat$p_t,p1=650,p2=600)
    dat$tsl0  = get_tsl(dat$plev,p_t=dat$p_t,p1=650,p2=600)
    dat$t2m0  = get_t2m(plev=dat$plev,p_t=dat$p_t,p_Z=dat$p_Z,zs=zs,p1=650,p2=600)
    dat$tcorr = dat$t2m0 - dat$t2m 
    dat$tsl   = dat$t2m  + lapse*zs 

    dat$q_sat   = dat$t2m*NA 
    dat$tcw_sat = dat$t2m*NA 

    for (m in 1:dim(dat$t2m)[3]) {

        # Calculate saturated water content 
        dat$q_sat[,,m] = calc_qsat(Ts=dat$t2m0[,,m],zs=dat$zs)

        # Calculate saturated total water content
        dat$tcw_sat[,,m] = dat$q_sat[,,m] * ( dat$rho_a *2.5e3 )

    }


    # Get relative humidity 
    dat$q_r = dat$tcw / dat$tcw_sat 

    return(dat)
}

load_racmo_clim <- function(filename)
{

    dat = list(dataset=filename)

    # Open netcdf file
    nc = open.ncdf(filename)

    # Load dimensions
    dat$x = get.var.ncdf(nc,"xc")
    dat$y = get.var.ncdf(nc,"yc")
    dat$month = get.var.ncdf(nc,"month")

    # Calculate dimension lengths
    nm = length(dat$month)
    nx = length(dat$x)
    ny = length(dat$y)

    # Load static fields
    dat$lon  = get.var.ncdf(nc,"lon2D")
    dat$lat  = get.var.ncdf(nc,"lat2D")
    dat$x2D  = get.var.ncdf(nc,"x2D")
    dat$y2D  = get.var.ncdf(nc,"y2D")
    dat$zs   = get.var.ncdf(nc,"zs") / 9.81

    dat$mask = get.var.ncdf(nc,"mask_ice")
    dat$mask[dat$mask==1] = 2   # ice 

    # Get gradient too 
    dzs = hgrad(dat$zs,dx=diff(dat$x[1:2])*1e3)
    dat$dzsdx  = dzs$x 
    dat$dzsdy  = dzs$y 
    dat$dzsdxy = dzs$xy 

    # Climate
    vnms = c("t2m","cc","evap","pr","qs","rf","rz","ru","smb","sf","me","subl","ts","uas",
             "vas","lhf","lwd","lwn","swd","swn","al") 
    for (q in 1:length(vnms)) {
        vnm = vnms[q]
        if (vnm %in% names(nc$var)) {
            dat[[vnm]] = get.var.ncdf(nc,vnm)
        }
    }

    close.ncdf(nc) 

    # Change some names 
    dat$q_s = dat$qs 
    dat$qs  = NULL 

    ### Calculate additional fields
    
    # Limit melt field to positive values 
    dat$me[dat$me<0] = NA 

    # Magnitude of horizontal wind field
    dat$uvas  = sqrt(dat$uas^2 + dat$vas^2)

    # Calculate air density
    dat$rho_a = calc_airdens(zs=dat$zs)

    dat$q_sat   = dat$t2m*NA 
    dat$tcw     = dat$t2m*NA 
    dat$tcw_sat = dat$t2m*NA 

    for (m in 1:dim(dat$t2m)[3]) {

        # Calculate saturated water content 
        dat$q_sat[,,m] = calc_qsat(Ts=dat$t2m[,,m],zs=dat$zs)

        # Calculate total water content
        dat$tcw[,,m]     = dat$q_s[,,m]   * ( dat$rho_a *2.5e3 )
        dat$tcw_sat[,,m] = dat$q_sat[,,m] * ( dat$rho_a *2.5e3 )

    }
        
    # Get relative humidity 
    dat$q_r = dat$tcw / dat$tcw_sat 
    
    # Get sea level temperature
    zs3D    = array(dat$zs,dim=dim(dat$t2m))
    dat$tsl = dat$t2m + lapse*zs3D 
    
    return(dat)
}

load_mar_clim <- function(filename)
{

    dat = list(dataset=filename)

    # Open netcdf file
    nc = open.ncdf(filename)

    # Load dimensions
    dat$x = get.var.ncdf(nc,"xc")
    dat$y = get.var.ncdf(nc,"yc")
    dat$month = get.var.ncdf(nc,"month")

    # Calculate dimension lengths
    nm = length(dat$month)
    nx = length(dat$x)
    ny = length(dat$y)

    # Load static fields
    dat$lon  = get.var.ncdf(nc,"lon2D")
    dat$lat  = get.var.ncdf(nc,"lat2D")
    dat$x2D  = get.var.ncdf(nc,"x2D")
    dat$y2D  = get.var.ncdf(nc,"y2D")
    dat$zs   = get.var.ncdf(nc,"zs")

    if ("msk" %in% names(nc$var)) {
        dat$mask0 = get.var.ncdf(nc,"mask")
        dat$msk   = get.var.ncdf(nc,"msk")

        dat$mask = dat$mask0*0          # Ocean
        dat$mask[dat$mask0==4] = 1      # Land
        dat$mask[dat$msk > 50] = 2      # Ice
        dat$mask[dat$mask0==2] = -1     # Sea ice
    } else {
        dat$mask = get.var.ncdf(nc,"mask")
    }

    # Get gradient too 
    dzs = hgrad(dat$zs,dx=diff(dat$x[1:2])*1e3)
    dat$dzsdx  = dzs$x 
    dat$dzsdy  = dzs$y 
    dat$dzsdxy = dzs$xy 

    # cutline = 0.7*dat$x2D - 575 
    # cutmask = dat$y2D-cutline
    # cutmask[cutmask>=0] = NA 
    # cutmask[cutmask<0]  = 1

    # Climate
    vnms = c("smb","ru","rz","me","ts","t3m","sf","rf","su","al","pdd","t3m_min","t3m_max",
             "swd","lwd","lwu","shf","lhf","sp","u","v","q","cc","SH0","SH2","SH3") 
    for (q in 1:length(vnms)) {
        vnm = vnms[q]
        if (vnm %in% names(nc$var)) {
            dat[[vnm]] = get.var.ncdf(nc,vnm)

            # # Delete regions not of interest
            # cutmask3D = array(cutmask,dim=dim(dat[[vnm]]))
            # dat[[vnm]] = dat[[vnm]]*cutmask3D
        }
    }

    close.ncdf(nc) 

    # Rename some fields for consistency 
    dat$t2m     = dat$t3m 
    dat$t2m_min = dat$t3m_min 
    dat$t2m_max = dat$t3m_max 
    dat$t3m = dat$t3m_min = dat$t3m_max = NULL 
    
    dat$uas     = dat$u 
    dat$vas     = dat$v 
    dat$u = dat$v = NULL 

    dat$q_s     = dat$q 
    dat$q       = NULL 

    ### Calculate additional fields
    
    # Limit melt field to positive values 
    dat$me[dat$me<0] = NA 

    # Magnitude of horizontal wind field
    dat$uvas  = sqrt(dat$uas^2 + dat$vas^2)

    # Total precip 
    dat$pr  = dat$sf + dat$rf 

    # Convert units
    dat$q_s     = dat$q_s * 1e-3       # g/kg => kg/kg
    dat$sp      = dat$sp  * 1e2        # hPa  => Pa
    dat$t2m     = dat$t2m + 273.15    # C    => K
    dat$t2m_max = dat$t2m_max + 273.15    # C    => K
    dat$t2m_min = dat$t2m_min + 273.15    # C    => K
    dat$ts      = dat$ts  + 273.15            # C    => K 

    dat$t2md = dat$t2m_max - dat$t2m_min 

    # Calculate saturated specific humidity quantities too 
    
    # Calculate air density
    dat$rho_a = calc_airdens(zs=dat$zs)

    dat$q_sat   = dat$t2m*NA 
    dat$tcw     = dat$t2m*NA 
    dat$tcw_sat = dat$t2m*NA 

    for (m in 1:dim(dat$t2m)[3]) {

        # Calculate saturated water content 
        dat$q_sat[,,m] = calc_qsat(Ts=dat$t2m[,,m],zs=dat$zs)

        # Calculate total water content
        dat$tcw[,,m]     = dat$q_s[,,m]   * ( dat$rho_a *2.5e3 )
        dat$tcw_sat[,,m] = dat$q_sat[,,m] * ( dat$rho_a *2.5e3 )

    }
        
    # Get relative humidity 
    dat$q_r = dat$tcw / dat$tcw_sat 
    
    # Get sea level temperature
    zs3D    = array(dat$zs,dim=dim(dat$t2m))
    dat$tsl = dat$t2m + lapse*zs3D 
    
    # Diagnose sensible heat exchange coefficient
    dat$cm = array(NA,dim=dim(dat$shf))
    for (m in 1:12) dat$cm[,,m] = 1e3*dat$shf[,,m] / (dat$rho_a*1e3*dat$uv[,,m]*(dat$ts[,,m]-dat$t2m[,,m]))

    return(dat)
}

calc_hsnow_equil <- function(dh,hmax=10,ny=60)
{
    dh[is.na(dh)] = 0 
    dh_ave = apply(dh,c(1,2),sum)

    hsnow = array(0,dim=dim(dh))
    hsnow[,,1] = hmax 
    hsnow[,,1][dh_ave<0] = 0 

    for (y in 1:ny) {
        for (k in 1:dim(dh)[3]) {
            k0 = k-1
            if (k0 == 0) k0 = 12 

            hsnow[,,k] = hsnow[,,k0] + dh[,,k]*1e-3*30
            hsnow[hsnow<1e-3] = 0 
            hsnow[hsnow>hmax] = hmax 
        }
    }

    return(hsnow)
}


load_mar_clim_daily <- function(filename)
{

    dat = list(dataset=filename)

    # Open netcdf file
    nc = open.ncdf(filename)

    # Load dimensions
    dat$x = get.var.ncdf(nc,"X10_692_59")
    dat$y = get.var.ncdf(nc,"Y21_1322_109")
    dat$time = get.var.ncdf(nc,"time")
    dat$day  = c(1:length(dat$time))

    # Calculate dimension lengths
    nt = length(dat$day)
    nx = length(dat$x)
    ny = length(dat$y)

    # Load static fields
    dat$lon  = get.var.ncdf(nc,"LON")
    dat$lat  = get.var.ncdf(nc,"LAT")
    dat$zs   = get.var.ncdf(nc,"SH")
    dat$mask = get.var.ncdf(nc,"SOL")

    # Get gradient too 
    dzs = hgrad(dat$zs,dx=diff(dat$x[1:2])*1e3)
    dat$dzsdx  = dzs$x 
    dat$dzsdy  = dzs$y 
    dat$dzsdxy = dzs$xy 

    # Climate
    vnms  = c("SMB","RU","ME","STT","TT","SF","RF","SU","AL1",
              "SWD","LWD","SHF","LHF","SP","UU","VV","QQ") 
    vnms1 = c("smb","ru","me","ts","t3m","sf","rf","su","al",
              "swd","lwd","shf","lhf","sp","u","v","q") 
    for (q in 1:length(vnms)) {
        vnm  = vnms[q]
        vnm1 = vnms1[q]
        dat[[vnm1]] = get.var.ncdf(nc,vnm)
        cat("Loaded: ",vnm,"\n")
    }

    close.ncdf(nc) 

    ### Calculate additional fields
    
    # Magnitude of horizontal wind field
    dat$uv  = sqrt(dat$u^2 + dat$v^2)

    # Total precip 
    dat$pp  = dat$sf + dat$rf 

    # Convert units
    dat$q_s = dat$q  * 1e-3       # g/kg => kg/kg
    dat$sp  = dat$sp * 1e2        # hPa  => Pa
    dat$t3m = dat$t3m + 273.15    # C    => K
    dat$ts  = dat$ts  + 273.15    # C    => K 

    # Calculate saturated specific humidity quantities too 
    
    # Calculate air density
    dat$rho_a = calc_airdens(zs=dat$zs)

    dat$q_sat   = dat$t3m*NA 
    dat$tcw     = dat$t3m*NA 
    dat$tcw_sat = dat$t3m*NA 

    for (m in 1:dim(dat$t3m)[3]) {

        # Calculate saturated water content 
        dat$q_sat[,,m] = calc_qsat(Ts=dat$t3m[,,m],zs=dat$zs)

        # Calculate total water content
        dat$tcw[,,m]     = dat$q_s[,,m]   * ( dat$rho_a *2.5e3 )
        dat$tcw_sat[,,m] = dat$q_sat[,,m] * ( dat$rho_a *2.5e3 )

    }
        
    # Get relative humidity 
    dat$q_r = dat$tcw / dat$tcw_sat 
    
    # Get sea level temperature
    zs3D = array(dat$zs,dim=dim(dat$t3m))
    dat$tsl = dat$t3m + lapse*zs3D 
    
    # Diagnose sensible heat exchange coefficient
    dat$cm = array(NA,dim=dim(dat$shf))
    for (m in 1:nt) dat$cm[,,m] = 1e3*dat$shf[,,m] / (dat$rho_a*1e3*dat$uv[,,m]*(dat$ts[,,m]-dat$t3m[,,m]))

    return(dat)
}

load_mar_clean <- function(filename)
{   # Load a "clean" mar dataset, produced with the fortran program

    dat = list(dataset=filename)

    # Open netcdf file
    nc = open.ncdf(filename)

    # Load dimensions
    dat$x = get.var.ncdf(nc,"x")
    dat$y = get.var.ncdf(nc,"y")
    dat$time = get.var.ncdf(nc,"time")
    dat$month = c(1:12)
    dat$day   = c(1:365)
    if (length(dat$time)>12) dat$day  = c(1:length(dat$time))

    # Calculate dimension lengths
    nt = length(dat$time)
    nx = length(dat$x)
    ny = length(dat$y)

    # Load static fields
    dat$lon  = get.var.ncdf(nc,"LON")
    dat$lat  = get.var.ncdf(nc,"LAT")
    dat$zs   = get.var.ncdf(nc,"SH")
    dat$srf  = get.var.ncdf(nc,"SRF")
    dat$msk  = get.var.ncdf(nc,"MSK")

    dat$mask = dat$srf*0
    dat$mask[dat$srf%in%c(1,2)] = 0
    dat$mask[dat$srf==4] = 1
    dat$mask[dat$msk>50] = 2 

    # Get gradient too 
    dzs = hgrad(dat$zs,dx=diff(dat$x[1:2])*1e3)
    dat$dzsdx  = dzs$x 
    dat$dzsdy  = dzs$y 
    dat$dzsdxy = dzs$xy 

    # Climate
    vnms  = c("SHSN0","SHSN2","SHSN3","SMB","SU","ME","RZ",
              "SF","RF","RU","UU","VV","TT","QQ","SP","RH","TTMIN","TTMAX",
              "UV","SWD","LWD","LWU","SHF","LHF","AL","CC","ST","PDD")

    vnms1 = c("shsn0","shsn2","shsn3","smb","su","me","rz",
              "sf","rf","ru","u","v","t3m","q","sp","rh","t3m_min","t3m_max",
              "uv","swd","lwd","lwu","shf","lhf","al","cc","ts","pdd")
    for (q in 1:length(vnms)) {
        vnm  = vnms[q]
        vnm1 = vnms1[q]
        dat[[vnm1]] = get.var.ncdf(nc,vnm)
        cat("Loaded: ",vnm,"\n")
    }

    close.ncdf(nc) 

    ### Calculate additional fields
    
    # Total precip 
    dat$pp  = dat$sf + dat$rf 

    # Convert units
    dat$q_s = dat$q  * 1e-3       # g/kg => kg/kg
    dat$sp  = dat$sp * 1e2        # hPa  => Pa
    dat$t3m = dat$t3m + 273.15    # C    => K
    dat$t3m_max = dat$t3m_max + 273.15    # C    => K
    dat$t3m_min = dat$t3m_min + 273.15    # C    => K
    dat$ts  = dat$ts  + 273.15    # C    => K 

    # Calculate saturated specific humidity quantities too 
    
    # Calculate air density
    dat$rho_a = calc_airdens(zs=dat$zs)

    dat$q_sat   = dat$t3m*NA 
    dat$tcw     = dat$t3m*NA 
    dat$tcw_sat = dat$t3m*NA 

    for (m in 1:dim(dat$t3m)[3]) {

        # Calculate saturated water content 
        dat$q_sat[,,m] = calc_qsat(Ts=dat$t3m[,,m],zs=dat$zs)

        # Calculate total water content
        dat$tcw[,,m]     = dat$q_s[,,m]   * ( dat$rho_a *2.5e3 )
        dat$tcw_sat[,,m] = dat$q_sat[,,m] * ( dat$rho_a *2.5e3 )

    }
        
    # Get relative humidity 
    dat$q_r = dat$tcw / dat$tcw_sat 
    
    # Get sea level temperature
    zs3D = array(dat$zs,dim=dim(dat$t3m))
    dat$tsl = dat$t3m + lapse*zs3D 
    
    # Diagnose sensible heat exchange coefficient
    dat$cm = array(NA,dim=dim(dat$shf))
    for (m in 1:nt) dat$cm[,,m] = 1e3*dat$shf[,,m] / (dat$rho_a*1e3*dat$uv[,,m]*(dat$ts[,,m]-dat$t3m[,,m]))

    return(dat)
}

load_mar4D <- function(filename)
{

    dat = list(dataset=filename,scenario=NA,timestep=NA)

    # Open netcdf file
    nc = open.ncdf(filename)

    # Load dimensions
    dat$x = get.var.ncdf(nc,"xc")
    dat$y = get.var.ncdf(nc,"yc")
    dat$time  = get.var.ncdf(nc,"time")
    dat$month = get.var.ncdf(nc,"month")

    # Calculate dimension lengths
    nt = length(dat$time)
    nm = length(dat$month)
    nx = length(dat$x)
    ny = length(dat$y)

    # Load static fields
    dat$lon  = get.var.ncdf(nc,"lon2D")
    dat$lat  = get.var.ncdf(nc,"lat2D")
    dat$x2D  = get.var.ncdf(nc,"x2D")
    dat$y2D  = get.var.ncdf(nc,"y2D")
    dat$zs   = get.var.ncdf(nc,"zs")
    dat$mask = get.var.ncdf(nc,"mask")

    # Climate
    vnms = c("smb","ru","me","ts","t3m","sf","rf","su","al",
             "swd","lwd","shf","lhf","sp","u","v","q","cc","SH3")
    for (q in 1:length(vnms)) {
        vnm = vnms[q]
        dat[[vnm]] = get.var.ncdf(nc,vnm)
    }

    close.ncdf(nc) 

    ### Calculate additional fields
    
    # Magnitude of horizontal wind field
    dat$uv  = sqrt(dat$u^2 + dat$v^2)

    # Total precip 
    dat$pp  = dat$sf + dat$rf 

    # Convert units
    dat$q_s = dat$Q  * 1e-3       # g/kg => kg/kg
    dat$sp  = dat$sp * 1e2        # hPa  => Pa
    dat$t3m = dat$t3m + 273.15    # C    => K
    dat$ts  = dat$ts  + 273.15    # C    => K 

    # Calculate saturated specific humidity quantities too 
    
    # Calculate air density
    dat$rho_a = calc_airdens(zs=dat$zs)

    # dat$q_sat   = dat$t3m*NA 
    # dat$tcw     = dat$t3m*NA 
    # dat$tcw_sat = dat$t3m*NA 

    # for (k in 1:dim(dat$t3m)[4]) {
    #     for (m in 1:dim(dat$t3m)[3]) {

    #         # Calculate saturated water content 
    #         dat$q_sat[,,m,k] = calc_qsat(Ts=dat$t3m[,,m,k],zs=dat$zs)

    #         # Calculate total water content
    #         dat$tcw[,,m,k]     = dat$q_s[,,m,k]   * ( dat$rho_a *2.5e3 )
    #         dat$tcw_sat[,,m,k] = dat$q_sat[,,m,k] * ( dat$rho_a *2.5e3 )

    #     }
    # }
        
    # # Get relative humidity 
    # dat$q_r = dat$tcw / dat$tcw_sat 
    
    # Get sea level temperature
    zs4D = array(dat$zs,dim=dim(dat$t3m))
    dat$tsl = dat$t3m + lapse*zs4D 
 
    return(dat)
}

# Load the 3D (climatological) CERES dataset (nx,ny,nmonth)
load_ceres_clim <- function(filename)
{
    # Open netcdf file
    nc = open.ncdf(filename)

    dat = list() 

    # Load dimensions
    dat$x = get.var.ncdf(nc,"xc")
    dat$y = get.var.ncdf(nc,"yc")
    dat$month = get.var.ncdf(nc,"month")

    # Load static fields
    dat$lon  = get.var.ncdf(nc,"lon2D")
    dat$lat  = get.var.ncdf(nc,"lat2D")
    dat$x2D  = get.var.ncdf(nc,"x2D")
    dat$y2D  = get.var.ncdf(nc,"y2D")

    # Load variable 
    dat$sw_all   = get.var.ncdf(nc,"toa_sw_all")
    dat$sw_clr   = get.var.ncdf(nc,"toa_sw_clr")
    dat$lw_all   = get.var.ncdf(nc,"toa_lw_all")
    dat$lw_clr   = get.var.ncdf(nc,"toa_lw_clr")
    dat$net_all  = get.var.ncdf(nc,"toa_net_all")
    dat$net_clr  = get.var.ncdf(nc,"toa_net_clr")

    dat$S        = get.var.ncdf(nc,"solar")

    # Close the netcdf file
    close.ncdf(nc) 

    # Calculate additional fields
    dat$Snet = dat$S-dat$sw_all
    dat$Snet[dat$Snet<0] = 0 

    sw_all = dat$sw_all 
    ii = which(sw_all > dat$S)
    sw_all[ii] = dat$S[ii]

    dat$ap   = sw_all / dat$S 
    dat$ap[!is.finite(dat$ap)] = NA 
    dat$ap[dat$S < 5] = NA 

    return(dat)
}


######################
#
#  Extra functions 
#
######################

# Get climatology
gen_clim <- function(dat,FUN)
{
    datc    = dat[c("x","y","month","plev","ph","lon","lat","x2D","y2D","zs","mask","rho_a")]
    nms = names(dat)
    for (i in 1:length(nms)) {
        nm = nms[i]
        ndim = length(dim(dat[[nm]]))
        if (ndim >= 4) {
            datc[[nm]] = apply(dat[[nm]],MARGIN=c(1:(ndim-1)),FUN=FUN,na.rm=TRUE)
        } 
    }

    return(datc)
}

# Get 1D dataset 
gen_1D <- function(dat)
{

    xnm = "x2D"
    if (! xnm %in% names(dat)) xnm = "xx"
    if (! xnm %in% names(dat)) xnm = "zs"
    if (! xnm %in% names(dat)) cat("gen_1D:: error: 2D variable not found: ",xnm,"\n")
    
    nx = dim(dat[[xnm]])[1]
    ny = dim(dat[[xnm]])[2]
    nm = 12 

    dat1D = data.frame(index=c(1:(nx*ny*nm)))

    allnms = names(dat)
    for (q in 1:length(allnms)) {

        vnm = allnms[q]


        if (length(dat[[vnm]]) %in% c(nx,nx*ny,nx*ny*nm)) {
            dat1D[[vnm]] = as.vector(array(dat[[vnm]],dim=c(nx,ny,nm)))
        } else if (length(dat[[vnm]] %in% c(ny))) {
            tmp = t(array(dat[[vnm]],dim=c(ny,nx)))
            dat1D[[vnm]] = as.vector(array(tmp,dim=c(nx,ny,nm)))
        } else if (length(dat[[vnm]] %in% c(nm))) {
            tmp12 = array(NA,dim=c(nx,ny,nm))
            for (m in 1:nm) tmp12[,,m] = dat[[vnm]][m]
            dat1D[[vnm]] = as.vector(array(tmp12,dim=c(nx,ny,nm)))
        }
    }

    return(dat1D)
}

gen_1D_combined <- function(grl,ant)
{
    grl1D = gen_1D(grl)
    ant1D = gen_1D(ant)

    grl1D$domain = 1 
    ant1D$domain = 2 

    dat1D = rbind(grl1D,ant1D)
    return(dat1D)
}


