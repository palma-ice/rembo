
g0      = 9.80665
sec_day = 8.64e4
omega   = 7.2921e-5

# Get the vector magnitude from two (or three) components
calc_magnitude <- function(u,v,w=0) return(sqrt(u^2+v^2+w^2))

# Get coriolis parameter (1/s)
calc_coriolis <- function(lat,omega=7.2921e-5) return(2*omega*sin(lat*pi/180))

# Convert geopotential into geopotential height, (m2/s2) => (m)
calc_geo_height <- function(Phi,g0=9.80665) return(Phi/g0)

# Get horizontal geostrophic wind component, u
calc_u_geo <- function(dZdy,f,g0=9.80665,dx=50e3)
{   
    ug = -(g0 / f) * dZdy 
    return(ug)
}

# Get horizontal geostrophic wind component, v
calc_v_geo <- function(dZdx,f,g0=9.80665,dx=50e3)
{   
    vg =  (g0 / f) * dZdx 
    return(vg)
}

calc_w <- function(u,v,dzdx,dzdy)
{
    w = u*dzdx + v*dzdy 
    return(w)
}

conv_omega_w <- function(omega,T,p,rgas=287.058,g0=9.80665)
{
    rgas = 287.058            # J/(kg-K) => m2/(s2 K)
    g    = 9.80665            # m/s2
    rho  = p/(rgas*T)         # density => kg/m3
    
    w    = -1*omega/(rho*g0)     

    return(w)
}

# calc_qsat <- function()
# {
#     qsat = 3.8e-3*exp(17.67*T/(T+243.5))
#     return(qsat)
# }

# Surface pressure ( Pa=> kg/(m s2) )
calc_spa    <- function(zs,p0=101325,g=9.80665,M=0.0289644,R=8.31447,T0=298.15) p0*exp(-g*M*zs/(R*T0))
# Air density (kg/m3) for given elevation
calc_airdens <- function(zs) 1.3 * exp(-zs/8.6e3)

calc_elev <- function(p,p0=1013.25,g=9.80665,M=0.0289644,R=8.31447,T0=288.15) -log(p/p0)*R*T0/(g*M)

calc_es <- function(T,e0=610.94,b=17.625,T0=273.15,T1=30.11) 
{   # August-Roche-Magnus approximation of the Clausius-Clapeyron equation
    # for equilibrium vapor pressure.
    #   Lawrence, M. G.: The Relationship between Relative Humidity and the 
    #   Dewpoint Temperature in Moist Air: A Simple Conversion 
    #   and Applications, Bull. Am. Meteorol. Soc., 86(2), 225–233, 
    #   doi:10.1175/BAMS-86-2-225, 2005.
    # Default parameter values for saturation vapor pressure from: 
    #   Alduchov, O. and Eskridge, R.: Improved Magnus form approximation
    #   of saturation vapor pressure, J. Appl. Meteorol., 35, 601–609, 1996.
    e_s = e0 * exp(b * (T-T0) / (T-T1))
    return(e_s)
}

calc_qs <- function(par,Ts,zs,eps=0.622)
{   # Specific humidity is ratio of the equilibrium vapor pressure 
    # with the pressure of the remaining gases in the air (scaled by molecular weights)
    # eps: ratio of the molecular weights of water and dry air
    psurf = calc_spa(zs)   # Surface pressure, Pa 
    e_s   = calc_es(T=Ts,e0=par[1],b=par[2])   # Effective saturation vapor pressure, Pa
    q_s   = eps * e_s / (psurf - (1-eps)*e_s)

    if (length(par)==3) q_s = q_s + par[3]

    return(q_s)
}

# Saturated specific humidity
calc_qsat <- function(Ts,zs) calc_qs(par=c(610.94,17.625),Ts,zs)


opt_tcw <- function(p,Ts,zs,rho_a,tcw_obs)
{
    tcw = calc_qs(par=p,Ts=Ts,zs=zs)* ( rho_a *2.5e3 )

    err = rmse(tcw-tcw_obs)
    return(err)
}

# get_gamma <- function(plev,p_t,p1=750,p2=500)
# {
#     ph = calc_elev(plev)
#     l1 = which(plev==p1)
#     l2 = which(plev==p2)
#     gamma = (p_t[,,l2,,]-p_t[,,l1,,])/(ph[l2]-ph[l1])
#     return(gamma)
# }
get_gamma <- function(plev,p_t,p1=750,p2=500)
{
    ph = calc_elev(plev)
    l1 = which(plev==p1)
    l2 = which(plev==p2)
    gamma = (p_t[l2,,,]-p_t[l1,,,])/(ph[l2]-ph[l1])
    return(gamma)
}


# get_t2m <- function(plev,p_t,zs,p1=750,p2=500)
# {
#     ph = calc_elev(plev)
#     l1 = which(plev==p1)
#     l2 = which(plev==p2)
#     t2m = p_t[,,l1,,] + (p_t[,,l2,,]-p_t[,,l1,,])/(ph[l2]-ph[l1])*(zs-ph[l1])
#     return(t2m)
# }

# get_t2m <- function(plev,p_t,p_Z,zs,p1=750,p2=500)
# {
#     l1 = which(plev==p1)
#     l2 = which(plev==p2)
#     t2m = p_t[,,l1,,] + (p_t[,,l2,,]-p_t[,,l1,,])/(p_Z[,,l2,,]-p_Z[,,l1,,])*(zs-p_Z[,,l1,,])

#     return(t2m)
# }

get_t2m <- function(plev,p_t,p_Z,zs,p1=750,p2=500)
{
    l1 = which(plev==p1)
    l2 = which(plev==p2)
    t2m = p_t[l1,,,] + (p_t[l2,,,]-p_t[l1,,,])/(p_Z[l2,,,]-p_Z[l1,,,])*(zs-p_Z[l1,,,])

    return(t2m)
}

# get_tsl <- function(plev,p_t,p1=750,p2=500)
# {
#     ph = calc_elev(plev)
#     l1 = which(plev==p1)
#     l2 = which(plev==p2)
#     tsl = p_t[,,l1,,] + (p_t[,,l2,,]-p_t[,,l1,,])/(ph[l2]-ph[l1])*(0-ph[l1])
#     return(tsl)
# }

get_tsl <- function(plev,p_t,p1=750,p2=500)
{
    ph = calc_elev(plev)
    l1 = which(plev==p1)
    l2 = which(plev==p2)
    tsl = p_t[l1,,,] + (p_t[l2,,,]-p_t[l1,,,])/(ph[l2]-ph[l1])*(0-ph[l1])
    return(tsl)
}

calc_P <- function(tcw,w,tau=5,kappa=10,sec_day=8.64e4) 
{   
    # tcw (kg/m2)==(mm)
    # all precipitates in 5 days, mm/s
    # multiply w sec_day to get mm/d
    P = tcw/(tau*sec_day)*(sec_day) * (1 + kappa*w)
    P[P<0] = 0 

    return(P)
}

calc_pp_roe <- function(Ts,w,a=2.5e-11,b=5.9e-9,sec_day=8.64e4,alpha=0.01)
{
    
    # Get the saturation vapor pressure (kg/m s2)
    esat = calc_esat(Ts=Ts)
    
    # Loop over normal distribution of winds (sum to get precip)
    pp0 = w*NA
    ii = which(!is.na(w))
    for (i in ii) {
        wtmp = seq(w[i]-3*alpha,w[i]+3*alpha,length.out=100)
        wdist = dnorm(wtmp,mean=w[i],sd=alpha)
        wdist = wdist / sum(wdist)
        
        linw = a + b*wtmp
        linw[linw<0] = 0 

        pp0[i] = sum( linw*wdist )

    }
    
    # Scale precip by vapor pressure (kg/m s2) * (m2 s / kg) => (m/s)
    pp = esat * pp0

    # Convert units (m/s) => (mm/day)
    pp = pp *sec_day *1e3
    
    # If working on a grid, make sure to dimensionalize pp
    dim(pp) = dim(Ts)

    return(pp)
}

calc_precip <- function(qq=NULL,Ts,w=0,u=0,uv=0,zs=0,dzs=0,kx=1,kz=1,kz1=0,tau=5,sec_day=8.64e4)
{

    if (is.null(qq)) {
        qsat  = calc_qsat(Ts=Ts,zs=zs)
        qq = qsat * 1               # Saturation scaled by relative humidity
    }

    # Get total water in column (kg/kg) => (kg/m2)==(mm)
    # kg/kg * kg/m3 * m => ? kg/m2 == mm
    rho_a = calc_airdens(zs=zs)
    tcw = qq * ( rho_a *2.5e3 )
    
    # Fc   = 0.5*(1+tanh(kz*w))
    # tau1 = tau*(1-kz1*Fc)
    # pp = tcw/(tau1*sec_day) + kz*dzs

    #pp = tcw/(tau*sec_day) * (1 + kz*dzs)
    pp = tcw/(tau*sec_day)  * (1 + kz*dzs + kx*w)
    #pp = tcw/(tau*sec_day) * ( 1 + kx*w) #*kz #*dzs
    #pp  = tcw/(tau*sec_day)

    pp[pp<0] = 0 
    
    # mm/s => mm/d
    pp = pp*sec_day 

    return(pp)
}

opt.pp <- function(par,pp,subset=c(1:length(pp)),type,Ts,w,uv=1,u=0,qq=NULL,zs=NULL,dzs=0,...)
{
    if (type=="roe") {
        pp_new = calc_pp_roe(Ts=Ts,w=w,a=par[1],b=par[2],alpha=par[3])
    } else {
        # pp_new = calc_precip(qq=qq,Ts=Ts,w=w,u=u,uv=uv,zs=zs,dzs=dzs,kz=par[1],kx=0,tau=5)
        pp_new = calc_precip(qq=qq,Ts=Ts,w=w,u=u,uv=uv,zs=zs,dzs=dzs,kz=par[1],kx=14,tau=par[2])
        # pp_new = calc_precip(qq=qq,Ts=Ts,w=w,u=u,uv=uv,zs=zs,dzs=dzs,kz=par[1],kx=par[2],tau=par[3]) 
        # pp_new = calc_precip(qq=qq,Ts=Ts,w=w,u=u,uv=uv,zs=zs,dzs=dzs,kz=par[1],kx=0,tau=par[2]) 
    } 
    return(rmse(pp_new[subset]-pp[subset]))
}

calc_precip <- function(p,data,subset=c(1:dim(data)[1]),opt=FALSE)
{
    # pp = (data$ccw/p[1])*(1+p[2]*data$dzs)
    # pp = data$ccw/(p[1]+data$dzs*p[2])
    # pp = data$ccw/(p[1]+data$tt*p[2]+data$ww*p[3])
    # pp = data$ccw/(p[1]+data$tt*p[2]+data$ww*0)
    # pp = data$ccw/(p[1]*data$tt)*(1+data$dtt*p[2]+data$ww*p[3])
    # pp = data$ccw/110*(1+p[1]*data$tt^0.25) #+0.2*data$ww) #+data$dzs*p[3])
    # pp = data$ccw * (p[1]+data$tt*p[2])
    # pp = data$ccw * (p[1]+data$tt*p[2]+data$zs*p[3])
    # pp = data$ccw/(p[1] + p[2]*data$tt)*dat$q_r
    # pp = data$ccw/p[1]*1/data$airdens
    # pp = data$ccw*(1 - p[1]*data$tt)/p[2]*1/data$airdens + p[3]*data$dzs
    # pp = data$ccw/(p[1])*(1+data$dttdxy*-2)
    # pp = data$ccw*(p[1]+data$Udtt*p[2])
    pp = data$ccw*(1 - p[1]*data$tt)/p[2]*1/data$airdens

    # pp = (data$ccw/p[1])*data$cc *1/data$airdens 
    # pp = (data$ccw/p[1])*(1/data$airdens) * (1+data$cc*p[2]+data$ww*p[3])

    # pp = (data$ccw/p[1]) * (p[2]+p[3]*dat$tcc+p[4]*dat$t2m)  /dat$rho_a

    if (opt) {
        err = rmse(pp-data$pp,ii=subset)
        return(err)
    } else {
        return(pp)
    }
}

calc_precip <- function(p,data,subset=c(1:dim(data)[1]),opt=FALSE)
{

    # pp = (p[1]*data$ccw)*(1 - p[2]*data$t3m) #*1/data$airdens
    
    # BEST:    pp = p[1]*data$ccw * (1 - p[2]*data$t3m) *data$cc
    # SIMPLER: pp = p[1]*data$ccw * (1 - p[2]*data$t3m) 
    pp = p[1]*data$ccw * (1 - p[2]*(data$t3m)) * data$cc

    if (opt) {
        err = rmse(pp-data$pp,ii=subset)
        return(err)
    } else {
        return(pp)
    }
}

calc_albedo <- function(p,data,subset=c(1:dim(data)[1]),opt=FALSE)
{

    amin = p[1]
    amax = p[2]
    afac = p[3]
    tmid = p[4]
    # e = p[5]
    
    # a1 = rep(a,length(data$mask))
    # a1[which(data$mask==2)] = e 
    alb = amin+(amax-amin)*(0.5*tanh(afac*(data$t2m-tmid))+0.5) 
    # alb[which(alb<amin)] = amin
    
    if (opt) {
        err = rmse(alb-data$al_s,ii=subset)
        return(err)
    } else {
        return(alb)
    }
}

calc_albedo_mario <- function(p,data,subset=c(1:dim(data)[1]),opt=FALSE,sec_day=86400)
{

    tmax     = p[1]
    tmin     = p[2]
    alb_smax = p[3]
    alb_smin = p[4]

    # flexible factor ensures continuous polynomial
    f = 1.0/(tmax-tmin)

    tm  = data$ts*0.0
    tm[which(data$ts > tmax)] = 1.0 
    ii = which(data$ts >= tmin & data$ts < tmax)
    tm[ii] = f*(data$ts[ii] - tmin)
    

    # In contrast to the formulation in their paper, I summed up alpha_nir
    # and alpha_nir immediately (fewer parameters: alb_smax and alb_smin).
    alb = alb_smax - (alb_smax - alb_smin)*tm^3

    if (opt) {
        err = rmse(alb-data$al,ii=subset)
        return(err)
    } else {
        return(alb)
    }
}

################

## CLOUDINESS ##

################

calc_cc = function(par,dat,dims=NULL,obs=NULL)
{
    k1 = par[1]
    k2 = par[2]
    k3 = par[3]
    k4 = par[4]
    k5 = par[5]
    k6 = par[6]

    # Fw = 0.5*(1-tanh(ww/k1))
    # cc = q_r^k2 * (k3 + k4*Fw)

    # Fw = 0.5*(1-tanh(k1*dat$ww))
    # cc = dat$q_r^1.5 * (k2 + k3*Fw)

    #cc = (1-k2*dat$rho_a) * (k3 + k4*Fw)
    # cc = k1 + k2*dat$q_r^k3 #* (k3+k4*Fw)
    # cc = k2*dat$q_r*(1-k3*dat$tcw) * (k4 + k5*Fw)

    # cc = 1 - ( (1-dat$q_r) / (1-k1) )^k2

    # cc =  0.5*(tanh( k1*dat$q_r+k2 )+1)   #+(k1+k2*dat$t2m)

    #cc =  ((k1+k2*dat$ccw)/dat$rho_a)^k3*(k4+k5*dat$t2m0)  # Good
    #cc =  ((k1+k2*dat$ccw)*(k4+k5*dat$t2m0)/dat$rho_a)^k3  # A little better
    cc =  (k1+k2*dat$ccw+k3*dat$t2m)  /dat$rho_a      # Similar, but cleaner, less parameters

    # Checking for Andrey
    # cc =  (k1+k2*dat$ccw+k3*dat$tcw_sat) / dat$rho_a #+ k3*dat$zs

    if (!is.null(obs)) {
        err = rmse(cc-obs)
        return(err)
    } else {   
        if (!is.null(dims)) dim(cc) = dims
        return(cc)
    }       
}


#####################################

#### DOWNWARD LONGWAVE RADIATION ####

#####################################

calc_ebs_sky   <- function(w,alpha,beta,m) 1 - (1+w/10)*exp(-(alpha+beta*w/10)^m)

calc_dlr_trigo <- function(tt,tcw,cc=1,mask=rep(tt,1),abm1,abm2)
{
    # Stefan-Boltzmann constant: 5.670373(21)×10−8 W m^−2 K^−4
    sigma = 5.670373e-8

    # Convert variables to Trigo format
    Tsky = tt
    w    = tcw
    n    = cc 

    # Get the skin temperature (2m temperature)
    # in our case, this is an input!!
    #Tsky = T_0 + (gamma*dT*d_0+delta)
    
    # Get emissivity
    ebs_sky_cloudy = calc_ebs_sky(w=w,alpha=abm1[1],beta=abm1[2],m=abm1[3])
    ebs_sky_clear  = calc_ebs_sky(w=w,alpha=abm2[1],beta=abm2[2],m=abm2[3])

    # Get flux downward
    Fd_cloudy = sigma * ebs_sky_cloudy * Tsky^4
    Fd_clear  = sigma * ebs_sky_clear  * Tsky^4
    
    dlr = Fd_cloudy*n + Fd_clear*(1-n)
    
    dim(dlr) = dim(tt)

    return(dlr)
}

opt.calc_dlr_trigo <- function(par,tt,tcw,cc,lwd,all=TRUE)
{
    if (all) {
        lwd.new = calc_dlr_trigo(tt=tt,tcw=tcw,cc=cc,abm1=par[1:3],abm2=par[4:6])
    } else {
        lwd.new = calc_dlr_trigo(tt=tt,tcw=tcw,cc=cc*0+1,abm1=par[1:3],abm2=c(0,0,0))
    }

    err = rmse(lwd.new-lwd)

    return(err)
}


## ADDITIONAL FUNCTIONS ##

# Diffuse a field to make it smoother
diffuse <- function(var,iter=1)
{
    nx = dim(var)[1]
    ny = dim(var)[2]

    for (q in 1:iter) {

        var1 = var 
        for (i in 2:(nx-1)) {
            for (j in 2:(ny-1)) {
                var1[i,j] = mean(var[(i-1):(i+1),(j-1):(j+1)],na.rm=TRUE)
            }
        }

        var = var1
    }

    return(var1)
}


mean.mask <- function(var,mask)
{
    ave = mean(var[mask],na.rm=T)
    return(ave)
}

