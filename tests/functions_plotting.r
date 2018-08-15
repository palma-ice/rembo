
plot.blank <- function()
{
    plot(c(0,1),c(0,1),type="n",axes=FALSE,ann=FALSE)

    return()
}


my.image <- function(var,dat,mask=var*0==0,alpha=100,ij=rep(TRUE,length(var)),
                     zlim=NULL,col=NULL,breaks=NULL,mask.cont=TRUE,
                     col.axis="grey30")
{
  
  # Get static fields from dat
  x <- dat$xc; y <- dat$yc

  # Limit variable to points of interest
  var[!ij] = NA
  
  # Get breaks and col if not provided
  if (is.null(zlim[1])) zlim = range(var,na.rm=TRUE)
  if (is.null(breaks[1])) {
    breaks = pretty(zlim,20)
    col    = colorRampPalette(jet.colors)(length(breaks)-1)
  }
  
  # Limit variable to contour breaks
  zlim = range(breaks)
  var[var>zlim[2]] = zlim[2]
  var[var<zlim[1]] = zlim[1]

  mask2D = array(0,dim=dim(mask))
  mask2D[which(mask)] = 1

  image(x=x,y=y,z=var,axes=FALSE,col=col,breaks=breaks,xlab="",ylab="")
  image(x=x,y=y,z=mask2D,add=TRUE,col=c(alpha("white",100-alpha),NA))

  ## Add land/ice mask
  if (mask.cont) contour(x=x,y=y,z=dat$mask,add=TRUE,drawlabels=FALSE,nlevels=4,lwd=1,col=alpha(1,50))
  contour(x=x,y=y,z=dat$z_srf,add=TRUE,drawlabels=FALSE,nlevels=10,lwd=0.5,lty=1,col=alpha(1,50))
  
  ## Add lat/lon
  contour(x=x,y=y,z=dat$lat,add=TRUE,drawlabels=FALSE,nlevels=4,lwd=0.5,lty=2,col=col.axis)
  contour(x=x,y=y,z=dat$lon,add=TRUE,drawlabels=FALSE,nlevels=4,lwd=0.5,lty=2,col=col.axis)
  
  box(col=col.axis)
}

## Plot seasonal cycle
plot_year <- function(mod,obs,top=obs,mask=top$z_srf*0==0,cex=1,alpha=30,
                      vnm="tt",onm="t2m",mnm="REMBOv2",long_name="Temp. (K)",ylim=NULL,fldr="plots",type="pdf")
{

    mask3D = array(mask,dim=dim(mod[[vnm]]))
    var_mod = mod[[vnm]]
    var_mod[!mask3D] = NA 
    var_obs = obs[[onm]]
    var_obs[!mask3D] = NA 

    x = seq(15,360,by=30)
    y1 = apply(var_mod,3,mean,na.rm=TRUE)
    y0 = apply(var_obs,3,mean,na.rm=TRUE)

    xlim = c(0,360)
    if (is.null(ylim)) ylim = range(y0,y1,na.rm=TRUE)

    lwd = c(6,2)
    col = c("grey70",1)
    lty = c(1,1)

    if (!is.null(type)) myfigure(fldr,paste("seasonal_",vnm,sep=""),date=TRUE,type=type,asp=2,pointsize=12)
    par(plt=c(0.1,0.97,0.1,0.95))
    plot(xlim,ylim,type='n',axes=F,ann=F)
    grid()
    axis(1,at=x,labels=month.abb)
    axis(2)
    mtext(side=2,line=2.2,las=0,long_name)
    lines(x,y0,lwd=lwd[1],col=col[1],lty=lty[1])
    lines(x,y1,lwd=lwd[2],col=col[2],lty=lty[2])
    box()
    legend("topleft",bty="n",inset=0.01,c("Obs.",mnm),lwd=lwd,col=col,lty=lty)
    if (!is.null(type)) graphics.off()

}

## Plot monthly fields
plot_months <- function(datc,vnm="tt",dataset="erac",long_name=vnm,
                        top=datc,mask=datc$z_srf*0==0,alpha=90,ncol=15,
                        zlim=range(datc[[vnm]],na.rm=TRUE),type="pdf")
{
    breaks = pretty(zlim,ncol)

    col    = colorRampPalette(jet.colors)(length(breaks)-1)
    if (vnm %in% c("uv_s","uv","wind","p750_uv")) col = alpha(col,70)

    if (!is.null(type)) myfigure("plots",paste("monthly2D_",dataset,"_",vnm,sep=""),asp=0.8,pointsize=12,type=type)
    layout(matrix(c(1:12,13,13,13,13),nrow=4,byrow=TRUE),heights=c(1,1,1,0.3))
    par(mar=c(0.5,0.5,2,0.5))

    for (k in 1:12) {
        my.image(datc[[vnm]][,,k],dat=top,mask=mask,alpha=alpha,breaks=breaks,col=col)
        title(month.abb[k])

        if (vnm %in% c("uv_s"))    quiver(datc$x2D,datc$y2D,datc$u_s[,,k],datc$v_s[,,k],scale=80,length=25,thin=2,col="grey20")
        if (vnm %in% c("uv"))      quiver(datc$x2D,datc$y2D,datc$u[,,k],datc$v[,,k],scale=80,length=25,thin=2,col="grey20")
        if (vnm %in% c("wind"))    quiver(datc$x2D,datc$y2D,datc$ua[,,k],datc$va[,,k],scale=80,length=25,thin=2,col="grey20")
        if (vnm %in% c("p750_uv")) quiver(datc$x2D,datc$y2D,datc$p750_u[,,k],datc$p750_v[,,k],scale=80,length=25,thin=2,col="grey20")

    }
    
    plot.blank()
    par(new=TRUE,plt=c(0.2,0.8,0.3,0.6))
    mylegend(breaks=breaks,col=col,vertical=FALSE,units="")
    title(long_name)
    if (!is.null(type)) graphics.off()

}


plot_comparison <- function(mod,obs,top=obs,mask=top$z_srf*0==0,months=c(1,7),plevel=750,cex=1,alpha=30,col=jet.colors,thin = 5,
                            vnm="tt",onm="t2m",long_name="Temp. (K)",zlim0=NULL,zlim=NULL,percent=FALSE,type="pdf",fldr="plots")
{

    mask[is.na(mask)] = FALSE 
    dim(mask) = dim(top$z_srf)
    mask12 = array(mask,dim=c(dim(mask),12))

    # Modeled and observed variable (range of observed variable)
    if (is.null(zlim0)) zlim0 = range(obs[[onm]][mask12],na.rm=TRUE)
    breaks1 = pretty(zlim0,15)
    if (vnm %in% c("pp","P","Pr")) breaks1 = pretty(range(obs[[onm]]*0.8,na.rm=TRUE),20) 
    col1    = colorRampPalette(col)(length(breaks1)-1)
    # col1    = tim.colors(length(breaks1)-1)

    pp = round(abs(breaks1)/max(abs(breaks1))*100)
    pp = pp + 50 
    pp[pp>100] = 100 
    for (i in 1:length(col1)) col1[i] = alpha(col1[i],pp[i])

    # if (vnm %in% c("uv_s","uv","wind")) col1 = alpha(col1,70)

    # Difference with model
    tmp = (mod[[vnm]][mask12]-obs[[onm]][mask12])
    if (percent) {
        tmp = tmp / obs[[onm]][mask12]*100
        tmp[!is.finite(tmp)] = NA 
    }
    if (is.null(zlim)) {
        zlim_pos = abs(range(obs[[onm]][mask12],na.rm=TRUE))*0.3
        zlim = c(-zlim_pos,zlim_pos)
    }

    breaks2 = pretty(zlim,15)
    col2    = colorRampPalette(jet.colors)(length(breaks2)-1)
    col2    = colorRampPalette(col=c("darkcyan","cyan","white","magenta","darkmagenta"))(length(breaks2)-1)
    i = floor(length(col2)/2)
    col2[c(i,(i+1))] = alpha(col2[c(i,(i+1))],50)

    # Colors based on 2D info 
    zs = top$z_srf
    zs[zs<0]=0
    col3 = get.col(zs,jet.colors,n=10,extend=0)
    col3$col[mask]  = alpha(col3$col[mask],80)
    col3$col[!mask] = alpha(col3$col[!mask],alpha)

    pch = rep(20,length(col3$col))
    pch[which(mask==FALSE)] = 21
    pt.cex = rep(0.8,length(col3$col))
    pt.cex[which(mask==FALSE)] = 0.6

    xlim = range(mod[[vnm]],obs[[onm]],na.rm=TRUE)
    if (!is.null(zlim0)) xlim = zlim0 
    ylim = xlim 

    nm = length(months) 

    if (nm == 1) { ## SUMMARY PLOT, ANNUAL
        myfigure(fldr,paste("summary_",plevel,"Mb_",vnm,sep=""),asp=1.2,pointsize=10,type=type)
        layout(matrix(c(1,1,1,1,2:5,6,6,7,8,9:12),ncol=4,byrow=TRUE),widths=c(1,1,1,1.8),heights=c(0.3,1,0.3))
    } else if (nm == 2) { ## SUMMARY PLOT
        myfigure(fldr,paste("summary_",plevel,"Mb_",vnm,sep=""),asp=1.2,pointsize=10,type=type)
        layout(matrix(c(1,1,1,1,2:9,10,10,11,12),ncol=4,byrow=TRUE),widths=c(1,1,1,1.8),heights=c(0.2,1,1,0.3))
    } else if (nm == 3) { ## SUMMARY PLOT
        myfigure(fldr,paste("summary_",plevel,"Mb_",vnm,sep=""),asp=0.9,pointsize=10,type=type)
        layout(matrix(c(1,1,1,1,2:13,14,14,15,16),ncol=4,byrow=TRUE),widths=c(1,1,1,1.8),heights=c(0.5,rep(1,2,3),0.3))
    } else if (nm == 12) {       ## ALL MONTHS PLOT
        myfigure(fldr,paste("monthly_",plevel,"Mb_",vnm,sep=""),asp=0.3,pointsize=10,type=type)
        layout(matrix(c(1,1,1,1,2:49,50,50,51,52),ncol=4,byrow=TRUE),widths=c(1,1,1,1.8),heights=c(0.4,rep(1,12),0.3))
    } else {
        cat("plot_comparison:: Number of months not handled: ",nm,"\n")
        return(1)
    }

    par(mar=c(0.5,0.5,0.5,0.5),cex=cex)
    plot.blank()
    
    q = 0 
    for (k in months) {
        q = q+1 

        if (k %in% c(1:12)) {
            var0 = obs[[onm]][,,k]
            var1 = mod[[vnm]][,,k]

            if (vnm %in% c("uv_s")) { uu_mod = mod$u_s[,,k]; vv_mod = mod$v_s[,,k] }
            if (vnm %in% c("uv"))   { uu_mod = mod$u[,,k];   vv_mod = mod$v[,,k]   }
            if (vnm %in% c("wind"))  { uu_mod = mod$ua[,,k];  vv_mod = mod$va[,,k]  }
            
            if (vnm %in% c("uv_s","uv")) { uu_obs = obs$u[,,k];       vv_obs = obs$v[,,k] }
            if (vnm %in% c("wind"))       { uu_obs = obs$p750_u[,,k];  vv_obs = obs$p750_v[,,k]  }
            
        } else {
            var0 = apply(obs[[onm]],c(1,2),mean)
            var1 = apply(mod[[vnm]],c(1,2),mean)

            if (vnm %in% c("uv_s")) { uu_mod = apply(mod$u_s,c(1,2),mean); vv_mod = apply(mod$v_s,c(1,2),mean) }
            if (vnm %in% c("uv"))   { uu_mod = apply(mod$u,c(1,2),mean);   vv_mod = apply(mod$v,c(1,2),mean)  }
            if (vnm %in% c("wind"))  { uu_mod = apply(mod$ua,c(1,2),mean);  vv_mod = apply(mod$va,c(1,2),mean)  }
            
            if (vnm %in% c("uv_s","uv")) { uu_obs = apply(obs$u,c(1,2),mean);       vv_obs = apply(obs$v,c(1,2),mean) }
            if (vnm %in% c("wind"))       { uu_obs = apply(obs$p750_u,c(1,2),mean);  vv_obs = apply(obs$p750_v,c(1,2),mean)  }
            
        }

        # var0[!mask] = NA 
        # var1[!mask] = NA 

        ## Panel 1: Model results (2D)
        par(mar=c(0.5,0.5,0.5,0.5),cex=cex)
        my.image(var1,dat=top,mask=mask,alpha=alpha,breaks=breaks1,col=col1)

        if (q==1) {
            par(xpd=NA)
            title("Model",line=0.2,cex.main=0.8,col=col.axis)
            par(xpd=FALSE)
        }
        # Add arrows to plot if it is wind
        if (vnm %in% c("uv_s","uv","wind")) quiver(mod$x2D,mod$y2D,uu_mod, vv_mod, scale=90,length=25,thin=thin,col="grey20")

        ## Panel 2: Observations (2D)
        par(mar=c(0.5,0.5,0.5,0.5),cex=cex)
        my.image(var0,dat=top,mask=mask,alpha=alpha,breaks=breaks1,col=col1)

        if (q==1) {
            par(xpd=NA)
            title("Observations",line=0.2,cex.main=0.8,col=col.axis)
            par(xpd=FALSE)
        }
        
        # Add arrows to plot if it is wind
        if (vnm %in% c("uv_s","uv","wind")) quiver(mod$x2D,mod$y2D,uu_obs, vv_obs, scale=90,length=25,thin=thin,col="grey20")

        ## Panel 3: Comparison with observations (2D)
        tmp = (var1-var0)
        if (percent) {
            tmp = tmp / var0*100
            tmp[!is.finite(tmp)] = NA 
        }
        par(mar=c(0.5,0.5,0.5,0.5),cex=cex)
        my.image(tmp,dat=top,mask=mask,alpha=alpha,breaks=breaks2,col=col2)
        if (q==1) {
            par(xpd=NA)
            title("Error (mod-obs)",line=0.2,cex.main=0.8,col=col.axis)
            par(xpd=FALSE)
        }
        ## Panel 4: Comparison with observations (scatter plot)
        par(mar=c(2,2,2,1),cex=cex*0.8)
        plot(xlim,ylim,type='n',ann=F)
        grid()
        points(var0,var1,pch=pch,cex=pt.cex,col=col3$col)
        abline(a=0,b=1,lwd=2)
        text(xlim[1],ylim[2]-diff(ylim)*0.05,pos=4,"Mod.",cex=1.2)
        text(xlim[2],ylim[1]+diff(ylim)*0.05,pos=2,"Obs.",cex=1.2)
        if (k %in% c(1:12)) title(month.abb[k])
        if (k == 13)        title("Ann")

        x = as.vector(var0)
        y = as.vector(var1)
        fit = lm(y~x,subset=which(mask))
        fitsum = summary(fit) 
        abline(a=fitsum$coefficients[1],b=fitsum$coefficients[2],col=2,lwd=2)
        txt = substitute(R^2==r2,list(r2=round(fitsum$r.squared,3)))
        text(mean(xlim),ylim[2]-diff(ylim)*0.02,txt,col=2,cex=1.1)
        box()
    }
    
    # Panel 1 legend
    plot(c(0,1),c(0,1),type='n',ann=F,axes=F)
    par(new=TRUE,plt=c(0.2,0.7,0.5,0.8))
    mylegend(breaks=breaks1,col=col1,vertical=FALSE,units="",cex=0.9)
    
    # Panel 2 legend 
    plot(c(0,1),c(0,1),type='n',ann=F,axes=F)
    par(new=TRUE,plt=c(0.05,0.95,0.5,0.8))
    mylegend(breaks=breaks2,col=col2,vertical=FALSE,units="",cex=0.9)
    abline(v=0)

    # Panel 3 legend
    plot(c(0,1),c(0,1),type='n',ann=F,axes=F)
    par(new=TRUE,plt=c(0.3,0.7,0.65,0.95)) 
    mylegend(breaks=col3$breaks,col=col3$palette,vertical=FALSE,units="",cex=0.9)
    par(new=TRUE,plt=c(0,1,0.12,0.95))
    plot(c(0,1),c(0,1),type='n',ann=F,axes=F)
    text(0.5,0.12,"Elevation (m)",cex=0.9)

    par(new=TRUE,fig=c(0.05,0.95,0.05,1),plt=c(0,1,0,0.92))
    plot(c(0,1),c(0,1),type='n',ann=F,axes=F)
    title(long_name,cex=2,line=1.5)

    graphics.off()

}

gen_comparison <- function(mod,obs,mask=mod[[vnm]]*0==0,vnm="t2m",onm="t3m",long_name="Temp. (K)",type="pdf")
{

    mask12 = array(mask,dim=dim(mod[[vnm]]))
    var0 = obs[[onm]]
    var0[!mask12] = NA 
    var1 = mod[[vnm]]
    var1[!mask12] = NA 

    err = var1-var0
    return(list(mod=var1,obs=var0,err=err,sd=sd(err,na.rm=TRUE)))

}

plot_resid <- function(x,y,xlim=range(x,y,na.rm=T),ylim=xlim,title,
                       xlab="Observed",ylab="Predicted",
                       col="grey50",bg=alpha(col,50),cex=0.5,lwd=0.5,units="W m-2") {
    plot(xlim,ylim,type="n",ann=F,axes=F)
    grid()
    mtext(side=1,line=1.6,xlab)
    mtext(side=2,line=2.1,las=0,ylab)
    axis(1)
    axis(2)
    axis(3,labels=FALSE)
    axis(4,labels=FALSE)
    points(x,y,pch=21,cex=cex,lwd=lwd,col=col,bg=bg)
    abline(a=0,b=1,col=1,lwd=2)
    
    fit = lm(y~x)
    abline(coef=fit$coefficients,lwd=1,col=1,lty=2)

    err = rmse(x-y)
    text(ylim[1],ylim[2]*0.95,pos=4,paste("RMSE =",round(err,3),units))
    title(line=0.7,title)
    box()
}

plot_basin_comparison <- function(mod,obs,top=obs,mask=top$mask,mask_basin=top$basin,basin=1,month=c(1,7),cex=1,alpha=30,
                                  vnm="tt",onm="t2m",long_name="Temp. (K)",xlim=NULL,type="pdf")
{

    # Basin is the whole domain
    if (basin == 0) {
        mask_basin[mask_basin==0] = -1 
        mask_basin[mask_basin>0]  = 0 
    }

    nm = length(months) 

    if (nm == 1) { ## SUMMARY PLOT, ANNUAL
        myfigure("plots",paste0("basin_",vnm,"_",basin),asp=1.1,pointsize=12,type=type)
    } else if (nm == 2) { ## SUMMARY PLOT
        myfigure("plots",paste0("basin_",vnm,"_",basin),asp=2.2,pointsize=14,type=type)
        par(mfrow=c(1,2))
    } else if (nm == 3) { ## SUMMARY PLOT
        myfigure("plots",paste0("basin_",vnm,"_",basin),asp=3.3,pointsize=16,type=type)
        par(mfrow=c(1,3))
    } else {
        cat("plot_basin_comparison:: Number of months not handled: ",nm,"\n")
        return(1)
    }


    for (m in months) {

        if (m %in% c(1:12)) {
            tmp0 = obs[[onm]][,,m]
            tmp1 = mod[[vnm]][,,m]
        } else {
            tmp0 = apply(obs[[onm]],c(1,2),mean,na.rm=TRUE)
            tmp1 = apply(mod[[vnm]],c(1,2),mean,na.rm=TRUE) 
        }
        plot_basin(tmp1,tmp0,subset=which(mask==2 & mask_basin==basin),mask_basin,long_name,month=m)

    }

    graphics.off()
}

plot_basin <- function(mod,obs,subset=c(1:length(mod)),mask_basin=NULL,basin=1,
                       long_name,xlim=NULL,histogram=FALSE)
{
    ii = subset 
    if (is.null(xlim)) xlim = range(mod[ii],obs[ii],na.rm=TRUE)

    if (histogram) {
        breaks = pretty(xlim,30)
        breaks = seq(xlim[1],xlim[2],length.out=30)
        breaks = NULL 
        h0 = hist(obs[ii],breaks=30,plot=FALSE)
        a = list(x=h0$mids,y=h0$density)
        h1 = hist(mod[ii],breaks=30,plot=FALSE)
        b = list(x=h1$mids,y=h1$density)
    } else {
        a = density(obs[ii],n=256)
        b = density(mod[ii],n=256)
    }
    
    ylim = range(a$y,b$y,na.rm=TRUE) 

    lwd = c(6,2)
    col = c("grey70",1)
    lty = c(1,1)

    par(plt=c(0.1,0.97,0.1,0.95))
    plot(xlim,ylim,type="n",ann=FALSE)
    grid()
    mtext(side=1,line=1.6,long_name)
    mtext(side=2,line=2,las=0,"Density")

    if (histogram) {
        my.hist(h0,freq=FALSE,border=col[1])
        my.hist(h1,freq=FALSE,border=col[2])
    } else {
        lines(a$x,a$y,lwd=lwd[1],col=col[1],lty=lty[1])
        lines(b$x,b$y,lwd=lwd[2],col=col[2],lty=lty[2])
    }
        
    box()
    legend("topleft",bty="n",inset=0.01,c("Obs.","REMBOv2"),lwd=lwd,col=col,lty=lty)

}

plot_basin2D <- function(mod,obs,mask_basin,long_name)
{


}

plot_density2D <- function(x,mod,obs,subset,long_name,xlim=NULL,ylim=NULL)
{   # x, mod, obs: vectors of the same length
    # z is the frequency of mod/obs

    ii = subset 

    # Make bins (eg, temp vs massbal)
    if (is.null(xlim)) xlim = range(x[ii],na.rm=TRUE)
    if (is.null(ylim))ylim = range(mod[ii],obs[ii],na.rm=TRUE)
    
    lwd = c(6,2)
    col = c("grey70",1)
    lty = c(1,1)

    par(plt=c(0.1,0.97,0.1,0.95))
    plot(xlim,ylim,type="n",ann=FALSE)
    grid()
    mtext(side=1,line=1.6,long_name[1])
    mtext(side=2,line=2,las=0,long_name[2])

    points(x[ii],obs[ii],pch=20,col=alpha(col[1],50))
    points(x[ii],mod[ii],pch=20,col=alpha(col[2],50))


}


myfigure <- function(fldr=".",file="Rplot",date=TRUE,type="pdf",engine="cairo",
                     width=NULL,height=NULL,units="mm",asp=1,pointsize=12,res=300,
                     cex=1,cex.lab=1,cex.axis=1,bg="white",onefile=FALSE)
{
  
    # Make filename
    file = paste(file,".",type,sep="")
    if (date == TRUE) file = paste(today,"_",file,sep="")
    file = file.path(fldr,file)

    host = system("hostname",intern=TRUE)
    os   = system("uname",intern=TRUE)

    # If running on a mac, make sure engine is quartz!
    if ( os == "Darwin" ) engine = "quartz"

    # Determine width/heights in inches
    if ( is.null(width) & is.null(height) ) {  # Use default height, determine win via asp
        width  = 189  # Default width for pointsize 12
        height = width/asp
    } else if ( is.null(height) ) {  # only width specified, determine height
        height = width/asp
    } else if ( is.null(width) ) {  # only height specified, determine width
        width = asp*height
    } else {                    # height and width specified, determine asp
        asp = width/height
    }

    # Convert quantities if input was not inches
    cat(type,":",file,"\n")
    cat("width=",width,", height=",height," (",units,".) \n",sep="")
    conv = 1.0
    if ( units == "mm" ) conv = 0.0393700787
    if ( units == "cm" ) conv = 0.393700787
    hin = height*conv
    win = width*conv 
    #cat("width =",win,", height =",hin," (in.)\n")

    if (FALSE & os %in% c("Darwin") & type %in% c("png","jpg","tiff","pdf","ps")) {
        
        cat("Quartz plotting","\n")
        quartz(file=file,type=type,width=win,height=hin,pointsize=pointsize,dpi=res)

    } else if ( type == "png" ) {

        cat("engine = ",engine,"\n")
        png(file,width=win,height=hin,units="in",pointsize=pointsize,res=res,type=engine)

    } else if ( type == "jpg" ) {

        jpeg(file,width=win,height=hin,units="in",pointsize=pointsize,res=res,type=engine)

    } else if ( type == "tiff" ) {

        tiff(file,width=win,height=hin,units="in",pointsize=pointsize,res=res,type=engine)

    } else if ( type == "pdf" ) {

        if (engine %in% c("cairo","cairo1")) {
            cat("**cairo_pdf","\n")
            cairo_pdf(file,width=win,height=hin,pointsize=pointsize,onefile=onefile)
        } else {
            pdf(file,width=win,height=hin,pointsize=pointsize,onefile=onefile)
        }

    } else if ( type == "svg" ) {

        svg(file,width=win,height=hin,pointsize=pointsize)

    } else if ( type == "fig" ) {

        xfig(file,width=win,height=hin,pointsize=pointsize)

    } else {

        cairo_ps(file,width=win,height=hin,pointsize=pointsize)

    }

    par(bg=bg,cex=cex,cex.axis=cex.axis,cex.lab=cex.lab,tcl=0.2,mgp=c(2.5,0.3,0),las=1)

    return(win)
}


mylegend <- function(breaks,col,labs=paste(breaks),units="mm",x=c(0,1),y=c(0,1),
                     xlab="",ylab="",xlim=NULL,ylim=NULL,zlim=range(breaks),at=NULL,
                     cex=1,cex.lab=1,new=TRUE,vertical=TRUE,line=1.8,
                     asp=1,mgp=c(3,0.5,0),col.axis="grey10",...)
{
    n      = length(breaks)    
    ynorm  = (breaks - min(breaks))
    ynorm  = ynorm / max(ynorm)
    y00    = ynorm[1:(n-1)]
    y11    = ynorm[2:n]
    x00    = rep(0,n)
    x11    = rep(1,n)

    if ( vertical ) {
      x0   = x00
      x1   = x11 
      y0   = y00
      y1   = y11 
      xlim = c(0,1)
      ylim = zlim 
      ax   = 4
    } else {
      x0   = y00
      x1   = y11
      y0   = x00
      y1   = x11
      xlim = zlim
      ylim = c(0,1)
      ax   = 1
    }

    xlim0 = range(x0,x1)
    ylim0 = range(y0,y1)

    par(new=new,xpd=NA,xaxs="i",yaxs="i",...)
    plot( xlim0,ylim0, type="n",axes=F,ann=F,cex=cex)
    rect(x0,y0,x1,y1,col=col,border=col,lwd=1)

    par(new=TRUE,xpd=NA,xaxs="i",yaxs="i",...)
    plot(xlim,ylim,type="n",axes=F,ann=F,cex=cex)
    axis(ax,at=at,mgp=mgp,tcl=-0.1,col=col.axis,col.axis=col.axis,cex.axis=cex)
    box(col="grey10")

    mtext(side=1,line=line,xlab,cex=cex.lab)
    mtext(side=2,line=line,ylab,cex=cex.lab)

    par(xpd=FALSE)
}
