library(secr)

# Assumes working directory is .../SCR-Book/
load("./analysis/objects/ERdetfun-objects.RData")

## this to create the above:
#load("./analysis/data/ERdetfun-data.RData")
## then fit stuff, then save as below:
savelist=c("kruger.capt","kruger.capt.bin","kruger.cams","kruger.cams.bin",
"kruger.mask",
"er.hn.fit","df.hn.fit","er.bin.hn.fit","df.bin.hn.fit",
"er.hr.fit","df.hr.fit","er.bin.hr.fit","df.bin.hr.fit")
save(list=savelist,file="./analysis/objects/ERdetfun-objects.RData")

# Fit to count data
# -----------------
detfn="HHN"
er.hn.model=list(lambda0~1,sigma~1)
er.hn.fit=secr.fit(kruger.capt,model=er.hn.model,mask=kruger.mask,detectfn=detfn)

detfn="HHR"
er.hr.model=list(lambda0~1,sigma~1)
er.hr.fit=secr.fit(kruger.capt,model=er.hr.model,mask=kruger.mask,detectfn=detfn)

detfn="HN"
df.hn.model=list(g0~1,sigma~1)
df.hn.fit=secr.fit(kruger.capt,model=df.hn.model,mask=kruger.mask,detectfn=detfn)
df.hn.fit

detfn="HR"
df.hr.model=list(g0~1,sigma~1)
df.hr.fit=secr.fit(kruger.capt,model=df.hr.model,mask=kruger.mask,detectfn=detfn)
df.hr.fit

AIC(er.hn.fit,df.hn.fit,er.hr.fit,df.hr.fit)

predict(er.hr.fit)[,c("estimate","lcl","ucl")]

abs(round(100*predict(er.hr.fit)[1,"estimate"]/predict(er.hn.fit)[1,"estimate"]-100))

er.hr.est=predict(er.hr.fit)
df.hr.est=predict(df.hr.fit)
df.hn.est=predict(df.hn.fit)
er.hn.est=predict(er.hn.fit)

er.hr.est["D","estimate"]/er.hn.est["D","estimate"]
100*esa(er.hr.fit)/esa(df.hn.fit)

# Fit to binary data
# ------------------
detfn="HHN"
er.bin.hn.model=list(lambda0~1,sigma~1)
er.bin.hn.fit=secr.fit(kruger.capt.bin,model=er.bin.hn.model,mask=kruger.mask,detectfn=detfn)
er.bin.hn.fit

detfn="HN"
df.bin.hn.model=list(g0~1,sigma~1)
df.bin.hn.fit=secr.fit(kruger.capt.bin,model=df.bin.hn.model,mask=kruger.mask,detectfn=detfn)
df.bin.hn.fit

detfn="HHR"
er.bin.hr.model=list(lambda0~1,sigma~1)
er.bin.hr.fit=secr.fit(kruger.capt.bin,model=er.bin.hr.model,mask=kruger.mask,detectfn=detfn)
er.bin.hr.fit

detfn="HR"
df.bin.hr.model=list(g0~1,sigma~1)
df.bin.hr.fit=secr.fit(kruger.capt.bin,model=df.bin.hr.model,mask=kruger.mask,detectfn=detfn)
df.bin.hr.fit

AIC(er.bin.hn.fit,df.bin.hn.fit,er.bin.hr.fit,df.bin.hr.fit)

100*esa(er.bin.fit)/esa(df.bin.fit)


er.est=predict(er.fit)
df.est=predict(df.fit)
er.bin.est=predict(er.bin.fit)
df.bin.est=predict(df.bin.fit)


# fiddle with distances, etc:
distances <- function (X, Y) {
  onerow <- function (xy) {
    d <- function(xy2) sqrt(sum((xy2 - xy)^2))
    apply(Y, 1, d)
  }
  t(apply(X, 1, onerow))
}
#' @title Draws histogram.
#'
#' @description
#'  Utility function to draw histograms with more options than \code{hist} allows.
#'  
#' @param height Height of histogram bars.
#' @param breaks Locations of boundaries of histogram bins (must be 1 longer than \code{height}).
#' @param lineonly If TRUE, uses \code{\link{lines}} to draw lines on current plot; else uses 
#' \code{\link{plot}} to draw lines on new plot.
#' @param outline If TRUE, draws only the outline (profile) of the histogram; else draws each 
#' complete bar.
#' @param fill If TRUE, uses polygon() to fill barsl in this case valid arguments to polygon() 
#' are passed via argument(s) "...". If fill==FALSE, valid arguments to plot() or lines() are 
#' passed via argument(s) "..."
#' @param ylim Range of y-axis.
#' @param xlab Label for x-axis.
#' @param ylab Label for y-axis.
#' @param ... See aargument \code{fill}.
histline=function(height,breaks,lineonly=FALSE,outline=FALSE,fill=FALSE,ylim=range(height),
                  xlab="x",ylab="y",...)
{
  n=length(height)
  if(length(breaks)!=(n+1)) stop("breaks must be 1 longer than height")
  if(outline) {
    y=c(0,rep(height,times=rep(2,n)),0)
    x=rep(breaks,times=rep(2,(n+1)))
  }   else {
    y=rep(0,4*n)
    x=rep(0,4*n+2)
    for(i in 1:n) {
      y[((i-1)*4+1):(i*4)]=c(0,rep(height[i],2),0)
      x[((i-1)*4+1):(i*4)]=c(rep(breaks[i],2),rep(breaks[i+1],2))
    }
    x=x[1:(4*n)]
  }
  if(lineonly) {
    if(!fill) lines(x,y,...)
    else polygon(x,y,...)
  } else {
    if(!fill) plot(x,y,type="l",ylim=ylim,xlab=xlab,ylab=ylab,...)
    else {
      plot(x,y,type="n",ylim=ylim,xlab=xlab,ylab=ylab)
      polygon(x,y,...)
    }
  }
}



plot.er.fit=function(scrfit,dmax,binwidth,lwd=2,xlab="x",ylab="y",
                     hcol="red",hlineonly=FALSE,houtline=FALSE,hfill=TRUE) {
  ch = scrfit$capthist
  n=summary(ch)$counts["n","Total"] # number of individuals
  loc=data.frame(x=rep(0,n),y=rep(0,n)) # data frame for locations of individuals
  for(i in 1:n) loc[i,] = fxi.mode(er.fit,i=i) # add estimated modes of invivs to data frame
  dists = distances(loc,kruger.cams) # distances between indiv modes and detectors
  capts = apply(kruger.capt,c(1,3),sum) # sum over occasions
  dvec = as.vector(dists) # turn into a vector of distances
  cvec = as.vector(capts) # turn into a vector of counts
  keep = dvec<=dmax # look only within given max dist
  
  breaks=seq(0,dmax,binwidth)
  nb=length(breaks)-1
  mean.n=rep(0,nb)
  for(i in 1:nb) {
    keep=which(breaks[i]<=dvec & dvec<breaks[i+1])
    mean.n[i]=mean(cvec[keep])
  }
  
  # gymnastics to avoid plotting while getting plot values in pl
  ff <- tempfile()
  png(filename=ff)
  pl=plot(scrfit,xval=seq(0,dmax,length=200),lwd=2)
  dev.off()
  unlink(ff)
  
  # now do the plot we want
  histline(mean.n,breaks,xlab=xlab,ylab=ylab,col=hcol,fill=hfill,lineonly=hlineonly,outline=houtline,
           ylim=c(0,max(c(mean.n,pl$y))))
  plot(scrfit,xval=seq(0,dmax,length=200),lwd=lwd,add=TRUE)
}


xlim=range(kruger.cams$x)
ylim=range(kruger.cams$y)
aspratio=abs(diff(xlim)/abs(diff(ylim)))
det2plot=c(1,3,5,7,9)
cols=c("gray","blue","green","red","purple")
pdf("./keepfigure/ERdetfun-locest.pdf",h=5,w=aspratio*5)
plot(kruger.cams$x,kruger.cams$y,col="white",xlim=xlim,ylim=ylim,xlab="Easting (m)",
     ylab="Northing (m)",asp=1)
for(i in 1:length(det2plot)) fxi.contour(er.hr.fit,i=det2plot[i],add=TRUE,col=cols[i])
points(kruger.cams$x,kruger.cams$y,col="black",pch=19,cex=0.75)
dev.off()

pdf("./keepfigure/ERdetfun-hrfit.pdf",h=3.5)
par(mar=c(4.5,5,4,2))
plot.er.fit(er.hr.fit,dmax=15000,binwidth=3000,lwd=3,xlab="Distance, d (m)",
            ylab=expression(hat(lambda)(d)))
plot(er.hn.fit,xval=seq(0,15000,length=200),lwd=3,add=TRUE,col="gray")
legend("top",legend=c("Gaussian (or half-normal)","Hazard rate"),col=c("gray","black"),lwd=3,bty="n",cex=0.8)
dev.off()

plot.er.fit(er.hr.fit,dmax=15000,binwidth=3000,lwd=3,xlab="Distance, d (m)",
            ylab=expression(hat(lambda)(d)))
plot(er.hn.fit,xval=seq(0,15000,length=200),lwd=3,add=TRUE,col="gray")
legend("top",legend=c("Gaussian (or half-normal)","Hazard rate"),col=c("gray","black"),lwd=3,bty="n",cex=0.8)

round(100*esa(erhr.fit)[1]/esa(er.fit)[1]-100)



