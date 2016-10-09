library(secr)

# Assumes working directory is .../SCR-Book/
load("./analysis/objects/ERdetfun-objects.RData")

## this to create the above:
#load("./analysis/data/ERdetfun-data.RData")
## then fit stuff, then save as below:
savedatalist=c(
  "kruger.capt","kruger.capt.bin","kruger.cams","kruger.cams.bin","kruger.mask"
)
savefitlist=c(
  "er.hn.fit","df.hn.fit","er.bin.hn.fit","df.bin.hn.fit",
  "er.hr.fit","df.hr.fit","er.bin.hr.fit","df.bin.hr.fit"
)
save(list=savedatalist,file="./scrmlebook/data/ERdetfun-data.RData")
save(list=savefitlist,file="./scrmlebook/data/ERdetfun-fits.RData")

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


round(100*esa(erhr.fit)[1]/esa(er.fit)[1]-100)


er.hr.plot = plot(er.hr.fit,xval=seq(0,20000,length=200),lwd=3,xlab="Distance, d (m)",ylab=expression(lambda(d)))
er.bin.hr.plot = plot(er.bin.hr.fit,xval=seq(0,20000,length=200),add=TRUE,lwd=3,lty=2,col="gray")


round(100*esa(er.hr.fit)[1]/esa(er.bin.hr.fit)[1]-100,3)
round(100*predict(er.hr.fit)[1,"estimate"]/predict(er.bin.hr.fit)[1,"estimate"]-100,3)
round(100*predict(er.hr.fit)[1,"SE.estimate"]/predict(er.bin.hr.fit)[1,"SE.estimate"]-100,3)
round(100*cv.predict(er.hr.fit)[1]/cv.predict(er.bin.hr.fit)[1]-100,3)

predict(er.hr.fit)[1,"SE.estimate"]/predict(er.hr.fit)[1,"estimate"]
predict(er.bin.hr.fit)[1,"SE.estimate"]/predict(er.bin.hr.fit)[1,"estimate"]

predict.cv=function(scrfit) {
  pred = predict(scrfit)
  cv = pred[,"SE.estimate"]/pred[,"estimate"]
  return(cv)
}

plot(er.bin.hr.plot$x,er.bin.hr.plot$y*er.bin.hr.plot$x,type="l")
lines(er.hr.plot$x,er.hr.plot$y*er.hr.plot$x,lty=2)


lines(er.bin.hr.plot$x,1-exp(-er.bin.hr.plot$y))
lines(er.bin.hr.plot$x,-log(1-er.bin.hr.plot$y),col="blue")
lines(er.hr.plot$x,1-exp(-er.hr.plot$y))
