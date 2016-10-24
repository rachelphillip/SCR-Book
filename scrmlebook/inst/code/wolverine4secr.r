woverine2chist=function(dat,type="count"){
  # get traps:
  if(type=="proximity") {
    d=dim(dat$wtraps)
    usage=rep("",d[1])
    for(i in 1:d[1]) usage[i]=paste(as.character(dat$wtraps[i,4:168]),collapse="")
    tdat=list(detectorID=as.character(dat$wtraps[,1]),x=dat$wtraps[,2],y=dat$wtraps[,3],usage=usage)
    write.table(tdat,file="wtrapsProximity.txt",sep=" ",quote=FALSE,row.names=FALSE,col.names=FALSE) 
    wtraps=read.traps("wtrapsProximity.txt",detector="proximity")
  } else if(type=="count"){
    usage=apply(dat$wtraps[,4:168],1,sum)
    tdat=list(detectorID=as.character(dat$wtraps[,1]),x=dat$wtraps[,2],y=dat$wtraps[,3],usage=usage)
    write.table(tdat,file="wtrapsCount.txt",sep=" ",quote=FALSE,row.names=FALSE,col.names=FALSE)   
    wtraps=read.traps("wtrapsCount.txt",detector="count",binary.usage=FALSE)
  } else stop("Invalid type. Only `proximity' and `count' allowed.")
  # get capture history
  sex=rep(NA,dim(dat$wcaps)[1])
  for(i in 1:length(sex)) sex[i]=dat$wsex[dat$wcaps[i,2]]
  cdat=as.data.frame(cbind(dat$wcaps,sex))
  if(type=="proximity") wch=make.capthist(cdat,wtraps,fmt="trapID",noccasions=165,covnames="sex")
  else {
    cdat[,3]=1 # set all occasions to 1
    wch=make.capthist(cdat,wtraps,fmt="trapID",covnames="sex")
  }
  return(wch)
}

library(scrbook)
library(secr)

# Commands below were used to create secr objects
# -----------------------------------------------
#data(wolverine)
#wch.prox=woverine2chist(wolverine,type="proximity")
#wch.count=woverine2chist(wolverine,type="count")
#save(wch.prox,file="wch.prox.rda")
#save(wch.count,file="wch.count.rda")
# -----------------------------------------------
load("wch.prox.rda")
load("wch.count.rda")
wtraps=traps(wch.prox)
wmask=make.mask(wtraps,buffer=20000,type="trapbuffer",nx=100,ny=100)
plot(wmask)
plot(wtraps,add=TRUE)
wfit.prox=secr.fit(wch.prox,mask=wmask)
wfit.count=secr.fit(wch.count,mask=wmask)

bhmodel=list(g0~b,sigma~b)
wfit.prox.b=secr.fit(wch.prox,mask=wmask,model=bhmodel)

AIC(wfit.prox,wfit.prox.b)

# Try ignoring usage

wfit.prox.nou=secr.fit(wch.prox,mask=wmask,details = list(ignoreusage = TRUE))
wfit.prox.b.nou=secr.fit(wch.prox,mask=wmask,model=bhmodel,details = list(ignoreusage = TRUE))
