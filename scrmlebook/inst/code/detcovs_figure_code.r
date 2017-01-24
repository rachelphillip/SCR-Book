ny=4; nx=4 # dimensions of mask
# set costs (NA is "hole" - nothing there & can't go there):
costs=c(100,100,100,100,NA,100,100,NA,1,NA,100,1,1,1,1,1) 
rmesh=data.frame(x=rep(1:nx,ny),y=rep(1:ny,rep(nx,ny)),noneuc=costs) # make data frame with coords and costs

rmask=read.mask(data=rmesh,columns="noneuc")  # make mask with covariate 'noneuc' containing cost

pdf("./keepfigure/igraphplot.pdf",h=6,w=12)
par(mfrow=c(1,2))
plotcovariate(rmask,covariate="noneuc",contour=FALSE,col=c("lightgray","gray"),key=FALSE) # look at the covariate surface
text(rmesh$x,rmesh$y,labels=c(13:16,9:12,5:8,1:4),cex=1.25)

layout=cbind(c(1:4,1,3,4,2,3,1:4),c(4,4,4,4,3,3,3,2,2,1,1,1,1))
vcol=c(rep("lightgray",5),"gray","lightgray",rep("gray",6))
ig=make_igraph(rmask,"noneuc")
layout=layout[order(as.numeric(names(V(ig)))),]
vcol=vcol[order(as.numeric(names(V(ig))))]
plot(ig, edge.label=signif(E(ig)$weight, 2), edge.label.cex=0.9,layout=layout,vertex.label.cex=1,vertex.size=25,vertex.color=vcol)
dev.off()


pdf("./keepfigure/lcpath.pdf",h=6,w=6)
from=c(4,1); to=c(2,4) # corrdinates of start and end points (Note: need not be exactly on mask coords)
dist = plot_lcpath(from,to,rmask,lwd=2,linecol="black",col=c("lightgray","gray"),key=FALSE) # calculate and plot the least cost path, returning cost
dev.off()

