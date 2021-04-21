set.seed(1)

rm(list=ls())

this.dir = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(this.dir)

library(INLA)
library(spData)
library("R.matlab")
library(graphics)
library(ggplot2)
library(maps)
library(nlsr)

load("INLAData.R")

# Data selection

minx = 65
maxx = 90
miny = 95
maxy = 120

cond = (dfS1$coord.x>=minx)&(dfS1$coord.x<=maxx)&(dfS1$coord.y>=miny)&(dfS1$coord.y<=maxy)
dfsR = dfS1[cond,]
coor = dfsR[,c("coord.x","coord.y")]

# Mesh creation

domain <- matrix(cbind(c(minx,maxx,maxx,minx), c(miny,miny,maxy,maxy)),ncol=2)
mesh <- inla.mesh.2d(loc=coor,loc.domain = domain,max.edge=c(0.8, 0.8), cutoff=3, offset = c(1,2))

# Model creation

spde <- inla.spde2.matern(mesh=mesh, alpha=2)
mesh.index <- inla.spde.make.index(name="spatial.field",n.spde=spde$n.spde)
A1 <- inla.spde.make.A(mesh, loc=as.matrix(coor))


stack.est2 <- inla.stack(data=list(Pixel_Value=dfsR$Pixel_Value),A=list(A1,2),
                         effects=list(c(mesh.index,
                                        list(intercept=1)),list(Substance=dfsR$Substance,Lambda=dfsR$Lambda)), tag="est2")

formula.csp1 = Pixel_Value ~ Substance + Lambda + coord.x + coord.y
formula.csp2 = Pixel_Value ~ -1 +  intercept + Substance + Lambda + f(spatial.field, model=spde)

output_stat1 = inla(formula.csp1, family="gaussian", control.predictor=list(compute=TRUE), data=dfsR,verbose=TRUE)

output_stat2 = inla(formula.csp1, family="gamma", control.predictor=list(compute=TRUE), data=dfsR,verbose=TRUE)

output_stat3 <- inla(formula.csp2,data=inla.stack.data(stack.est2,spde=spde),
                     family="gaussian",control.family=list(link="identity"),
                     control.predictor=list(A=inla.stack.A(stack.est2),compute=TRUE,link=1),
                     verbose=TRUE)

output_stat4 <- inla(formula.csp2,data=inla.stack.data(stack.est2,spde=spde),
                     family="gamma",control.family=list(link="log"),
                     control.predictor=list(A=inla.stack.A(stack.est2),compute=TRUE,link=1),
                     verbose=TRUE)


# Plotting

x = readMat(con = "Crisco-1234-1.mat")
y = apply(x$ImageStack[,1:206,],MARGIN=c(1,2),FUN=sum)-500*x$Background[,1:206]
mask = y/y-1
data = array(0,dim=c(156,206,500))
for (i in 1:500){
  data[,,i] = x$ImageStack[,1:206,i]-x$Background[,1:206]+mask
}
CriscoS11170 = data
x = CriscoS11170[70,70,]
xaxis=seq(13,500)
xr = x[13:500]
mod = nlxb(xr~b*(1-exp(-k*xaxis))+a,start =list(a=250,b=500,k=0.005))
s = summary(mod)
p = c(s$param[1,1],s$param[2,1],s$param[3,1])
pupper = p+c(s$param[1,2],s$param[2,2],s$param[3,2])*1.96
plower = p-c(s$param[1,2],s$param[2,2],s$param[3,2])*1.96
yfit = p[2]*(1-exp(-p[3]*xaxis))+p[1]
yupper = pupper[2]*(1-exp(-pupper[3]*xaxis))+pupper[1]
ylower = plower[2]*(1-exp(-plower[3]*xaxis))+plower[1]
plot(x, type = 'n',xlab = "Time (Frames)",ylab = "Light Intensity (arbitrary unit)",main = "Light Intensity over Time",)
polygon(c(rev(xaxis), xaxis), c(rev(yupper), ylower), col = 'grey80', border = NA)
lines(x,type="l",lwd=2,col="black")
arrows(x0 = c(100,100,500),y0 = c(0,100,100),x1 = c(15,20,500),y1 = c(0,200,400),lwd=2,col="blue")
text(x=c(105,105),y=c(105,0),labels = c("Reflectance Value","Laser is turned on at this instant"),adj = c(0,0))
text(495,100,"Increase due to heating of the sample",adj = 1)
lines(xaxis,yfit,col="red",lwd=2)
lines(xaxis,yupper,col="red",lwd=2,lty='dashed')
lines(xaxis,ylower,col="red",lwd=2,lty='dashed')

dfra = df1[df1$Lambda==0,]
dfra = dfra[!is.na(dfra$Pixel_Value),]
par(mar = c(2, 2, 2, 2))
hist(dfra$Pixel_Value,breaks = 100,freq = F,xlim=c(0,2.5),main="Histogram of Pixel Value",xlab="Pixel Value")
lines(seq(0, 3, by=0.005), dnorm(seq(0, 3, by=0.005),mean(dfra$Pixel_Value), sd(dfra$Pixel_Value)), col="red",lwd="3")
lines(density(dfra$Pixel_Value), col="blue",lwd="3")
legend(1.75,1,c("Empirical Distribution","Gaussian Distribution"),col=c("blue","red"),lty = 1,lwd=3)

plot(mesh)
points(coor, pch=21, bg=1, col="white", cex=1.8)

plot(ecdf(dfsR$Pixel_Value),lwd=3,xlim=c(0,2.5),xlab="Pixel Value",main="Experimental CDF")
lines(ecdf(output_stat1$summary.linear.predictor[,1]),col="blue",lwd=3)
lines(ecdf(exp(output_stat2$summary.linear.predictor[,1])),col="red",lwd=3)
legend(1.5,0.6,c("Original Data","Gaussian","Gamma"),col=c("black","blue","red"),lty = 1,lwd=3)

plot(ecdf(dfsR$Pixel_Value),lwd=3,xlim=c(0,2.5),xlab="Pixel Value",main="Experimental CDF")
lines(ecdf(output_stat3$summary.linear.predictor[,1]),col="blue",lwd=3)
lines(ecdf(exp(output_stat4$summary.linear.predictor[,1])),col="red",lwd=3)
legend(1.5,0.6,c("Original Data","Gaussian","Gamma"),col=c("black","blue","red"),lty = 1,lwd=3)