library('EMD')

### Empirical Mode Decomposition
ndata <- 3000
tt2 <- seq(0, 9, length=ndata)

Xarma1<-arima.sim(list(order=c(1,0,1), ar=c(-.5), ma=c(-.4)), n=length(tt2), rand.gen = rt, df=4)
xt3 <- sin(8 * tt2) + sin(4 * tt2) 
#+ 0.1*Xarma1
try <- emd(xt3, tt2, boundary="wave")
### Ploting the IMF's
par(mfrow=c(3, 1), mar=c(2,1,2,1))
rangeimf <- range(try$imf)
plot(tt2, xt3,  type="l", main="Señal", xaxt = "n", yaxt = "n")
plot(tt2, try$imf[,1], type="l", xlab="", ylab="", ylim=rangeimf,
       main=paste(1, "er IMF", sep=""), xaxt = "n", yaxt = "n"); abline(h=0)
lines(tt2, sin(8 * tt2), col="red")
plot(tt2, try$imf[,2], type="l", xlab="", ylab="", ylim=rangeimf,
     main=paste(2, "do IMF", sep=""), xaxt = "n", yaxt = "n"); abline(h=0)
lines(tt2, sin(4 * tt2), col="red")
#plot(tt2, try$residue, xlab="", ylab="", main="Residuo", type="l", axes=FALSE); box()
