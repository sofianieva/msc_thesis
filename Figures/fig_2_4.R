library('EMD')

### Empirical Mode Decomposition
ndata <- 3000
tt2 <- seq(0, 9, length=ndata)
slice <- 1:2000
t_short <- tt2[slice]


# FGURA 2.4 (no guardé la semilla exacta que genera la figura de la tesis, pero esta tiene un error muy similar)
set.seed(47)
WN<-arima.sim(list(order=c(0,0,0), ar=c(), ma=c()), n=length(tt2))
xt3 <- sin(8 * tt2) + sin(4 * tt2) + 0.1*WN 
try <- emd(xt3, tt2, boundary="wave")

### Ploting the IMF's
par(mfrow=c(try$nimf+2, 1), mar=c(1,1,1,1))
rangeimf <- range(try$imf)
plot(tt2, xt3,  type="l", main="Señal", xaxt = "n", yaxt = "n")
for(i in 1:4) {
  plot(tt2, try$imf[,i], type="l", xlab="", ylab="", ylim=rangeimf,
       main=paste(i, "° IMF", sep=""), xaxt = "n", yaxt = "n"); abline(h=0)
}
plot(tt2, try$imf[,5], type="l", xlab="", ylab="", ylim=rangeimf,
     main=paste(5, "° IMF", sep=""), xaxt = "n", yaxt = "n"); abline(h=0)
lines(tt2, sin(8 * tt2), col="red")
plot(tt2, try$imf[,6], type="l", xlab="", ylab="", ylim=rangeimf,
     main=paste(6, "° IMF", sep=""), xaxt = "n", yaxt = "n"); abline(h=0)
lines(tt2, sin(4 * tt2), col="red")
for(i in 7:try$nimf) {
  plot(tt2, try$imf[,i], type="l", xlab="", ylab="", ylim=rangeimf,
       main=paste(i, "° IMF", sep=""), xaxt = "n", yaxt = "n"); abline(h=0)
}
plot(tt2, try$residue, xlab="", ylab="", main="Residuo", type="l", xaxt = "n", yaxt = "n")
