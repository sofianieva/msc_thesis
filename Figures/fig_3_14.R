library('SynchWave')

t <- seq(0,10*pi, (1/30))

s <- function(t) sin(t)
r <- function(t) 2*(arima.sim(list(order=c(1,0,1), ar=c(-.5), ma=c(-.4)), n=length(t), rand.gen = rt, df=4))

set.seed(7)
f <- function(t) s(t) + r(t)*(5*pi<t)*(t<7*pi)
ft <- f(t)
dt <- t[2]-t[1]
nv <- 32
opt <- list(type = "mhat")

# transformada wavelet continua 
cwtfit <- cwt_fw(ft, opt$type, nv, dt, opt)
par(mfrow=c(2,3), oma = c(2, 2, 2, 2), 
    mar = c(2, 2, 2, 2),
    mgp = c(2, 1, 0),   
    xpd = NA)
    
    #mai=c(0.5,0.5,0.5,0.5))

R1 = Re(cwtfit$Wx)
image.plot(list(x=t, y=cwtfit$asc, z=t(R1)), log="y",
           xlab="", ylab="", main="a*=0",
           col=designer.colors(64, c("grey90", "grey65", "grey35", "grey10")),
           ylim=c(10, 0.04), xlim=c(3, 27))

thresh2 = 0.15
indices2 = cwtfit$asc[which(cwtfit$asc < thresh2)]
R2 = Re(cwtfit$Wx)
R2[1:length(indices2), 1:943] <- 0.0
image.plot(list(x=t, y=cwtfit$asc, z=t(R2)), log="y",
           xlab="", ylab="", main="a*=0.15",
           col=designer.colors(64, c("grey90", "grey65", "grey35", "grey10")),
           ylim=c(10, 0.04), xlim=c(3, 27))
lines(t, rep(thresh2, length(t)), type='l', col='red')

thresh3 = 0.4
indices3 = cwtfit$asc[which(cwtfit$asc < thresh3)]
R3 = Re(cwtfit$Wx)
R3[1:length(indices3), 1:943] <- 0.0
image.plot(list(x=t, y=cwtfit$asc, z=t(R3)), log="y",
           xlab="", ylab="", main="a*=0.4",
           col=designer.colors(64, c("grey90", "grey65", "grey35", "grey10")),
           ylim=c(10, 0.04), xlim=c(3, 27))
lines(t, rep(thresh3, length(t)), type='l', col='red')

fcwt1 <- cwt_iw(R1, opt$type, opt);
plot(t, fcwt1, type="l", xlab='', ylab='')

fcwt2 <- cwt_iw(R2, opt$type, opt);
plot(t, fcwt2, type="l", xlab='', ylab='')

fcwt3 <- cwt_iw(R3, opt$type, opt);
plot(t, fcwt3, type="l", xlab='', ylab='')
