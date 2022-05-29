library('SynchWave')
library('EMD')
library('Metrics')
library('huxtable')
library('dplyr') 

# Definimos las los componentes y el dominio
t <- seq(0,10, (1/100))
t2 <- seq(-4,14, (1/100))

A1 <- function(t) (1+0.1*cos(t))*(atan(1132/87 - 200*t/87))/2 + 2
A2 <- function(t) (3.5*(t<=7.5)+2*(t>7.5))

phi1 <- function(t) (t+0.1*sin(t))
phi2 <- function(t) (3.4*t-0.02*(abs(t)^(2.3)))

s21 <- function(t) (A1(t)*cos(2*pi*phi1(t)))
s22 <- function(t) (A2(t)*cos(2*pi*phi2(t)))
s2 <- function(t) s21(t) + s22(t)

T1 <- function(t) (8*(1/(1+(t/5)^2)+exp(-t/10)))

f <- function(t) s2(t) + T1(t)


# FIGURA 5.3
par(mfrow=c(2,2), mai=c(0.4,0.7,0.3,0.3))
plot(t, s2(t)+T1(t), type = "l", ylim=c(-5,25))
plot(t, s21(t), type = "l", ylim=c(-4, 4))
plot(t, T1(t), type = "l", ylim=c(0,20))
plot(t, s22(t), type = "l", ylim=c(-4, 4))
par(mfrow = c(1, 1))


# Definimos los parámetros de Synchrosqueezing
set.seed(7)
ft <- f(t2)
dt <- t[2]-t[1]
nv <- 32
wavelet <- "hhhat"
opt <- list(type = wavelet)

# Synchrosqueezed wavelet transform
sstfit <- synsq_cwt_fw(t2, ft, nv, opt)

# Extraccion automática de curvas
lambda <- 1e+04
nw <- 16
nc <- 2
imtfit <- curve_ext_multi(sstfit$Tx, log2(sstfit$fs), nc, lambda, nw)

# Reconstruction componentes estacionales
curvefit <- curve_ext_recon(sstfit$Tx, sstfit$fs, imtfit$Cs, opt, nw)
s21_hat <- curvefit[,2][401:1401]
s22_hat <- curvefit[,1][401:1401]

# Reconstruction tendencia
cwtfit <- cwt_fw(ft, opt$type, nv, dt, opt)
ft_sinT1 <- cwt_iw(cwtfit$Wx, opt$type, opt)
T1_hat <- ft - ft_sinT1  
T1_hat <- T1_hat[401:1401]

# FIGURA 5.4
par(mfrow=c(3,1), mai=c(0.3,1,0.3,1))
plot(t, s21_hat, type="l", col="red")
lines(t, s21(t), lty=1)
plot(t, s22_hat, type="l", col="red")
lines(t, s22(t), lty=1)
plot(t, T1(t), type="l", ylim=c(4,17))
lines(t, T1_hat, col="red", lty=1)

par(mfrow=c(1,1), mai=c(1.5,1.5,1,1))
image.plot(list(x=t2, y=sstfit$fs, z=t(abs(sstfit$Tx))),
           xlab="Tiempo", ylab="Frecuencia", main="Representación tiempo-frecuencia aplicando SST",
           col=designer.colors(64, c("grey100", "grey75", "grey50", "grey0")), ylim=c(0.5, 5), xlim=c(0, 10))

dphi1 <- function(t) (1+0.1*cos(t))
dphi2 <- function(t) (3.4-0.046*abs(t)^(1.3))

plot(t2, dphi1(t2), type = "l", ylim=c(0.5, 5), xlim=c(0, 10), xlab="Tiempo", ylab="Frecuencia")
lines(t2, dphi2(t2), lty=1)
lines(t2, sstfit$fs[imtfit$Cs[,1]], col="red", lty=1)
lines(t2, sstfit$fs[imtfit$Cs[,2]], col="red", lty=1)


# Empirical Mode Decomposition
try <- emd(f(t), t, boundary="wave")

# FIGURA 5.5
par(mfrow=c(4, 1), mai=c(0.3,1,0.3,1))
rangeimf <- range(try$imf)
plot(t, try$imf[,1], type="l", xlab="", ylab="", ylim=rangeimf,
     main=paste(1, "er IMF", sep=""), col="red")
lines(t, s22(t))
plot(t, try$imf[,2], type="l", xlab="", ylab="", ylim=rangeimf,
     main=paste(2, "do IMF", sep=""), col="red")
lines(t, s21(t))
plot(t, try$imf[,3], type="l", xlab="", ylab="", ylim=rangeimf,
     main=paste(3, "er IMF", sep=""), col="red")
plot(t, try$residue, xlab="", ylab="", main="Residuo", type="l", col="red")
lines(t, T1(t), lty=1)


# Errores
error_s21_SST = rmse(s21(t), s21_hat)
error_s22_SST = rmse(s22(t), s22_hat)
error_T1_SST = rmse(T1(t), T1_hat)
error_s21_EMD = rmse(s21(t), try$imf[,2])
error_s22_EMD = rmse(s22(t), try$imf[,1])
error_T1_EMD = rmse(T1(t), try$residue)


# TABLA 1
errores <- hux(
  Metodo  = c("SST", "EMD"),
  s21 = c(round(error_s21_SST, 3), round(error_s21_EMD, 3)),
  s22 = c(round(error_s22_SST, 3), round(error_s22_EMD, 3)),
  T1 = c(round(error_T1_SST, 3), round(error_T1_EMD, 3))
)
errores%>% 
  set_all_padding(4) %>% 
  set_outer_padding(0) %>% 
  set_bold(row = 1, col = everywhere) %>% 
  set_bottom_border(row = 1, col = everywhere) %>% 
  set_width(0.4) %>% 
  set_markdown_contents(1, 1, "Método") %>% 
  set_align(1, everywhere, "center") %>%
  set_right_border(everywhere, 1, brdr(2))