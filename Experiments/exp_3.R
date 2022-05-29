library('SynchWave')
library("forecast")
library("fGarch")
library('Metrics')

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
T2 <- function(t) (2*t+10*exp(-(t-4)^2/6))

Xarma1f <- function(t) arima.sim(list(order=c(1,0,1), ar=c(-.5), ma=c(-.4)), n=length(t), rand.gen = rt, df=4)
sigma <- function(t) (1+.1*cos(pi*t))
X1f <- function(t) (2*sigma(t)*Xarma1f(t))
Xarma2f <- function(t) arima.sim(list(order=c(1,0,1), ar=c(.2), ma=c(-.51)), n=length(t), rand.gen = rt, df=4)
X2f <- function(t) sigma(t)*(4*Xarma1f(t)*(t<=5)+Xarma2f(t)*(t>5)) 
l <- list(omega = c(1), alpha = c(0.2), beta=c(0.2, 0.3))
spec = garchSpec(model = l)
Xgarchf <- function(t) garchSim(spec, n = length(t))
X3f <- function(t) 2*as.numeric(Xgarchf(t))

# Acá se cambia la señal asignando distintos valores a trend, Xf y sigma0
trend <- T1
Xf <- X3f
sigma0 <- 0.5
# Hacemos 101 iteraciones en lugar de 100 para que la mediana 
# sea uno de los valores presentes en el vector de errores y no el promedio de dos valores
it <- 101

# Simulación SST
errores_s21 <- c()
errores_s22 <- c()
errores_T <- c()
errores_r <- c()
tiempos <- c()

wavelet <- "hhhat"
# "bump", "hhhat" y "morlet" funcionan

for (i in 1:it){
  set.seed(i)
  X <- Xf(t2)
  
  Y0 <- function(t) s2(t) + trend(t) + sigma0*X
  Y0t <- Y0(t2)
  dt <- t[2]-t[1]
  nv <- 32
  opt <- list(type = wavelet)
  
  # inicializo el "cronometro"
  start.time <- Sys.time()
  
  # Synchrosqueezed wavelet transform
  sstfit <- synsq_cwt_fw(t, Y0t, nv, opt)
  # Ridge extraction
  lambda <- 1e+04
  nw <- 16
  imtfit <- curve_ext_multi(sstfit$Tx, log2(sstfit$fs), 2, lambda, nw)
  # Reconstruction
  curvefit <- curve_ext_recon(sstfit$Tx, sstfit$fs, imtfit$Cs, opt, nw)
  
  #temp <- synsq_filter_pass(sstfit$Wx, sstfit$asc, 0.1, 40);
  cwtfit <- cwt_fw(Y0t, opt$type, nv, dt, opt)
  #opt$gamma <- est_riskshrink_thresh(cwtfit$Wx, nv)
  Y0_sinT <- cwt_iw(cwtfit$Wx, opt$type, opt)
  T_hat <- Y0t - Y0_sinT  
  r_hat = Y0t - curvefit[,1] - curvefit[,2] - T_hat
  
  # finalizo e "cronometro"
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  
  error_s21 = min(rmse(s21(t), curvefit[,1][401:1401]), rmse(s21(t), curvefit[,2][401:1401]))
  error_s22 = min(rmse(s22(t), curvefit[,2][401:1401]), rmse(s22(t), curvefit[,1][401:1401]))
  error_T = rmse(trend(t),  T_hat[401:1401])
  error_r = rmse(X[401:1401],  r_hat[401:1401])
  errores_s21 <- c(errores_s21, error_s21)
  errores_s22 <- c(errores_s22, error_s22)
  errores_T <- c(errores_T, error_T)
  errores_r <- c(errores_r, error_r)
  tiempos <- c(tiempos,  time.taken)
}

ms21 = mean(errores_s21)
sds21 = sd(errores_s21)
ms22 = mean(errores_s22)
sds22 = sd(errores_s22)
mT = mean(errores_T)
sdT = sd(errores_T)
mr = mean(errores_r)
sdr = sd(errores_r)
mtime = mean(tiempos)
sdtime = sd(tiempos)

erroresY0 = errores_s21 + errores_s22 + errores_T + errores_r
error <- erroresY0
m = median(error)

# Vuelvo a generar la descomposición de la realización que obtuvo la mediana de los errores para graficarla
j <- which(error == m)[1]
set.seed(j)
SX <- Xf(t2)

Y0 <- function(t) s2(t) + trend(t) + sigma0*SX
Y0t <- Y0(t2)
dt <- t[2]-t[1]
nv <- 32
opt <- list(type = wavelet)

# Ridge extraction automatico

# Synchrosqueezed wavelet transform
sstfit <- synsq_cwt_fw(t, Y0t, nv, opt)
# Ridge extraction
lambda <- 1e+04
nw <- 16
imtfit <- curve_ext_multi(sstfit$Tx, log2(sstfit$fs), 2, lambda, nw)
# Reconstruction
curvefit <- curve_ext_recon(sstfit$Tx, sstfit$fs, imtfit$Cs, opt, nw)

#temp <- synsq_filter_pass(sstfit$Wx, sstfit$asc, 0.1, 40);
cwtfit <- cwt_fw(Y0t, opt$type, nv, dt, opt)
Y0_sinT <- cwt_iw(cwtfit$Wx, opt$type, opt)
T_hat <- Y0t - Y0_sinT  
r_hat = Y0t - curvefit[,1] - curvefit[,2] - T_hat



# Simulación TBATS
Terrores_s21 <- c()
Terrores_s22 <- c()
Terrores_T <- c()
Terrores_r <- c()
Ttiempos <- c()
no_fallo <- 0

for (i in 1:it){
  set.seed(i)
  X <- Xf(t)
  
  Y0 <- function(t) s2(t) + trend(t) + sigma0*X
  TY0t <- Y0(t)
  
  
  # inicializo el "cronometro"
  start.time <- Sys.time()
  
  #TBATS en Y0
  fit <- tbats(TY0t, seasonal.periods=c(100, 100/pi), use.box.cox=TRUE, use.trend=FALSE) 
  components <- tbats.components(fit)
  if (length(components)==4004) r_hat_TBATBS <- TY0t - components[,'level'] - components[,'season1'] - components[,'season2'] else next
  
  # finalizo e "cronometro"
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  
  error_s21 = min(rmse(s21(t), components[,'season1']), rmse(s21(t), components[,'season2']))
  error_s22 = min(rmse(s22(t), components[,'season1']), rmse(s22(t), components[,'season2']))
  error_T = rmse(trend(t), components[,'level'])
  error_r = rmse(X, r_hat_TBATBS)
  Terrores_s21 <- c(Terrores_s21, error_s21)
  Terrores_s22 <- c(Terrores_s22, error_s22)
  Terrores_T <- c(Terrores_T, error_T)
  Terrores_r <- c(Terrores_r, error_r)
  Ttiempos <- c(Ttiempos,  time.taken)
  no_fallo <- no_fallo + 1
}

Tms21 = mean(Terrores_s21)
Tsds21 = sd(Terrores_s21)
Tms22 = mean(Terrores_s22)
Tsds22 = sd(Terrores_s22)
TmT = mean(Terrores_T)
TsdT = sd(Terrores_T)
Tmr = mean(Terrores_r)
Tsdr = sd(Terrores_r)
Tmtime = mean(Ttiempos)
Tsdtime = sd(Ttiempos)

TerroresY0 = Terrores_s21 + Terrores_s22 + Terrores_T + Terrores_r
Terror <- TerroresY0
Tm = median(Terror)

# Vuelvo a generar la descomposición de la realización que obtuvo la mediana de los errores para graficarla
k <- which(Terror == Tm)[1]
set.seed(k)
TX <- Xf(t)

Y0 <- function(t) s2(t) + trend(t) + sigma0*TX
TY0t <- Y0(t)

#TBATS en Y0
fit <- tbats(TY0t, seasonal.periods=c(100, 100/pi), use.box.cox=TRUE, use.trend=FALSE) 
components <- tbats.components(fit)
r_hat_TBATBS <- TY0t - components[,'level'] - components[,'season1'] - components[,'season2']


# FIGURA 
par(mfrow=c(5,2), mai=c(0.3,0.3,0.1,0.3))

plot(t, TY0t, type="l", xaxt = "n")
plot(t, Y0t[401:1401], type="l", xaxt = "n")

plot(t,as.double(components[,'season2']), type = "l", ylab="s11",col="red", xaxt = "n")
lines(t, s21(t), lty=1)
plot(t, curvefit[,2][401:1401], type="l", col="red", xaxt = "n")
lines(t, s21(t), lty=1)

plot(t,as.double(components[,'season1']), type = "l", ylab="s12",col="red", xaxt = "n", ylim = c(-4,4))
lines(t, s22(t), lty=1)
plot(t, curvefit[,1][401:1401], type="l", col="red", xaxt = "n")
lines(t, s22(t), lty=1)

plot(t,as.double(components[,'level']), type = "l", ylab="tendencia",col="red", xaxt = "n")
lines(t, trend(t), lty=1)
plot(t, trend(t), type="l", xaxt = "n")
lines(t, T_hat[401:1401], col="red", lty=1)

plot(t, TX, type="l")
lines(t, r_hat_TBATBS, col="red", lty=1)
plot(t, SX[401:1401], type="l")
lines(t, r_hat[401:1401], col="red", lty=1)
