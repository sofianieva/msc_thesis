library('SynchWave')
library("forecast")
library('Metrics')

# Definimos las los componentes y el dominio
t <- seq(0,10, (1/100))
t2 <- seq(-4,14, (1/100))

s11 <- function(t) (2.5*cos(2*pi*t))
s12 <- function(t) (3*cos(2*pi^2*t))
s1 <- function(t) s11(t) + s12(t)

T1 <- function(t) (8*(1/(1+(t/5)^2)+exp(-t/10)))

# Hacemos 101 iteraciones en lugar de 100 para que la mediana 
# sea uno de los valores presentes en el vector de errores y no el promedio de dos valores

# Simulación SST
errores_s11 <- c()
errores_s12 <- c()
errores_T1 <- c()
errores_r <- c()
tiempos <- c()

wavelet <- "hhhat"
# "bump", "hhhat" y "morlet" funcionan

for (i in 1:101){
  set.seed(i)
  Xarma1f <- function(t) arima.sim(list(order=c(1,0,1), ar=c(-.5), ma=c(-.4)), n=length(t), rand.gen = rt, df=4)
  sigma <- function(t) (1+.1*cos(pi*t))
  X1f <- function(t) (2*sigma(t)*Xarma1f(t))
  X1 <- X1f(t2)
  
  Y0 <- function(t) s1(t) + T1(t) + X1
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
  opt$gamma <- est_riskshrink_thresh(cwtfit$Wx, nv)
  Y0_sinT <- cwt_iw(cwtfit$Wx, opt$type, opt)
  T1_hat <- Y0t - Y0_sinT  
  r_hat = Y0t - curvefit[,1] - curvefit[,2] - T1_hat
  
  # finalizo e "cronometro"
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  
  error_s11 = min(rmse(s11(t), curvefit[,1][401:1401]), rmse(s11(t), curvefit[,2][401:1401]))
  error_s12 = min(rmse(s12(t), curvefit[,2][401:1401]), rmse(s12(t), curvefit[,1][401:1401]))
  error_T1 = rmse(T1(t),  T1_hat[401:1401])
  error_r = rmse(X1[401:1401],  r_hat[401:1401])
  errores_s11 <- c(errores_s11, error_s11)
  errores_s12 <- c(errores_s12, error_s12)
  errores_T1 <- c(errores_T1, error_T1)
  errores_r <- c(errores_r, error_r)
  tiempos <- c(tiempos,  time.taken)
}

ms11 = mean(errores_s11)
sds11 = sd(errores_s11)
ms12 = mean(errores_s12)
sds12 = sd(errores_s12)
mT1 = mean(errores_T1)
sdT1 = sd(errores_T1)
mr = mean(errores_r)
sdr = sd(errores_r)
mtime = mean(tiempos)
sdtime = sd(tiempos)

erroresY0 = errores_s11 + errores_s12 + errores_T1 + errores_r
error <- round(erroresY0, 5)
m = median(error)

# Vuelvo a generar la descomposición de la realización que obtuvo la mediana de los errores para graficarla
j <- which(error == m)[1]
set.seed(j)
Xarma1f <- function(t) arima.sim(list(order=c(1,0,1), ar=c(-.5), ma=c(-.4)), n=length(t), rand.gen = rt, df=4)
sigma <- function(t) (1+.1*cos(pi*t))
X1f <- function(t) (2*sigma(t)*Xarma1f(t))
X1 <- X1f(t2)

Y0 <- function(t) s1(t) + T1(t) + X1
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
opt$gamma <- est_riskshrink_thresh(cwtfit$Wx, nv)
Y0_sinT <- cwt_iw(cwtfit$Wx, opt$type, opt)
T1_hat <- Y0t - Y0_sinT  
r_hat = Y0t - curvefit[,1] - curvefit[,2] - T1_hat


# Simulación TBATS
Terrores_s11 <- c()
Terrores_s12 <- c()
Terrores_T1 <- c()
Terrores_r <- c()
Ttiempos <- c()

for (i in 1:101){
  set.seed(i)
  Xarma1f <- function(t) arima.sim(list(order=c(1,0,1), ar=c(-.5), ma=c(-.4)), n=length(t), rand.gen = rt, df=4)
  sigma <- function(t) (1+.1*cos(pi*t))
  X1f <- function(t) (2*sigma(t)*Xarma1f(t))
  X1 <- X1f(t)
  
  Y0 <- function(t) s1(t) + T1(t) + X1
  TY0t <- Y0(t)
  
  
  # inicializo el "cronometro"
  start.time <- Sys.time()
  
  #TBATS en Y0
  fit <- tbats(TY0t, seasonal.periods=c(100, 100/pi), use.box.cox=TRUE, use.trend=FALSE) 
  components <- tbats.components(fit)
  r_hat_TBATBS <- TY0t - components[,'level'] - components[,'season1'] - components[,'season2']
  
  # finalizo e "cronometro"
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  
  error_s11 = min(rmse(s11(t), components[,'season1']), rmse(s11(t), components[,'season2']))
  error_s12 = min(rmse(s12(t), components[,'season1']), rmse(s12(t), components[,'season2']))
  error_T1 = rmse(T1(t), components[,'level'])
  error_r = rmse(X1, r_hat_TBATBS)
  Terrores_s11 <- c(Terrores_s11, error_s11)
  Terrores_s12 <- c(Terrores_s12, error_s12)
  Terrores_T1 <- c(Terrores_T1, error_T1)
  Terrores_r <- c(Terrores_r, error_r)
  Ttiempos <- c(Ttiempos,  time.taken)
}

Tms11 = mean(Terrores_s11)
Tsds11 = sd(Terrores_s11)
Tms12 = mean(Terrores_s12)
Tsds12 = sd(Terrores_s12)
TmT1 = mean(Terrores_T1)
TsdT1 = sd(Terrores_T1)
Tmr = mean(Terrores_r)
Tsdr = sd(Terrores_r)
Tmtime = mean(Ttiempos)
Tsdtime = sd(Ttiempos)

TerroresY0 = Terrores_s11 + Terrores_s12 + Terrores_T1 + Terrores_r
Terror <- round(TerroresY0, 4)
Tm = median(Terror)

# Vuelvo a generar la descomposición de la realización que obtuvo la mediana de los errores para graficarla
k <- which(Terror == Tm)[1]
set.seed(k)
TX1 <- X1f(t)

Y0 <- function(t) s1(t) + T1(t) + TX1
TY0t <- Y0(t)

#TBATS en Y0
fit <- tbats(TY0t, seasonal.periods=c(100, 100/pi), use.box.cox=TRUE, use.trend=FALSE) 
components <- tbats.components(fit)
r_hat_TBATBS <- TY0t - components[,'level'] - components[,'season1'] - components[,'season2']


# FIGURA 5.6
par(mfrow=c(5,2), mai=c(0.3,0.3,0.1,0.15))

plot(t, TY0t, type="l", xaxt = "n")
plot(t, Y0t[401:1401], type="l", xaxt = "n")

plot(t,as.double(components[,'season2']), type = "l", ylab="s11",col="red", xaxt = "n")
lines(t, s11(t), lty=1)
plot(t, curvefit[,1][401:1401], type="l", col="red", xaxt = "n")
lines(t, s11(t), lty=1)

plot(t,as.double(components[,'season1']), type = "l", ylab="s12",col="red", xaxt = "n")
lines(t, s12(t), lty=1)
plot(t, curvefit[,2][401:1401], type="l", col="red", xaxt = "n")
lines(t, s12(t), lty=1)

plot(t,as.double(components[,'level']), type = "l", ylab="tendencia",col="red", xaxt = "n")
lines(t, T1(t), lty=1)
plot(t, T1(t), type="l", ylim=c(4,17), xaxt = "n")
lines(t, T1_hat[401:1401], col="red", lty=1)

set.seed(j)
X1 <- X1f(t2)
plot(t, TX1, type="l")
lines(t, r_hat_TBATBS, col="red", lty=1)
plot(t, X1[401:1401], type="l")
lines(t, r_hat[401:1401], col="red", lty=1)

# TABLA 2 - medias
errores <- hux(
  Metodo  = c("SST", "TBATS"),
  s11 = c(round(ms11, 3), round(Tms11, 3)),
  s12 = c(round(ms12, 3), round(Tms12, 3)),
  T1 = c(round(mT1, 3), round(TmT1, 3)),
  r = c(round(mr, 3), round(Tmr, 3)),
  Tiempo = c(round(mtime, 3), round(Tmtime, 3))
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
