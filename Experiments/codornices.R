library('SynchWave')

# Leemos los archivos .txt disponibles en https://figshare.com/articles/dataset/High_resolution_locomotor_time_series_in_Japanese_quail_in_a_home_cage_environment_over_a_6_5_day_period_ALL_SERIES_/1424729
# Poner el path correspondiente para leerlo
data<-read.csv(file="R/Experimentos/Quail_5_group2.txt", header=FALSE)[, 1]
ft<-data

# para liberar memoria sobreescribimos data
data<-c(1)  
t<-seq(1,length(ft-1))

# se tienen 2 muestras por segundo y queremos promediar a 6 
# minutos => prmediamos de a 720 datos
largo<-720
data_prom<-c()

for (i in 1:floor(length(ft)/largo)){
  data_prom<-c(data_prom, mean(ft[((i-1)*largo+1):(i*largo)],na.rm = TRUE))
}

time<-seq(1,length(data_prom-1))*0.1
mov<-data_prom*100

# FIGURA 5.9
par(mfrow=c(1,1), mai=c(1,1,1,1))
plot(time ,mov, type = "l", xaxp = c(0, 144, 6),  main = "Serie promediada a 6 minutos", xlab = "Tiempo (h)", ylab = "Movimiento (%)")

# se necesita modificar el tiempo ya que se tiene un dato cada 6 miutos
tnew<-seq(1,length(data_prom-1))*60*6 
nv <- 32
opt <- list(type = "bump")

sstfit <- synsq_cwt_fw(tnew, data_prom, nv, opt)

# Figuras 5.10 y 5.11 (cambiar el archivo que se lee para generar los distintos gráficos)
par(mfrow=c(1,2), mai=c(1,1,1,1)) 

# Para poder ver los períodos se deben hacer algunas modificaciones, pero la 
# relación escala <-> período es lineal

per<-(sstfit$asc/3000)
image.plot(list(x=tnew/3600, y=per, z=t(abs(sstfit$Wx))), #log="y",
           ylim=c(0,30),
           xlab="Tiempo (h)", ylab="Período (h) ", 
           yaxp = c(0, 30,5),
            # AMARILLOS Y NARANJAS
           col=designer.colors(64, c("grey100", "lightyellow", "goldenrod1", "tan1", "orange", "tomato", "orangered", "orangered2", "orangered3", "red4", "firebrick4", "grey32")))
# PARA AGREGAR LINEAS
#c1 <- 24 * rep(1, length((tnew/3600)))
#lines((tnew/360), c1, lwd = 3, col = "red", lty = 3)
#c2 <- 12 * rep(1, length((tnew/3600)))
#lines((tnew/360), c2, lwd = 3, col = "red", lty = 3)
#c3 <- 8 * rep(1, length((tnew/3600)))
#lines((tnew/360), c3, lwd = 3, col = "red", lty = 3)
#c4 <- 6 * rep(1, length((tnew/3600)))
#lines((tnew/360), c4, lwd = 3, col = "red", lty = 3)

# Para la SST en periodos se usa que frecuencia = 1/ período y se realizan 
# algunas modificaciones
t2<-tnew/3600  # para tener en horas el tiempo
per<-((rev(1/sstfit$fs))/3600)  # para pasar la frecuencia a periodo
mat<-t(abs(sstfit$Tx))  # la matrix
mat <- mat[ ,ncol(mat):1 ]  # la matriz con las columnas al reves
image.plot(list(x=t2, y=per, z=mat), 
           ylim=c(0,30),
           xlab="Tiempo (h)", ylab="Período (h)", 
           yaxp = c(0, 30,5),
           # AMARILLOS Y NARANJAS
           col=designer.colors(64, c("grey100", "lightyellow", "goldenrod1", "tan1", "orange", "tomato", "orangered", "orangered2", "orangered3", "red4", "firebrick4", "grey32")))
# PARA AGREGAR LINEAS
c1 <- 24 * rep(1, length((t2)))
lines((t2), c1, lwd = 2,  lty = 2)
c2 <- 12 * rep(1, length((tnew/3600)))
lines((tnew/360), c2, lwd = 2,  lty = 2)
c3 <- 8 * rep(1, length((tnew/3600)))
lines((tnew/360), c3, lwd = 2,  lty = 2)
#c4 <- 6 * rep(1, length((tnew/3600)))
#lines((tnew/360), c4, lwd = 2, lty = 2)

# Ridge extraction automatico
lambda <- 1e+04
nw <- 16
imtfit <- curve_ext_multi(sstfit$Tx, log2(sstfit$fs), 3, lambda, nw)
# Reconstruction
curvefit <- curve_ext_recon(sstfit$Tx, sstfit$fs, imtfit$Cs, opt, nw)


# FIGURA 5.12 
par(mfrow=c(1,2), mai=c(1,0.8,1,0.3)) 
image.plot(list(x=t2, y=per, z=mat), 
           ylim=c(0,30),
           xlab="Tiempo (h)", ylab="Período (h)", 
           yaxp = c(0, 30,5),
           col=designer.colors(64, c("grey100", "lightyellow", "goldenrod1", "tan1", "orange", "tomato", "orangered", "orangered2", "orangered3", "red4", "firebrick4", "grey32")))
image.plot(list(x=t2, y=per, z=mat), 
           ylim=c(0,30),
           xlab="Tiempo (h)", ylab="Período (h)", 
           yaxp = c(0, 30,5),
           col=designer.colors(64, c("grey100", "lightyellow", "goldenrod1", "tan1", "orange", "tomato", "orangered", "orangered2", "orangered3", "red4", "firebrick4", "grey32")))
# Lineas extraidas automáticamente
lines(t2, 1/sstfit$fs[imtfit$Cs[,1]]/3600, lty=2, lwd=2)
lines(t2, 1/sstfit$fs[imtfit$Cs[,2]]/3600,  lty=2, lwd=2)
lines(t2, 1/sstfit$fs[imtfit$Cs[,3]]/3600,  lty=2, lwd=2)


# FIGURA 5.13
par(mfcol=c(4,2), mai=c(0.3,0.3,0.1,0.1))

plot(t2, curvefit[,1], type="l", col="red", xaxp = c(0, 144, 6),  xlab = "", ylab = "")

plot(t2, curvefit[,2], type="l", col="red", xaxp = c(0, 144, 6),  xlab = "", ylab = "")

plot(t2, curvefit[,3], type="l", col="red", xaxp = c(0, 144, 6),  xlab = "", ylab = "")

plot(t2, (curvefit[,1] + curvefit[,2] + curvefit[,3]),  type="l", col="red", xaxp = c(0, 144, 6), xlab = "Tiempo (h)")


plot(time ,mov, type = "l", xaxp = c(0, 144, 6), xlab = "", ylab = "", xaxp = c(0, 144, 6))
lines(t2, curvefit[,1]*100-min(curvefit[,1]*100), lty=1, col="red", xaxp = c(0, 144, 6), lwd=2)

plot(time ,mov, type = "l", xaxp = c(0, 144, 6), xlab = "", ylab = "", xaxp = c(0, 144, 6))
lines(t2, curvefit[,2]*100-min(curvefit[,2]*100), lty=1, col="red", xaxp = c(0, 144, 6), lwd=2)

plot(time ,mov, type = "l", xaxp = c(0, 144, 6), xlab = "", ylab = "", xaxp = c(0, 144, 6))
lines(t2, curvefit[,3]*100-min(curvefit[,3]*100)-2, lty=1, col="red", xaxp = c(0, 144, 6), lwd=2)

plot(time ,mov, type = "l", xaxp = c(0, 144, 6), xlab = "Tiempo (h)", ylab = "")
lines(t2, (curvefit[,1] + curvefit[,2] + curvefit[,3])*100-min((curvefit[,1] + curvefit[,2] + curvefit[,3])*100)-2, lty=1, col="red", xaxt = "n", xaxp = c(0, 144, 6), lwd=2)