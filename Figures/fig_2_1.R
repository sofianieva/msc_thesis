library('EMD')

#Primer ejemplo
# Definamos las siguientes funciones para modelar la estacionalidad

# Primero hay que definir el dominio
t<-seq(0,10, (1/100))
n<-length(t)

s11t<-2.5*cos(2*pi*t)
s12t<-3*cos(2*pi^2*t)
s1t<-s11t+s12t

a1t<-2+0.5*(1+0.1*cos(t))*atan(t-13)

# para a2t primero armo un vector vacio.
a2t<-c()
#luego le voy agregando a este vector 3.5 si 0<=t<=7.5 o 2 si 7.5<t<=10
for (i in t){
  if(i<=7.5){a2t<-c(a2t, 3.5)}
  else{a2t<-c(a2t, 2)}
}

phi1t<-t+0.1*sin(t)
phi2t<-3.4*t-0.02^(2.3)

s21t<- a1t*cos(2*pi*phi1t)
s22t<-a2t*cos(2*pi*phi2t)

s2t<-s21t+s22t

# Ahora consideremos las siguientes dos funciones de tendencia

t1t<-8*(1/(1+(t/5)^2)+exp(-t/10))
t2t<-2*t+10*exp(-(t-4)^2/6)

# Armamos Xarma1 y Y0
Xarma1<-arima.sim(list(order=c(1,0,1), ar=c(-.5), ma=c(-.4)), n=n, rand.gen = rt, df=4)
sigmat<-1+.1*cos(pi*t)
X1<-2*sigmat*Xarma1

Y0<-s1t+t1t+X1

try <- emd(Y0, t, boundary="wave")

# Figura 2.1 (no guarde la semilla, con lo cual no va a ser exacta)
par(mfrow=c(1, 1), mar=c(3,1,3,1))
plot(t, try$imf[,4], type="l", xlab="", ylab="", ylim=c(-6,6),
     xaxt = "n", yaxt = "n"); abline(h=0)
