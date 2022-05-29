t <- seq(-2,3, (1/30))

A <- function(t) -0.1*t^4+t^2+2
s <- function(t) A(t)*sin(20*t)

par(mfrow=c(1,1))
plot(t, s(t), type="l", xlab='', ylab='', xaxt = "n", yaxt = "n")
lines(t, A(t), col='red')
lines(t, -A(t), col='blue')
