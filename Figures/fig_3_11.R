library('SynchWave')

t <- seq(0,10*pi, (1/30))
s <- function(t) sin(t)
st <- s(t)

dt <- t[2]-t[1]
nv <- 32
opt <- list(type = "mhat")

# transformada wavelet continua 
cwtfit <- cwt_fw(st, opt$type, nv, dt, opt)

# Set plot layout
par(mai=c(0.5,1,0.5,1))
layout(mat = matrix(c(1, 2), 
                    nrow = 2, 
                    ncol = 1),
       heights = c(1, 2),    # Heights of the two rows
       widths = c(3))     # Widths of the two columns

plot(t, st, type="l", xlab='', ylab='')

R1 = Re(cwtfit$Wx)
image.plot(list(x=t, y=cwtfit$asc, z=t(R1)), log="y",
           xlab="", ylab="", main="",
           col=designer.colors(64, c("grey100", "grey65", "grey35", "grey0")),
           ylim=c(7, 0.2), xlim=c(3, 27))

# La imagen presente en la tesis fue posteriormente retocada en paint 