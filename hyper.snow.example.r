
##Simple example of snow phase transition
library(scales)

coeffs <- list(a=-49.49,b=0.5628,c=1.5,d=1.0209)
tas <- seq(-10,10,0.05)
frac <- coeffs$a*(tanh(coeffs$b*(tas-coeffs$c))-coeffs$d)

x <- seq(-5,7,0.01)
frac.poly <- coeffs$a*(tanh(coeffs$b*(x-coeffs$c))-coeffs$d)

plot.file <- '/storage/data/projects/rci/data/winter_sports/ncc_2019/plots/hyper_snow_phase.2019.png'

png(file=plot.file,width=5,height=2,units='in',res=600,pointsize=6,bg='white')
par(mar=c(5,5,2,1))
plot(tas,frac,type='l',lwd=2,xlab="Daily Average Temperature (\u00B0C)",ylab='Percent Snow (%)',
              xaxs='i',yaxs='i',cex.lab=2,cex.axis=2)

polygon(x=c(x,rev(x)),y=c(rep(0,length(x)),rev(frac.poly)),col=alpha('gray',0.25))
polygon(x=c(x,rev(x)),y=c(rep(100,length(x)),rev(frac.poly)),col=alpha('blue',0.25))

rect(-10,0,-5,100,col=alpha('gray',0.25))
##rect (-5,0, 7,100,col=alpha('blue',0.25))
rect(  7,0,10,100,col=alpha('blue',0.25))

box(which='plot')

dev.off()