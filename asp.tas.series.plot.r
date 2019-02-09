##Script to plot temperature series comparisons of the snow pillows vs ERA GCM-PRISM simulations

library(zoo)

sites <- c('spuzzum_creek',
           'upper_squamish',
           'chilliwack_river',
           'tenquille_lake')
           
site.names <- c('Spuzzum Creek','Upper Squamish','Chilliwack River','Tenquille Lake')

plot.dir <- '/storage/data/projects/rci/data/winter_sports/plots/'
png(file=paste0(plot.dir,'asp.ncep2.tas.series.png'),width=1200,height=1200)

par(mfrow=c(4,1))
par(mar=c(5.1,5.1,1.1,2.1))
for (i in seq_along(sites)) {
    site <- sites[i]
    print(site)

    era.file <- paste0('/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/snow/snow_sites/',site,'_NCEP2_800m_data.csv')
    era.data <- read.csv(era.file,header=T,as.is=T)
    era.pr <- era.data$Pr
    era.tasmax <- era.data$Tasmax
    era.tasmin <- era.data$Tasmin
    era.tas <- era.data$Tas
    era.dates <- era.data$Dates

    ##Snow Pillow Data
    pillow.file <- paste('/storage/data/projects/rci/data/assessments/snow_model/snow_pillow/',site,'_asp.csv',sep='')
    pillow.data <- read.csv(pillow.file,header=T,as.is=T)
    pillow.dates <- format(as.Date(pillow.data[,2]),'%Y-%m-%d')
    pillow.tasmax <- pillow.data[,3]
    pillow.tasmin <- pillow.data[,5]
    pillow.tas <- (pillow.tasmax + pillow.tasmin)/2
    pillow.precip <- pillow.data[,7]##mm
    pillow.swe <- pillow.data[,11] ##mm
    pillow.pack <- pillow.data[,13] ##cm

    ix <- era.dates %in% pillow.dates
    
    timex <- 1:4748
    era.slct <- era.tas[ix][timex]
    asp.slct <- pillow.tas[timex]
    flag <- is.na(asp.slct) 

    plot(as.Date(pillow.dates[timex][!flag]),era.slct[!flag]-asp.slct[!flag],type='l',
         xlim=c(as.Date('1992-06-01'),as.Date('2011-01-01')),ylim=c(-10,10),  
         main='',col='gray',xaxs='i',xlab='Date',ylab='Mean Temp. (\u00B0C)',cex.lab=1.75,cex.axis=1.75)
    lines(as.Date(pillow.dates[timex][!flag]),rollmean(era.slct[!flag]-asp.slct[!flag],7,fill='extend'),lwd=1)                
    abline(v=as.Date(pillow.dates[timex][flag]),col='white')
    bias <- round(mean(era.slct[!flag]-asp.slct[!flag]),2)
    text(as.Date('1995-01-01'),8,paste0(site.names[i],' (',bias,' \u00B0C)'),cex=2.5)
    
    abline(h=0)
    box(which='plot',lwd=2)
    if (i==4) {
       legend('bottomleft',legend=c('Daily (Model-Obs.)','7-Day Smoothed'),col=c('gray','black'),pch=15,cex=1.75)
    }

}

dev.off()