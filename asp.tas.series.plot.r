##Script to plot temperature series comparisons of the snow pillows vs ERA GCM-PRISM simulations

library(zoo)

sites <- c('blackwall_peak_pillow',
           'spuzzum_creek',
           'upper_squamish',
           'chilliwack_river',
           'tenquille_lake',
           'wahleach_lake')
           
site.names <- c('Blackwall Peak','Spuzzum Creek','Upper Squamish','Chilliwack River','Tenquille Lake','Wahleach_Lake')
site.model <- c('PNWNAmet','PNWNAmet','PNWNAmet','PNWNAmet','PNWNAmet','PNWNAmet')
site.xlims <-         rbind(as.Date(c('1967-01-01','1979-12-31')),
                            as.Date(c('1998-01-01','2011-12-31')),
                            as.Date(c('1989-01-01','2012-12-31')),
                            as.Date(c('1991-01-01','2012-12-31')),
                            as.Date(c('2000-01-01','2012-12-31')),
                            as.Date(c('2005-01-01','2012-12-31')))
                            
plot.dir <- '/storage/data/projects/rci/data/winter_sports/plots/'
plot.file <- paste0(plot.dir,'asp.era5.tas.series.png')
png(file=plot.file,width=8,height=8,units='in',res=600,pointsize=6,bg='white')

par(mfrow=c(3,2))
par(mar=c(3.1,5.1,1.1,2.1))
for (i in seq_along(sites)) {

    site <- sites[i]
    print(site)

    era.file <- paste0('/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/snow/snow_sites/',site,'_PNWNAmet_800m_data.csv')
    era.data <- read.csv(era.file,header=T,as.is=T)
    era.pr <- era.data$Pr
    era.tasmax <- era.data$Tasmax
    era.tasmin <- era.data$Tasmin
    era.tas <- era.data$Tas
    era.dates <- era.data$Dates

    ##Snow Pillow Data
    pillow.file <- paste('/storage/data/projects/rci/data/winter_sports/obs/snow_pillow/',site,'.csv',sep='')
    pillow.data <- read.csv(pillow.file,header=T,as.is=T)
    pillow.dates <- format(as.Date(pillow.data[,2]),'%Y-%m-%d')
    pillow.tasmax <- pillow.data[,3]
    pillow.tasmin <- pillow.data[,5]
    pillow.tas <- (pillow.tasmax + pillow.tasmin)/2
    pillow.precip <- pillow.data[,7]##mm
    pillow.swe <- pillow.data[,11] ##mm
    pillow.pack <- pillow.data[,13] ##cm

    ix <- era.dates %in% pillow.dates
    
    era.slct <- era.tas[ix]
    asp.slct <- pillow.tas
    flag <- is.na(asp.slct) 

    plot(as.Date(pillow.dates[!flag]),era.slct[!flag]-asp.slct[!flag],type='l',
         xlim=site.xlims[i,],ylim=c(-15,15),  
         main='',col='gray',xaxs='i',xlab='',ylab='Mean Temp. Diff. (\u00B0C)',cex.lab=2.25,cex.axis=2.25)
    if (i==1) {
       axis(1,at=as.Date(c('1968-01-01','1970-01-01','1972-01-01','1974-01-01','1976-01-01','1978-01-01')),
              label=c('1968','1970','1972','1974','1976','1978'),cex.axis=2.25)
    }
    lines(as.Date(pillow.dates[!flag]),rollmean(era.slct[!flag]-asp.slct[!flag],7,fill='extend'),lwd=1,col='red')
    abline(v=as.Date(pillow.dates[flag]),col='white')
    bias <- round(mean(era.slct[!flag]-asp.slct[!flag]),2)
    text(site.xlims[i,1]+0.4*diff(site.xlims[i,]),14,paste0(site.names[i],' (',bias,' \u00B0C)'),cex=2.5)
    
    abline(h=0)
    box(which='plot',lwd=2)
    if (i==6) {
       legend('bottomright',legend=c('Daily (Model-Obs.)','7-Day Smoothed'),col=c('gray','red'),pch=15,cex=2.25)
    }

}

dev.off()