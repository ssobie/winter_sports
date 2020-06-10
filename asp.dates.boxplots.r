##SCript to compute the start and end dates and lengths of the snow season for each year

source('/storage/home/ssobie/code/repos/winter_sports/site.coordinates.r',chdir=T)

library(ncdf4)
library(plotrix)

get_snow_season_dates <- function(input.snow,dates) {

  years <- format(as.Date(dates),'%Y')
  jdays <- format(as.Date(dates),'%j')
  jdays.list <- tapply(jdays,as.factor(years),list)
  peaks <- tapply(input.snow,as.factor(years),function(x){which.max(x[1:300])})                      
                    
  sim.snow <- input.snow
  sim.snow[sim.snow>0] <- 100
                      
  patterns <- rle(sim.snow)
  pt     <- patterns$lengths
  pt.sum <- cumsum(patterns$lengths)
  all.zeroes <- which(patterns$values==100)
  zeroes <- all.zeroes[patterns$lengths[all.zeroes] > 100]
  top.ix <- pt.sum[zeroes]
  bot.ix <- pt[zeroes]

  pt.diff <- c(pt[1],top.ix[1], top.ix[-1] - bot.ix[-1], tail(top.ix,1))
  ix <- sort(unique(c(pt.diff,top.ix))) 
  if (ix[1] > 100) {
     ix <- sort(unique(c(1,pt.diff,top.ix)))
  }   
  ##plot(as.Date(dates),input.snow,xlim=c(as.Date('1967-01-01'),as.Date('1980-01-01')))
  ##abline(v=as.Date(dates[ix]),col='red')
  len <- length(ix)

  lengths <- ix[seq(2,len,2)] - ix[seq(1,len-1,2)]
  starts <- as.numeric(format(as.Date(dates[ix[seq(1,len-1,2)]]),'%j'))
  ends <- as.numeric(format(as.Date(dates[ix[seq(2,len,2)]]),'%j'))

  rv <- list(lengths=lengths,
             starts =starts,
             peaks = peaks,
             ends = ends)

  return(rv)               

}

##Slope and Aspect Values

model <- 'ERA5'
type <- 'SWE'

sites <- c('blackwall_peak_pillow','spuzzum_creek','upper_squamish','chilliwack_river','tenquille_lake','wahleach_lake')
site.names <- c('Blackwall\nPeak','Spuzzum\nCreek','Upper\nSquamish','Chilliwack\nRiver','Tenquille\nLake','Wahleach\nLake')
site.locs <- seq(1-1/12,by=-1/6,length.out=6)
snow.file <- paste0('/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/snow_series/spuzzum_creek_',model,'_PRISM_TPS_snow_model_data.csv')
snow.data <- read.csv(snow.file,header=T,as.is=T)
snow.dates <- snow.data$Dates
snow.years <- unique(format(as.Date(snow.dates),'%Y'))
ylen <-  length(snow.years)

slen <- length(sites)
sim.lengths <- matrix(0,nrow=slen,ncol=ylen)
sim.starts <- matrix(0,nrow=slen,ncol=ylen)
sim.peaks <- matrix(0,nrow=slen,ncol=ylen)
sim.ends <- matrix(0,nrow=slen,ncol=ylen)

results.matrix <- matrix(0,nrow=20,ncol=length(sites))


##Create a plot to display the season dates
##Start on September 1st


plot.dir <- '/storage/data/projects/rci/data/winter_sports/plots/'
plot.file <- 'asp.start.end.boxplots.png'

png(file=paste0(plot.dir,plot.file),width=4,height=3,units='in',res=600,pointsize=6,bg='white')
par(mfrow=c(6,1))
par(mar=c(0,0,0,0),oma=c(6,10,5,3))

##plot(c(),xlim=c(0,366),ylim=c(0,2),xlab='',ylab='',
##         cex.lab=1.5,cex.axis=1.5,axes=F)
##axis(1,at=seq(0,366,1),label=seq(0,366,1),cex.axis=1.5)
##axis(2,at=seq(1,6,1),label=site.names,cex.axis=1.5)

##Loop over sites
for (i in seq_along(sites)) {

    site <- sites[i]
    print(site)
    snow.file <- paste0('/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/snow_series/',site,'_',model,'_PRISM_TPS_snow_model_data.csv')
    snow.data <- read.csv(snow.file,header=T,as.is=T)
    snow.dates <- snow.data$Dates
    snow.swe <- snow.data$SWE*1000
    snow.seasons <- get_snow_season_dates(snow.swe,snow.dates)

    sim.lengths <- snow.seasons$lengths
    sim.starts <- snow.seasons$starts
    sim.peaks <- snow.seasons$peaks
    sim.ends <- snow.seasons$ends

    sim.lengths[sim.lengths==0] <- NA
    sim.starts[sim.starts==0] <- NA
    sim.peaks[sim.peaks==0] <- NA
    sim.ends[sim.ends==0] <- NA

    dates.file <- paste0('/storage/data/projects/rci/data/winter_sports/plots/snow_seasons/',site,'_season_dates.csv') 
    season.dates <- read.csv(dates.file,header=T,as.is=T)
    pillow.starts <- as.numeric(format(as.Date(season.dates$Start),'%j'))
    pillow.peaks <- as.numeric(format(as.Date(season.dates$Peak),'%j'))
    pillow.ends <- as.numeric(format(as.Date(season.dates$End),'%j'))
    pillow.lengths <- as.numeric(season.dates$Length)

    snow.match <- as.numeric(snow.years) %in% as.numeric((format(as.Date(season.dates$Start),'%Y')))
    pillow.match <- as.numeric((format(as.Date(season.dates$Start),'%Y'))) %in% as.numeric(snow.years) 
    
    lengths.diff <- season.dates$Length[pillow.match] - sim.lengths[snow.match]
    starts.diff <-  pillow.starts[pillow.match] - sim.starts[snow.match]
    peaks.diff <-  pillow.peaks[pillow.match] - sim.peaks[snow.match]
    ends.diff <-  pillow.ends[pillow.match] - sim.ends[snow.match]

    rv <- c(mean(pillow.starts[pillow.match]),mean(pillow.peaks[pillow.match]),mean(pillow.ends[pillow.match]),mean(pillow.lengths[pillow.match]),
            mean(sim.starts[snow.match]),mean(sim.peaks[snow.match]),mean(sim.ends[snow.match]),mean(sim.lengths[snow.match]),
            sd(sim.starts[snow.match]),sd(sim.peaks[snow.match]),sd(sim.ends[snow.match]),sd(sim.lengths[snow.match]),
            mean(starts.diff),mean(peaks.diff),mean(ends.diff),mean(lengths.diff),
            sd(starts.diff),sd(peaks.diff),sd(ends.diff),sd(lengths.diff))
    results.matrix[,i] <- round(rv,1)


    starts.df <- data.frame(pillow=pillow.starts[pillow.match],
                            sim=sim.starts[snow.match])
    if (i==1) {                                                        
       starts.df <- data.frame(pillow=pillow.starts[pillow.match][-1],
                               sim=sim.starts[snow.match][-1])
    }
    peaks.df <- data.frame(pillow=pillow.peaks[pillow.match]+121,
                            sim=sim.peaks[snow.match]+121)
    ends.df <- data.frame(pillow=pillow.ends[pillow.match]+121,
                          sim=sim.ends[snow.match]+121)

    boxplot(starts.df$pillow-244,starts.df$sim-244,at=c(1.5,0.5),names=c('Pillow','Model'),las=2,boxwex=0.5,
            col=c('gray','lightblue'),horizontal=TRUE,axes=FALSE,ylim=c(1,366),xlim=c(0,2),yaxs='i')
    abline(v=c(1,31,62,92,123,154,182,213,243,274,304,335),col='lightgray')            
    if (i==1) {
        text(y=1.5,x=350,'ASP',cex=1.5)
        text(y=0.5,x=350,'Model',cex=1.5)
    }
    boxplot(starts.df$pillow-244,starts.df$sim-244,at=c(1.5,0.5),boxwex=0.5,
            col=c('gray','lightblue'),horizontal=TRUE,axes=FALSE,add=TRUE)
    mtext(site.names[i],at=site.locs[i],las=2,
          side=2,outer=TRUE,cex=1.25,line=1)

    boxplot(peaks.df$pillow,peaks.df$sim,at=c(1.5,0.5),names=c('pillow','sim'),las=2,boxwex=0.5,
            col=c('gray','lightblue'),horizontal=TRUE,axes=FALSE,add=TRUE)
    boxplot(ends.df$pillow,ends.df$sim,at=c(1.5,0.5),names=c('pillow','sim'),las=2,boxwex=0.5,
            col=c('gray','lightblue'),horizontal=TRUE,axes=FALSE,add=TRUE)
    if (i==6) {
       axis(1,at=c(1,31,62,92,123,154,182,213,243,274,304,335),
            label=c('Sep','Oct','Nov','Dec','Jan','Feb','Mar','Apr','May','Jun','Jul','Aug'),
            cex.axis=1.5)}
    box(which='plot')

    
}

mtext('Date',side=1,outer=TRUE,cex=1.25,line=3.25)
##mtext(site.names,at=c(1:6),las=2,
##       side=2,outer=TRUE,cex=1.25,line=1)

mtext(c('Starts','Peaks','Ends'),at=c(0.1,0.62,0.84),
       side=3,outer=TRUE,cex=1.25,line=1)

dev.off()
    results.matrix <- cbind(c('Site','ASP Start','ASP Peak','ASP End','ASP Length',
                              'Sim Start','Sim Peak','Sim End','Sim Length',
                              'Start SD' ,'Peak SD','End SD' ,'Length SD' ,
                              'Start Diff','Peak Diff','End Diff','Length Diff',      
                              'SDiff SD'  ,'PDiff SD','EDiff SD','LDiff SD'),rbind(site.names,results.matrix))

##write.table(results.matrix,file='/storage/data/projects/rci/data/winter_sports/asp.snow.seas.pnw.data.csv',
##            sep=',',quote=FALSE,row.name=F,col.name=F)

