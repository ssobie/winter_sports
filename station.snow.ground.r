##Script to evaluate YVR temperature stats to deal with
##the EPW file year selection issues


##-------------------------------------------------------------------
##Read EC downloaded data

read_ec_obs_file <- function(var.name,interval,dir) {

   prefix <- 'eng-daily-'
   ec.file <- paste0(dir,prefix,interval,'.csv')
   ec.raw <- read.csv(ec.file,skip=24,header=T)

   ix <- grep(var.name,names(ec.raw))
  
   data <- ec.raw[,ix]
   dates <- as.Date(paste(ec.raw[,2],sprintf('%02d',ec.raw[,3]),sprintf('%02d',ec.raw[,4]),sep='-'))

   rv <- list(data=data,dates=dates)
   return(rv)
}

all_years_from_ec_obs <- function(var.title,dir) {

    ec.2004 <- read_ec_obs_file(var.title,'01012004-12312004',dir)    
    ec.2005 <- read_ec_obs_file(var.title,'01012005-12312005',dir)    
    ec.2006 <- read_ec_obs_file(var.title,'01012006-12312006',dir)    
    ec.2007 <- read_ec_obs_file(var.title,'01012007-12312007',dir)    
    ec.2008 <- read_ec_obs_file(var.title,'01012008-12312008',dir)    
    ec.2009 <- read_ec_obs_file(var.title,'01012009-12312009',dir)
    ec.2010 <- read_ec_obs_file(var.title,'01012010-12312010',dir)
    ec.2011 <- read_ec_obs_file(var.title,'01012011-12312011',dir)
    ec.2012 <- read_ec_obs_file(var.title,'01012012-12312012',dir)
    ec.2013 <- read_ec_obs_file(var.title,'01012013-12312013',dir)
    ec.2014 <- read_ec_obs_file(var.title,'01012014-12312014',dir)
    ec.2015 <- read_ec_obs_file(var.title,'01012015-12312015',dir)
    ec.2016 <- read_ec_obs_file(var.title,'01012016-12312016',dir)
    ec.2017 <- read_ec_obs_file(var.title,'01012017-12312017',dir)
    ec.2018 <- read_ec_obs_file(var.title,'01012018-12312018',dir)

    ec.data <- c(ec.2004$data,ec.2005$data,
                 ec.2006$data,ec.2007$data,ec.2008$data,ec.2009$data,
                 ec.2010$data,ec.2011$data,ec.2012$data,ec.2013$data,
                 ec.2014$data,ec.2015$data,ec.2016$data,ec.2017$data,
                 ec.2018$data)
    ec.dates <- c(ec.2004$dates,ec.2005$dates,
                 ec.2006$dates,ec.2007$dates,ec.2008$dates,ec.2009$dates,
                 ec.2010$dates,ec.2011$dates,ec.2012$dates,ec.2013$dates,
                 ec.2014$dates,ec.2015$dates,ec.2016$dates,ec.2017$dates,
                 ec.2018$dates)

    rv <- list(data=ec.data,dates=ec.dates)
    return(rv)
}

##-------------------------------------------------------------------
##Read CDCD file

read_cdcd_file <- function(var.name,dir) {
   file <- paste0(dir,'1048898_',var.name,'.csv')
   cdcd <- read.csv(file,header=T)
   years <- cdcd[,1]
   jdays <- cdcd[,2]
   dates <- as.Date(paste(years,jdays,sep='-'),'%Y-%j')
   
   rv <- list(data=cdcd[,3],dates=dates)
   return(rv)
}

##-------------------------------------------------------------------
##Create time series of variable


##-------------------------------------------------------------------
##

whi.dir <- '/storage/data/projects/rci/data/winter_sports/obs/Whistler/'

snow.cdcd <- read_cdcd_file('SNOW_ON_THE_GROUND',whi.dir)

snow.ec <- all_years_from_ec_obs('Snow.on.Grnd..cm.',whi.dir)


##snow.series <- c(snow.cdcd$data[1:5130],snow.ec$data)
##snow.dates <- c(snow.cdcd$dates[1:5130],snow.ec$dates)

snow.series <- c(snow.cdcd$data[1:8380],snow.ec$data)
snow.dates <- c(snow.cdcd$dates[1:8380],snow.ec$dates)


snow.yr.fac <- as.factor(format(snow.dates,'%Y'))
snow.mn.fac <- as.factor(format(snow.dates,'%m'))
seasons <-  c("DJF", "DJF", "MAM", "MAM", "MAM", "JJA", "JJA", "JJA", "SON", "SON", "SON", "DJF")
seasonal.fac <- factor(seasons[snow.mn.fac], levels=c('DJF', 'MAM', 'JJA', 'SON'))
avg.fac <- list(snow.yr.fac,seasonal.fac)

snow.peak <- tapply(snow.series,snow.yr.fac,max,na.rm=T)
snow.seas.mean <- tapply(snow.series,avg.fac,mean,na.rm=T)
snow.days <- snow.series < 5
snow.seas.sums <- tapply(snow.days,avg.fac,sum,na.rm=T)

snow.djf <- snow.seas.mean[,1]

x <- 1980:2018

