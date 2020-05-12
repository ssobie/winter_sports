##Find snow course information

read_course_obs <- function(site) {

   obs.file <- paste('/storage/data/projects/rci/data/winter_sports/obs/snow_courses/',site,'_snow_course.csv',sep='')
   obs.data <- read.csv(obs.file,header=F,as.is=T)
   obs.dates <- format(as.Date(obs.data[,1]),'%Y-%m-%d')
   obs.swe <- obs.data[,3] ##mm
   obs.na <- is.na(obs.swe)
   obs.swe <- obs.swe[!obs.na]
   obs.dates <- as.Date(obs.dates[!obs.na])

   rv <- list(dates=obs.dates,swe=obs.swe)
   return(rv)
}


cal.courses <- c('brookmere','callaghan','dickson_lake','disappointment_lake','dog_mountain','duffey_lake',
                    'gnawed_mountain','grouse_mountain','hamilton_hill','highland_valley',
                    'klesilkwa','lightning_lake','mcgillivray_pass','nahatlatch','orchid_lake','palisade_lake',
                    'shovelnose_mountain','stave_lake','sumallo_river_west','wahleach')
eval.courses <- c('blackwall_peak_course','boston_bar_lower','boston_bar_upper','burwell_lake',
                  'chapman_creek','cornwall_hills','diamond_head','edwards_lake','garibaldi_lake',
                  'hollyburn','hope',
                  'loch_lomond','lytton','mount_seymour','new_tashme',
                  'ottomite','pavilion_mountain',##'shalalth',
                      'sumallo_river','tenquille_course','whistler_mountain','wolverine_creek' )


sites <- sort(c(cal.courses,eval.courses))

site.info <- matrix(0,nrow=length(sites),ncol=4)

for (s in seq_along(sites)) {
   site <- sites[s]

   site.data <- read_course_obs(site)
   site.info[s,1] <- site
   site.info[s,2] <- length(site.data$swe)     
   site.info[s,3] <- format(site.data$dates[1],'%Y')
   site.info[s,4] <- format(tail(site.data$dates,1),'%Y')
}

write.table(site.info,file='/storage/data/projects/rci/data/winter_sports/snow_course_counts.csv',sep=',',quote=F,row.name=F,col.name=F)