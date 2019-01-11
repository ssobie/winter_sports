source('/storage/data/projects/rci/assessments/code/eco_hydrology.R')

snow.melt <- function(precip_mm, Tmax_C, Tmin_C, Date, lat_deg,slope=0, aspect=0, tempHt=1, windHt=2, groundAlbedo=0.25,
 		SurfEmissiv=0.95, windSp=1, forest=0, startingSnowDepth_m=0, startingSnowDensity_kg_m3=450){

  Tmax.check <- Tmax_C <  Tmin_C
  Tmax_C[Tmax.check] <- Tmin_C[Tmax.check]+1
      
   ## Constants :
	WaterDens <- 1000			# kg/m3
        SnowDens <- 600                         # kg/m3
        lambda <- 3.35*10^5			# latent heat of fusion (kJ/m3)
	lambdaV <- 2500				# (kJ/kg) latent heat of vaporization
	SnowHeatCap <- 2.1			# kJ/kg/C
	LatHeatFreez <- 333.3		        # kJ/kg
	Cw <- 4.2*10^3				# Heat Capacity of Water (kJ/m3/C)
	
##	Converted Inputs :
	Tav <- (Tmax_C+Tmin_C)/2		# degrees C
	precip_m <- precip_mm*0.001	 	# precip in m 
	R_m <- precip_m					# (m) depth of rain
	R_m[which(Tav < 0)] <- 0		# ASSUMES ALL SNOW at < 0C
	NewSnowDensity <- 50+3.4*(Tav+15)		# kg/m3
	NewSnowDensity[which(NewSnowDensity < 50)] <- 50
	NewSnowWatEq <- precip_m				# m
	NewSnowWatEq[which(Tav >= 0)] <- 0			# No new snow if average temp above or equals 0 C
        ##--------------------
        if (1==0) {
         tas.pos <- Tav
         tas.pos[tas.pos<0]<-0
         frac <- exp(-tas.pos)
         model.data <- precip_m
         model.data[Tav > 10] <- 0
         warm <- Tav <= 10 & Tav >=0
         model.data[warm] <- precip_m[warm]*frac[warm]
         NewSnowWatEq <- model.data
         R_m <- precip_m
         R_m[Tav < 0] <- 0
         R_m[warm] <- precip_m[warm]*(1-frac[warm]) 
 
         if (!is.null(Snow)) {
           NewSnowWatEq <- Snow$snow
           R_m <- Snow$rain
         }
        }
        ##------------ 
        tas <- Tav               # degrees C
        coeffs <- list(a=-49.49,b=0.5628,c=1.5,d=1.0209)
        frac <- coeffs$a*(tanh(coeffs$b*(tas-coeffs$c))-coeffs$d)
        sample <- runif(length(tas),min=0,max=100)
        test <- sample > frac
        high.temp <- tas > 12
        test[high.temp] <- TRUE
        snow.type <- rep(TRUE,length(tas))
        snow.type[test] <- FALSE

        NewSnowWatEq <- precip_m
        NewSnowWatEq[!snow.type] <- 0
        R_m <- precip_m
        R_m[snow.type] <- 0

        ##-----------------------------

        NewSnow <- NewSnowWatEq*WaterDens/NewSnowDensity		# m
	JDay <- strptime(Date, format="%Y-%m-%d")$yday+1
	lat <- lat_deg*pi/180		#	latitude in radians
	##rh 	<- log((windHt+0.001)/0.001)*log((tempHt+0.0002)/0.0002)/(0.41*0.41*windSp*86400)	# (day/m) Thermal Resistance	 
        rh      <- log((10+0.001)/0.001)*log((10+0.0002)/0.0002)/(0.41*0.41*windSp*86400)       # (day/m) Thermal Resistance
	if (length(windSp)==1) rh <- rep(rh,length(precip_mm))									##	creates a vector of rh values
	cloudiness 		<- EstCloudiness(Tmax_C,Tmin_C)
	AE 				<- AtmosphericEmissivity(Tav, cloudiness)	# (-) Atmospheric Emissivity

#  New Variables	:
	SnowTemp 		<- rep(0,length(precip_m)) 		# Degrees C
	rhos 			<- SatVaporDensity(SnowTemp)	# 	vapor density at surface (kg/m3)
	rhoa 			<- SatVaporDensity(Tmin_C)		#	vapor density of atmoshpere (kg/m3) 
	SnowWaterEq 	<- vector(length=length(precip_mm))		#  (m) Equiv depth of water
	TE 				<- rep(SurfEmissiv,length(precip_mm))	#	(-) Terrestrial Emissivity
	DCoef 			<- rep(0,length(precip_mm))				#   Density Coefficient (-) (Simplified version)
	SnowDensity 	<- rep(SnowDens,length(precip_mm))			#  (kg/m3)  Max density is 500
	SnowDepth 		<- vector(length=length(precip_mm))		#  (m)
        OldSnowDepth            <- vector(length=length(precip_mm))             #  (m)
	SnowMelt 		<- rep(0,length(precip_mm))				#  (m)
        melt            <- rep(0,length(precip_mm))                             #  (m)
        DensityPerc             <- rep(0,length(precip_mm))
        Albedo 			<- rep(groundAlbedo,length(precip_mm)) 	#  (-) This will change for days with snow
	
##	Energy Terms
	H 		<- vector(length=length(precip_mm))	#	Sensible Heat exchanged (kJ/m2/d)
	E 		<- vector(length=length(precip_mm))	#	Vapor Energy	(kJ/m2/d)
	S 		<- vector(length=length(precip_mm))	#	Solar Radiation (kJ/m2/d)
	La 		<- Longwave(AE, Tav)					#	Atmospheric Longwave Radiation (kJ/m2/d)
	Lt 		<- vector(length=length(precip_mm))	#	Terrestrial Longwave Radiation (kJ/m2/d)
	G 		<- 173								#	Ground Condution (kJ/m2/d) 
	P 		<- Cw * R_m * Tav					# 	Precipitation Heat (kJ/m2/d)
	Energy 	<- vector(length=length(precip_mm))	# Net Energy (kJ/m2/d)

##  Initial Values.  
	SnowWaterEq[1] 	<- startingSnowDepth_m * startingSnowDensity_kg_m3 / WaterDens		
	SnowDepth[1] 	<- startingSnowDepth_m			
        OldSnowDepth[1]         <- startingSnowDepth_m
	Albedo[1] <- ifelse(NewSnow[1] > 0, 0.98-(0.98-0.50)*exp(-4*NewSnow[1]*10),ifelse(startingSnowDepth_m == 0, groundAlbedo, max(groundAlbedo, 0.5+(groundAlbedo-0.85)/10)))  # If snow on the ground or new snow, assume Albedo yesterday was 0.5
	S[1] <- Solar(lat=lat,Jday=JDay[1], Tx=Tmax_C[1], Tn=Tmin_C[1], albedo=Albedo[1], forest=forest, aspect=aspect, slope=slope, latUnits='degrees')
	H[1] <- 1.29*(Tav[1]-SnowTemp[1])/rh[1] 
	E[1] <- lambdaV*(rhoa[1]-rhos[1])/rh[1]
	if(startingSnowDepth_m>0) TE[1] <- 0.97 
	Lt[1] <- Longwave(TE[1],SnowTemp[1])
	Energy[1] <- S[1] + La[1] - Lt[1] + H[1] + E[1] + G + P[1]
	SnowDensity[1] <- ifelse((startingSnowDepth_m+NewSnow[1])>0, min(SnowDens, (startingSnowDensity_kg_m3*startingSnowDepth_m + NewSnowDensity[1]*NewSnow[1])/(startingSnowDepth_m+NewSnow[1])), SnowDens)
	SnowMelt[1] <- max(0,	min((startingSnowDepth_m/10+NewSnowWatEq[1]),  # yesterday on ground + today new  
				      (Energy[1]-SnowHeatCap*(startingSnowDepth_m/10+NewSnowWatEq[1])*WaterDens*(0-SnowTemp[1]))/(LatHeatFreez*WaterDens) ) )
        melt[1] <- 0
	SnowDepth[1] <- max(0,(startingSnowDepth_m/10 + NewSnowWatEq[1]-SnowMelt[1])*WaterDens/SnowDensity[1])
	SnowWaterEq[1] <- max(0,startingSnowDepth_m/10-SnowMelt[1]+NewSnowWatEq[1])	
        DensityPerc[1] <- 10
	
##  Snow Melt Loop	
	for (i in 2:length(precip_m)){
		if (NewSnow[i] > 0){ 
			Albedo[i] <- 0.98-(0.98-Albedo[i-1])*exp(-4*NewSnow[i]*10)
		} else if (SnowDepth[i-1] < 0.1){ 
			Albedo[i] <- max(groundAlbedo, Albedo[i-1]+(groundAlbedo-0.85)/10)
		} else Albedo[i] <- 0.35-(0.35-0.98)*exp(-1*(0.177+(log((-0.3+0.98)/(Albedo[i-1]-0.3)))^2.16)^0.46)

		S[i] <- Solar(lat=lat,Jday=JDay[i], Tx=Tmax_C[i], Tn=Tmin_C[i], albedo=Albedo[i-1], forest=forest, aspect=aspect, slope=slope, printWarn=FALSE)

		if(SnowDepth[i-1] > 0) TE[i] <- 0.97 	#	(-) Terrestrial Emissivity
		if(SnowWaterEq[i-1] > 0 | NewSnowWatEq[i] > 0) {
			DCoef[i] <- 9.45  ##DCoef.val ##9.5
			if(SnowMelt[i-1] == 0){ 
				SnowTemp[i] <- max(min(0,Tmin_C[i]),min(0,(SnowTemp[i-1]+min(-SnowTemp[i-1],Energy[i-1]/((SnowDensity[i-1]*
					SnowDepth[i-1]+NewSnow[i]*NewSnowDensity[i])*SnowHeatCap*1000)))))
			}
		}

		rhos[i] <- SatVaporDensity(SnowTemp[i])
		H[i] <- 1.29*(Tav[i]-SnowTemp[i])/rh[i] 
		E[i] <- lambdaV*(rhoa[i]-rhos[i])/rh[i]
		Lt[i] <- Longwave(TE[i],SnowTemp[i])
		Energy[i] <- S[i] + La[i] - Lt[i] + H[i] + E[i] + G + P[i]
                ##print(i)
                ##print(Energy[i])
                if (is.na(Energy[i]))
                  browser()
		if (Energy[i]>0) k <- 1 else k <- 0.5
		
		SnowDensity[i] <- ifelse((SnowDepth[i-1]+NewSnow[i])>0, min(SnowDens, 
			((SnowDensity[i-1]+k*SnowDens*(SnowDens-SnowDensity[i-1])*exp(-DCoef[i]))*SnowDepth[i-1] + NewSnowDensity[i]*NewSnow[i])/(SnowDepth[i-1]+NewSnow[i])), SnowDens)

		SnowMelt[i] <- max(0,	min( (SnowWaterEq[i-1]+NewSnowWatEq[i]),  # yesterday on ground + today new
				      (Energy[i]-SnowHeatCap*(SnowWaterEq[i-1]+NewSnowWatEq[i])*WaterDens*(0-SnowTemp[i]))/(LatHeatFreez*WaterDens) )  )
                melt[i] <- SnowHeatCap*(SnowWaterEq[i-1]+NewSnowWatEq[i])*WaterDens*(0-SnowTemp[i])

                OldSnowDepth[i] <- max(0,(SnowWaterEq[i-1]+NewSnowWatEq[i]-SnowMelt[i])*WaterDens/SnowDensity[i])
		SnowDepth[i] <- max(0,(SnowWaterEq[i-1]+NewSnowWatEq[i]-SnowMelt[i])*WaterDens/SnowDensity[i])
		SnowWaterEq[i] <- max(0,SnowWaterEq[i-1]-SnowMelt[i]+NewSnowWatEq[i])	# (m) Equiv depth of water
                DensityPerc[i] <- SnowWaterEq[i]/SnowDepth[i]*100                

          if (grepl('19988-01-30',Date[i])) {
            print(i)
            print(Date[i])
            print(paste('New Snow (cm): ',NewSnow[i]*100,sep=''))
            print(paste('Snow Depth (cm): ',SnowDepth[i]*100,sep=''))
            print(paste('Old Snow Depth (cm): ',OldSnowDepth[i]*100,sep=''))
            print(paste('SWE (cm): ',SnowWaterEq[i-1]*WaterDens/SnowDensity[i]*100,sep=''))
            print(paste('New SWE (cm): ',NewSnowWatEq[i]*WaterDens/SnowDensity[i]*100,sep=''))
            print(paste('Snow Melt (cm): ',SnowMelt[i]*WaterDens/SnowDensity[i]*100,sep=''))
            print(paste('Snow Density: ',SnowDensity[i],sep=''))
            print(paste('Energy: ',Energy[i],sep=''))
            print(paste('Tmax: ',Tmax_C[i],sep=''))
            print(paste('Tmean: ',Tav[i],sep=''))
            print(paste('Tmin: ',Tmin_C[i],sep=''))
            print(paste('Snow Temp: ',SnowTemp[i],sep=''))
            print(paste('Solar: ',S[i],sep=''))
            print(paste('AtmLong: ',La[i],sep=''))
            print(paste('TerrLong: ',-1*Lt[i],sep=''))
            print(paste('Sensible: ',H[i],sep=''))
            print(paste('Vapour: ',E[i],sep=''))
            print(paste('PrecipHeat: ',P[i],sep=''))
            print(i)
            s <- i-30
            e <- i
            data <- (cbind(Date[s:e],S[s:e],La[s:e]-Lt[s:e],H[s:e],E[s:e],G,P[s:e]))

            par(mfrow=c(3,2))
            ##plot(as.Date(Date[s:e]),round(precip_mm[s:e],1),type='l',lwd=3,col='blue',main='Precip (mm)',cex.axis=1.5)
            plot(as.Date(Date[s:e]),Tmax_C[s:e]-Tmin_C[s:e],type='l',lwd=3,col='red',main='Tavg',cex.axis=1.5)
            plot(as.Date(Date[s:e]),round(NewSnow[s:e]*100,1),type='l',lwd=3,col='blue',main='New Snow (cm)',cex.axis=1.5)
            plot(as.Date(Date[s:e]),Tav[s:e],type='l',lwd=3,col='red',main='Tavg',cex.axis=1.5)
            plot(as.Date(Date[s:e]),SnowWaterEq[s:e]*1000,type='l',lwd=3,col='blue',main='SWE (mm)',cex.axis=1.5)
            plot(as.Date(Date[s:e]),SnowMelt[s:e]*100,type='l',lwd=3,col='orange',main='SnowMelt (cm)',cex.axis=1.5)
            a <- barplot(t(data[,2:7]),legend=c('Solar','Longwave','Sensible','Vapour','Ground','Precip'),cex.axis=1.5,
            args.legend=list(x='bottomleft'),col=c('red','orange','yellow','green','blue','purple'))
            axis(1,at=a,label=as.Date(Date[s:e]))
            box(which='plot')

            par(mfrow=c(4,2))
            plot(as.Date(Date[s:e]),round(Energy[s:e]-mean(Energy[s:e]),1),type='l',lwd=3,col='red',main='Energy',cex.axis=1.5)
            abline(h=0)
            plot(as.Date(Date[s:e]),round(S[s:e]-mean(S[s:e]),1),type='l',lwd=3,col='orange',main='Solar',cex.axis=1.5)
            abline(h=0)
            plot(as.Date(Date[s:e]),round(La[s:e]-mean(La[s:e]),1),type='l',lwd=3,col='red',main='Atm Long Wave',cex.axis=1.5)
            abline(h=0)
            plot(as.Date(Date[s:e]),round(-Lt[s:e]-mean(-Lt[s:e]),1),type='l',lwd=3,col='orange',main='Terr Long Wave',cex.axis=1.5)
            abline(h=0)
            plot(as.Date(Date[s:e]),round(La[s:e]-Lt[s:e],1),type='l',lwd=3,col='green',main='Long Wave Diff',cex.axis=1.5)
            abline(h=0)
            plot(as.Date(Date[s:e]),round(H[s:e]-mean(H[s:e]),1),type='l',lwd=3,col='orange',main='H',cex.axis=1.5)
            abline(h=0)
            plot(as.Date(Date[s:e]),round(E[s:e]-mean(E[s:e]),1),type='l',lwd=3,col='orange',main='E',cex.axis=1.5)
            abline(h=0)
            plot(as.Date(Date[s:e]),round(P[s:e]-mean(P[s:e]),1),type='l',lwd=3,col='orange',main='P',cex.axis=1.5)
            abline(h=0)
                                    
                       
            

            browser() 
          }

	}
	##Results<-data.frame(Date, Tmax_C, Tmin_C, precip_mm, R_m*1000, NewSnowWatEq,SnowMelt*100, NewSnow, SnowDepth, SnowWaterEq*1000,DensityPerc)
	##colnames(Results)<-c("Date", "MaxT_C", "MinT_C", "Precip_mm", "Rain_mm", "SnowfallWatEq_m", "SnowMelt_cm", "NewSnow_m",
        ##                     "SnowDepth_m", "SnowWaterEq_mm","DensityPerc")
        Results <- list(snowdepth=SnowDepth,swe=SnowWaterEq,snowfall=NewSnow,snowdense=DensityPerc,swe=SnowWaterEq)
	return(Results)
}
