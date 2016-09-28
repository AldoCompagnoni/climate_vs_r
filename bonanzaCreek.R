setwd("C:/Users/ac79/MEGA/Projects/LTER/")
#Functions
source("C:/Users/ac79/Documents/CODE/climate_vs_r/functions.R")


#DATA#####################################################################

#Meteorological data------------------------------------------------------
#FIRST 15 lines provide information related to the data set
meteo=read.csv("Data/BonanzaCreek/MCKINLEY PARK_Meteo.csv",skip=15) 

#Demographic data---------------------------------------------------------
#Beetle data, 
beetle1=read.table("Data/BonanzaCreek/35_BNZ_Beetles_Werner_1975-2012.txt",
                   header=T,sep=",")
#more detailed version of the above
beetle2=read.table("Data/BonanzaCreek/504_BNZ_Beetles_Kruse_1975-2013.txt",
              header=T,sep=",")
#Defoliating insects: samping twice each summer.
#samples are from 2 branches taken from 5 species
defoliate1=read.table("Data/BonanzaCreek/37_BNZ_Defoliating_Insects_Werner_1975_2012.txt",
                      header=T,sep=",")
#more detailed version of the above
defoliate2=read.table("Data/BonanzaCreek/505_BNZ_Defoliating_Insects_Kruse_1975_2014.csv")

#Snowshoe hare
hare=read.table("Data/BonanzaCreek/55_Hare_Data_2012.txt",sep=",",header=T)



#FORMAT DATA#####################################################################

#Meteorological quantities-------------------------------------------------
#substitute "S" and "T" (T=trace)
meteo$Min.Temperature[meteo$Min.Temperature=="S"]=NA
meteo$Max.Temperature[meteo$Max.Temperature=="S"]=NA
meteo$Precipitation=as.character(meteo$Precipitation)
meteo$Precipitation[which(meteo$Precipitation=="T")]="0"
meteo$Precipitation[which(meteo$Precipitation=="S")]=NA
#make fields numeric
meteo$Max.Temperature=as.numeric(as.character(meteo$Max.Temperature))
meteo$Min.Temperature=as.numeric(as.character(meteo$Min.Temperature))
meteo$Precipitation=as.numeric(meteo$Precipitation)
#derived field
meteo$Mean.Temp=(meteo$Min.Temperature+meteo$Max.Temperature)/2

#Climate summaries######################################
dates=formatMeteoDate(as.character(meteo$Day),c("year","month","day"),"-")
meteo=cbind(dates,meteo)

#Year
yrClim=summarise(group_by(meteo,year), yrTemp = mean(Mean.Temp,na.rm=T),
                 yrPrec = sum(Precipitation,na.rm=T))
#Growing season temp
gs=summarise(group_by(subset(meteo,month> 5 & month < 10),year), 
             gsTemp = mean(Mean.Temp,na.rm=T),gsPrec = sum(Precipitation,na.rm=T))
#lag (previous year) effects
gsTm1     <- mutate(gs, year=year+1)     ; gsTm1 = rename(gsTm1,gsTempTm1=gsTemp,gsPrecTm1=gsPrec)
yrClimTm1 <- mutate(yrClim, year=year+1) ; yrClimTm1 = rename(yrClimTm1,yrTempTm1=yrTemp,yrPrecTm1=yrPrec)
climSummaries=list(yrClim,gs,gsTm1,yrClimTm1)

#Meteorological data
meteoData1=Reduce(function(...) merge(...),climSummaries)
meteoData2=meteoData1[,-1]^2
colnames(meteoData2)=paste0(names(meteoData1)[-1],"_2")
meteoData=cbind(meteoData1,meteoData2)


#Demographic quantities----------------------------------------------

#beetle
beetleGr=melt(beetle1,id.var="Year")
beetleGr$Site=1 #there is only one site!
names(beetleGr)=c("year","species","abund","site")
beetleGr=formatDemographicSite(beetleGr[,c("year","site","species","abund")])

#defoliate
defoliateGr=melt(defoliate1,id.var="Year")
defoliateGr$Site=1 #there is only one site!
names(defoliateGr)=c("year","species","abund","site")
defoliateGr=formatDemographicSite(defoliateGr[,c("year","site","species","abund")])

insectGr=rbind(beetleGr,defoliateGr)

#Hares
hare=cbind(formatMeteoDate(as.character(hare$date),c("year","month","day"),"-"),hare)
hare$l_ear[hare$l_ear=="-"]=NA #just 5 missing values
hare$r_ear[hare$r_ear=="-"]=NA #many more missing values
#fix slight mistake in grid names
hare$grid[hare$grid=="bonrip "]="bonrip"
hare$grid[hare$grid=="bonmat "]="bonmat"
#individuals observed per year and grid (site)
indivXYr=distinct(select(hare,year,grid,l_ear))
hareS=summarise(group_by(indivXYr,year,grid),count=n())
hareS$species="Snow_hare"
names(hareS)=c("year","site","abund","species")
hareGr=hareS[,c("year","site","species","abund")]
hareGr=formatDemographicSite(hareGr)


###MODELS#########################################################################################

#Insects------------------------------------------------------------------------------------------
#No analyses possible for insects, no multiple sites
#insectD=merge(insectGr,meteoData)

#tiff("Results1/BonanzaCreek_insect.tiff",unit="in",height=8,width=6.3,res=500,compression="lzw")

#par(mfrow=c(3,3),mar=c(2.5,2.5,1,0.1),mgp=c(1,0.5,0))
#resInsect=analysisRandom(pop=insectD,speciesL=unique(insectD$species),
#                      meteoVars=names(meteoData)[2:9],measure="count",
#                      organism="arthropod",LTERsite="bonanzaCreek")
#dev.off()


#Hares------------------------------------------------------------------------------------------
hareD=merge(hareGr,meteoData)

tiff("Results1/BonanzaCreek_hare.tiff",unit="in",height=3.15,width=3.15,res=500,compression="lzw")

par(mfrow=c(1,1),mar=c(2.5,2.5,1,0.1),mgp=c(1.3,0.5,0))
resHare=analysisRandom(pop=hareD,speciesL="Snow_hare",
                       meteoVars=names(meteoData)[2:9],measure="count",
                       organism="mammal",LTERsite="bonanzaCreek")
dev.off()

#Writeup results
results=rbind(resHare)
write.csv(results,"Results1/resultsBonanzaCreek.csv",row.names=F)



#OLD CODE############################################################################

#Beetle--------------------------------------------------------------------------------------
#beetleTmp=beetle1
#defoliateTmp=defoliate1
#for(i in 2:ncol(beetleTmp)){ beetleTmp[,i]=logGr(beetleTmp[,i]) }
#for(i in 2:ncol(defoliateTmp)){ defoliateTmp[,i]=logGr(defoliateTmp[,i]) }

#beetleGr=melt(beetleTmp,id.vars = "Year")
#defoliateGr=melt(defoliateTmp,id.vars = "Year")
#insectGr=rbind(beetleGr,defoliateGr)
#names(insectGr)=c("year","species","logGr")
#insectGr=insectGr[-which(is.na(insectGr$logGr)),]

#Find species with at least 10 years of data
#tmp=summarise(group_by(insectGr,species),count=n())
#insectList=as.character(tmp$species[tmp$count>10])

#Hares-----------------------------------------------------------------------------------------
#hare=cbind(formatMeteoDate(as.character(hare$date),c("year","month","day"),"-"),hare)
#hare$l_ear[hare$l_ear=="-"]=NA #just 5 missing values
#hare$r_ear[hare$r_ear=="-"]=NA #many more missing values
##fix slight mistake in grid names
#hare$grid[hare$grid=="bonrip "]="bonrip"
#hare$grid[hare$grid=="bonmat "]="bonmat"
#individuals observed per year and grid (site)
#indivXYr=distinct(select(hare,year,grid,l_ear))
#hareS=summarise(group_by(indivXYr,year,grid),count=n())
#hareWide=dcast(hareS,year ~ grid,value.var="count")
#for(i in 2:ncol(hareWide)){ hareWide[,i]=logGr(hareWide[,i]) }

#Format data
#hareGr=melt(hareWide,id.vars="year")
#hareGr=hareGr[-which(is.na(hareGr$value)),][,-2]
#names(hareGr)=c("year","logGr")
#hareGr$species="hare" ; hareGr=hareGr[,c("year","species","logGr")]
