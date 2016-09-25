#This script runs preliminary analyses on LTER data
setwd("C:/Users/ac79/MEGA/Projects/RICE/LTER/")
source("Analysis1/functions.R")

#DATA
#Climate data
troutClimate=read.csv("Data/NorthTemperateLakes/TroutLakeClimate.csv")
mendotaClimate=read.csv("Data/NorthTemperateLakes/MendotaLakeClimate.csv")
slClima=read.csv("Data/NorthTemperateLakes/NTL_daily_meteorological_sparkling_lake_raft.csv")
#airport=read.csv("Data/NorthTemperateLakes/NTL_daily_meteorological_woodruff_airport.csv")
minocqua=read.csv("Data/NorthTemperateLakes/minocqua_UCC_ghcn_USC00475516_2016.csv",skip=15)


#Demography TO review, go to: https://lter.limnology.wisc.edu/datacatalog/search
#Then choose "NTL Core Datasets" for project and "population" as theme
#Fishes
fish=read.csv("Data/NorthTemperateLakes/NTL_fish_abundance_StanleyEmily.csv")
#Benthic macroinvertebrates
ben=read.csv("Data/NorthTemperateLakes/NTL_benthic_macroinvertebrates_StanleyEmily.csv")
#Crayfish
cray=read.csv("Data/NorthTemperateLakes/NTL_crayfish_abundance_StanleyEmily.csv")
#Pelagic
pela=read.csv("Data/NorthTemperateLakes/NTL_pelagic_macroinvertabrate_abundance.csv")
#Plankton
pla1=read.csv("Data/NorthTemperateLakes/NTL_phytoplankton_madisonLakes.csv")
#Only use Maidon lakes, the second does not have 10 years of data
#pla2=read.csv("Data/NorthTemperateLakes/NTL_phytoplankton_troutLake.csv")
#Zooplankton
zoo1=read.csv("Data/NorthTemperateLakes/NTL_zooplankton_madisonLakes.csv")
zoo2=read.csv("Data/NorthTemperateLakes/NTL_zooplankton_troutLake.csv")

##################################################################################################
#FORMAT
##################################################################################################

#DEMOGRAPHIC DATA------------------------------------------------------------------------------

#Fishes
fish$density=fish$total_caught / fish$effort
names(fish)[c(1,2,5,7)]=c("site","year","species","abund")
fishGr=formatDemographicSite(fish[,c("year","site","species","abund")])

#Benthic macroinvertebrates
names(ben)[c(1,2,3,5,7)]=c("site","year","sp_1","species","abund")
benGr=formatDemographicSite(ben[,c("year","site","species","abund")])

#Crayfish
cray$density=cray$total_caught / cray$effort
names(cray)[c(1,2,5,7)]=c("site","year","species","abund")
crayGr=formatDemographicSite(cray[,c("year","site","species","abund")])

#Pelagic macroinvertebrates
names(pela)[c(1,2,5,7)]=c("site","year","species","abund")
pelaGr=formatDemographicSite(pela[,c("year","site","species","abund")])

#Phytoplankton
names(pla1)[c(1,2,10,15)]=c("site","year","abund","species")
pla1=pla1[-which(pla1$species==""),]
plaGr=formatDemographicSite(pla1[,c("year","site","species","abund")])

#Zooplankton
names(zoo1)[c(1,2,7,8)]=c("site","year","species","abund")
names(zoo2)[c(1,2,6,7)]=c("site","year","species","abund")
zoo=rbind(zoo2[,c("site","year","species","abund")],zoo1[,c("site","year","species","abund")])
zooGr=formatDemographicSite(zoo)


#Climate data---------------------------------------------------------------------------------- 
#format date
minocqua=cbind(formatMeteoDate(minocqua$Day,c("year","month","day"),"-"),minocqua)
#introduce NAs
minocqua$Precipitation[minocqua$Precipitation=="T"]=NA
minocqua$Precipitation[minocqua$Precipitation=="M"]=NA
minocqua$Max.Temperature[minocqua$Max.Temperature=="M"]=NA
minocqua$Min.Temperature[minocqua$Min.Temperature=="M"]=NA
#"as numeric"
for(i in 5:7){ minocqua[,i]=as.numeric(minocqua[,i])}
minocqua$AvgT=(minocqua$Min.Temperature + minocqua$Max.Temperature)/2

#Summaries
yRClim=yrClimTm1=summarise(group_by(minocqua,year),yrTemp=mean(AvgT,na.rm=T),
                           yrPrec=sum(Precipitation,na.rm=T))
yrClimTm1$year=yrClimTm1$year+1
names(yrClimTm1)=c("year","yrTempTm1","yrPrecTm1")
meteoData1=merge(yRClim,yrClimTm1)
meteoData2=meteoData1[,-1]^2
colnames(meteoData2)=paste0(names(meteoData1[,-1]),"_2")
meteoData=cbind(meteoData1,meteoData2)


########################################################################################################################################################
#ANALYSIS############################################################################
########################################################################################################################################################

#Fishes----------------------------------------------------------------------------------------
fishD=merge(fishGr,meteoData)

tiff("Results1/NTL_Fishes.tiff",unit="in",width=6.3,height=8,res=500,compression="lzw")
par(mfrow=c(7,5),mar=c(2,2,1,0.2),mgp=c(1,0.5,0))
resFish=analysisRandom(pop=fishD,
                     meteoVars=names(meteoData)[2:5],measure="density",
                     organism="fish",LTERsite="NorthTemperateLakes")
dev.off()

#Benthic Macroinvertebrates---------------------------------------------------------------------------------
benD=merge(benGr,meteoData)

tiff("Results1/NTL_BenthicInvertebrates.tiff",unit="in",width=6.3,height=8,res=500,compression="lzw")
par(mfrow=c(5,2),mar=c(2,2,1,0.2),mgp=c(1,0.5,0))
resBen=analysisRandom(pop=benD,
                     meteoVars=names(meteoData)[2:5],measure="count",
                     organism="macroinv",LTERsite="NorthTemperateLakes")
dev.off()

#Crayfishes--------------------------------------------------------------------------
crayD=merge(crayGr,meteoData)

tiff("Results1/NTL_Crayfish.tiff",unit="in",width=6.3,height=6.3,res=500,compression="lzw")
par(mfrow=c(1,2),mar=c(2,2,1,0.2),mgp=c(1,0.5,0))
resCray=analysisRandom(pop=crayD,
                    meteoVars=names(meteoData)[2:5],measure="density",
                    organism="crayfish",LTERsite="NorthTemperateLakes")
dev.off()

#Pelagic macroinvertebrates------------------------------------------------------------------------------
pelaD=merge(pelaGr,meteoData)

tiff("Results1/NTL_Pelagic.tiff",unit="in",width=6.3,height=9,res=500,compression="lzw")
par(mfrow=c(2,2),mar=c(2,2,1,0.2),mgp=c(1,0.5,0))
resPela=analysisRandom(pop=pelaD,
                     meteoVars=names(meteoData)[2:5],measure="count",
                     organism="macroinv",LTERsite="NorthTemperateLakes")
dev.off()

#Phytoplankton------------------------------------------------------------------------------
plaD=merge(plaGr,meteoData)

tiff("Results1/NTL_Plankton.tiff",unit="in",width=6.3,height=9,res=500,compression="lzw")
par(mfrow=c(6,5),mar=c(2,2,1,0.2),mgp=c(1,0.5,0))
resPla=analysisRandom(pop=plaD,
                     meteoVars=names(meteoData)[2:5],measure="density",
                     organism="plankton",LTERsite="NorthTemperateLakes")
dev.off()

#Zooplankton
zooD=merge(zooGr,meteoData)

tiff("Results1/NTL_Zooplankton1.tiff",unit="in",width=6.3,height=9,res=500,compression="lzw")
par(mfrow=c(7,5),mar=c(2,2,1,0.2),mgp=c(1,0.5,0))
resZoo=analysisRandom(pop=zooD,speciesL=unique(zooD$species)[1:35],
                      meteoVars=names(meteoData)[2:5],measure="density",
                      organism="zooplankton",LTERsite="NorthTemperateLakes")
dev.off()
tiff("Results1/NTL_Zooplankton2.tiff",unit="in",width=6.3,height=9,res=500,compression="lzw")
par(mfrow=c(6,5),mar=c(2,2,1,0.2),mgp=c(1,0.5,0))
resZoo=analysisRandom(pop=zooD,speciesL=unique(zooD$species)[36:70],
                      meteoVars=names(meteoData)[2:5],measure="density",
                      organism="zooplankton",LTERsite="NorthTemperateLakes")
dev.off()


#write out
results=rbind(resFish,resBen,resCray,resPela,resPla,resZoo)
write.csv(results,"Results1/resultsNorthTempLakes.csv",row.names=F)
  