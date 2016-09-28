#THE PLAN:
#Climate data: http://sev.lternet.edu/data/sev-1
#Demographic data: http://sev.lternet.edu/data/sev-106
setwd("C:/Users/ac79/MEGA/Projects/LTER/")
source("C:/Users/ac79/Documents/CODE/climate_vs_r/functions.R")

#Data
#Meteo
climateSEV=read.csv("Data/Sevilleta/sev1_meteorology_20140111.csv")

#Demographic Data
#Grasshopppers
hoppers=read.table("Data/Sevilleta/sev106_hopperdynamics_20150826.txt",header=T,sep=",")
#Coyote
coyote=read.table("Data/Sevilleta/sev112_coyotedens_01062009_RobertParmenter.txt",header=T,sep=",")
#Rabbits
rabbit=read.table("Data/Sevilleta/sev113_rabbitdens_20040226_RobertParmenter.txt",header=T,sep=",")
rabbits=read.table("Data/Sevilleta/sev023_rabbitpopns_20150310_RobertParmenter.txt",header=T,sep=",")
#Ground arthopod
arthrop=read.csv("Data/Sevilleta/sev029_arthropop_02162009_0_DavidLightfoot.csv",header=T,sep=",")
#Small mammals
mamm=read.csv("Data/Sevilleta/sev008_rodentpopns_20160701_SethNewsome.csv")
#Ants


###################################################################################################
#Format#########################################
###################################################################################################

#Species count data---------------------------------------------------------------------------------
hoppers$YEAR=as.numeric(substr(hoppers$PER,1,4))
hoppers$SEASON=substr(hoppers$PER,5,5)
hoppers=rename(hoppers,year=YEAR,site=SITE,species=SPECIES,abund=CNT)
hoppersGr=formatDemographicSite(hoppers)

#Coyote------------------------------------------------------------
#coyote$species="coyote"
#coyoteS=summarise(group_by(coyote,year,species),density=mean(density,na.rm=T))
#coyoteS$density=logGr(coyoteS$density)
#names(coyoteS)[3]="logGr"
#coyoteS=coyoteS[-which(is.na(coyoteS$logGr)),]

#Rabbits------------------------------------------------------------
#rabbit$species="rabbit"
#rabbitS=summarise(group_by(rabbit,year,species),density=mean(density_estimate,na.rm=T))
#rabbitS$density=logGr(rabbitS$density)
#names(rabbitS)[3]="logGr"
#rabbitS=rabbitS[-which(is.na(rabbitS$logGr)),]

#Rabbits2--------------------------------------------------------------------------------
rabbits <- rename(rabbits, year = Year, site = leg)
rabbits_cnt <- summarise(group_by(rabbits,year,site,species),abund=n())
rabbits_cnt <- subset(rabbits_cnt , species != "")
rabbitsGr   <- formatDemographicSite(rabbits_cnt)

# Ground arthropods----------------------------------------------------------------------------
arthrop$Count[arthrop$Count==-888]=NA
# Only species at the Genus level
arthrop=arthrop[-which(arthrop$Genus==-888),]

# Spp replication withint Site
arthropS=summarise(group_by(arthrop,Site,Year,Genus,Species),count=sum(Count,na.rm=T))
arthropS$species=paste(arthropS$Genus,arthropS$Species,sep="_")
arthropS <- rename(arthropS, year = Year, site = Site, abund = count)
arthGr   <- formatDemographicSite(arthropS)

# Small mammals-------------------------------------------------------------------------
mammS  <- summarise(group_by(mamm,year,location,species),count=n())
mammS  <- rename(mammS, site = location, abund = count)
mammGr <- formatDemographicSite(mammS) 


#Climate data---------------------------------------------------------------------------------
climateSEV$Temp_C[climateSEV$Temp_C==-999]=NA
climateSEV$Temp_C[climateSEV$Temp_C < -39.9]=NA
climateSEV$Precip[climateSEV$Precip==-999]=NA
climateSEV$Precip[climateSEV$Precip<0]=NA
climateSEV=subset(climateSEV,Station_ID == 42 | Station_ID == 49 | Station_ID == 50)#exclude sites we don't care about!

#Daily means/sums
dailyClim=summarise(group_by(climateSEV,Year,Jul_Day,Station_ID),
                    max_temp=max(Temp_C,na.rm=T),min_temp=min(Temp_C,na.rm=T),
                    precip=sum(Precip,na.rm=T))
dailyClim$avg_Temp=(dailyClim$max_temp + dailyClim$min_temp)/2 

#Spring and summer temp./precip.
spr=summarise(group_by(subset(dailyClim, Jul_Day<152),Year,Station_ID),
              sprTemp=mean(avg_Temp,na.rm=T),sprPrec=sum(precip,na.rm=T))
summ=summarise(group_by(subset(dailyClim, Jul_Day<274),Year,Station_ID),
               summTemp=mean(avg_Temp,na.rm=T),summPrec=sum(precip,na.rm=T))
yrClim=summarise(group_by(dailyClim,Year,Station_ID),yrTemp=mean(avg_Temp,na.rm=T),yrPrec=sum(precip,na.rm=T))
yrClim$Year=yrClim$Year+1      #Associate this data with NEXT year!

climateList=list(spr,summ,yrClim)
annualClimate=Reduce(function(...) merge(...),climateList)
#I use the PINION-JUNIPER data (Cerro Montoso)
pjClimate1=subset(annualClimate, Station_ID==42)[,-2]
names(pjClimate1)[1]="year"

#Create squared values
pjClimate2=pjClimate1[,-1]^2
colnames(pjClimate2)=paste0(names(pjClimate1)[-1],"_2")
pjMeteoData=cbind(pjClimate1,pjClimate2)


########################################################################################################################################################
#ANALYSIS############################################################################
########################################################################################################################################################
hoppersD=merge(hoppersGr,pjMeteoData)

#Summary output information (this will eventually include sensitivity)
tiff("Results1/sevilleta_Grasshoppers.tiff",unit="in",height=6.3,width=6.3,res=500,compression="lzw")

par(mfrow=c(4,4),mar=c(2.5,2.4,1,0.1),mgp=c(1.3,0.5,0),lwd=1)
resHoppers=analysisRandom(pop=hoppersD,
                     meteoVars=names(pjMeteoData)[2:6],measure="count",
                     organism="arthropod",LTERsite="Sevilleta")
dev.off()


#coyotes--------------------------------------------
#coyoteD=merge(coyoteS,pjMeteoData)

#tiff("Results1/sevilleta_Mammals.tiff",unit="in",height=8,width=6.3,res=500,compression="lzw")
#par(mfrow=c(4,4),mar=c(2.5,2.4,1,0.1),mgp=c(1.3,0.5,0),lwd=1)

#Small rodents---------------------------------------------------------------
mammD=merge(mammGr,pjMeteoData)

tiff("Results1/sevilleta_Mammals.tiff",unit="in",height=8,width=6.3,res=500,compression="lzw")

par(mfrow=c(5,3),mar=c(2.5,2.4,1,0.1),mgp=c(1.3,0.5,0),lwd=1)
resMamm=analysisRandom(pop=mammD,
                     meteoVars=names(pjMeteoData)[2:7],measure="count",
                     organism="mammal",LTERsite="Sevilleta")
dev.off()

#Coyote-----------------------------------------------------------------------
#resCoyote=analysisRandom(pop=coyoteD,speciesL="coyote",
#                        meteoVars=names(pjMeteoData)[2:7],measure="density",
#                        organism="mammal",LTERsite="Sevilleta")

#Rabbits----------------------------------------------------------------
#rabbitD=merge(rabbitS,pjMeteoData)

#resRabb=analysisRandom(pop=rabbitD,speciesL="rabbit",
#                       meteoVars=names(pjMeteoData)[2:7],measure="density",
#                       organism="mammal",LTERsite="Sevilleta")

#Rabbits2----------------------------------------------------------------
rabbit2D=merge(rabbitsGr,pjMeteoData)

par(mfrow=c(1,2))
resRabb2=analysisRandom(pop=rabbit2D,
                     meteoVars=names(pjMeteoData)[2:7],measure="count",
                     organism="mammal",LTERsite="Sevilleta")


#Ground arthropods----------------------------------------------------------------
arthD=merge(arthGr,pjMeteoData)

tiff("Results1/sevilleta_Arthr.tiff",unit="in",height=8,width=6.3,res=500,compression="lzw")
par(mfrow=c(1,2),mar=c(2.5,2.4,1,0.1),mgp=c(1.3,0.5,0),lwd=1)
resArth=analysisRandom(pop=arthD,#speciesL=arthList,
                      meteoVars=names(pjMeteoData)[2:7],measure="count",
                      organism="arthropod",LTERsite="Sevilleta")
dev.off()


#Store summary results
results=rbind(resHoppers,resRabb2,resArth,resMamm)
write.csv(results,"Results1/resultsSevilleta.csv",row.names=F)
