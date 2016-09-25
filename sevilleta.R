#THE PLAN:
#Climate data: http://sev.lternet.edu/data/sev-1
#Demographic data: http://sev.lternet.edu/data/sev-106
setwd("C:/Users/ac79/MEGA/Projects/RICE/LTER/")
source("Analysis1/functions.R")

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

#Find species with 1. more than 9 years of data, and 2. continuous presence across the series
spp=names(table(hoppers$SPECIES)[order(table(hoppers$SPECIES),decreasing=T)])[1:40]
#eschew species with less than 21 years
noSpp=NULL
zeroI=NULL
for(i in 1:length(spp)){
  oneSpp=subset(hoppers,SPECIES==spp[i])
  oneSpp=aggregate(CNT ~ YEAR + SITE,sum,data=oneSpp)
  siteD=dcast(oneSpp,YEAR ~ SITE,sum,value.var="CNT")
  if(contiguousY(siteD$YEAR)==F | length(siteD$YEAR)<11) noSpp=c(noSpp,i)
}
sppList=spp[setdiff(1:40,noSpp)]
sppList=sppList[-10] #Number 10 is NONE, which is NOT a species!

#Given the file 'siteVsGap', I suggest to aggregate across sites,
#except for species 1 (TRPA)
hoppersS=summarise(group_by(hoppers,YEAR,SITE,SPECIES),count=sum(CNT,na.rm=T))
hoppersGr=list()
for(i in 1:length(sppList)){
 
  oneSpp=subset(hoppersS,SPECIES==sppList[i])
  siteD=dcast(oneSpp,YEAR ~ SITE,sum,value.var="count")
  for(s in 2:ncol(siteD)) { siteD[,s]=logGr(siteD[,s])}
  sppLong=melt(siteD,id.vars="YEAR")
  sppLong$species=sppList[i]
  hoppersGr[[i]]=sppLong
}
hoppersGr=Reduce(function(...) rbind(...),hoppersGr) 
hoppersGr=hoppersGr[-which(is.na(hoppersGr$value)),]
names(hoppersGr)=c("year","site","logGr","species")
hoppersGr=hoppersGr[,c("year","site","species","logGr")]


#Coyote------------------------------------------------------------
coyote$species="coyote"
coyoteS=summarise(group_by(coyote,year,species),density=mean(density,na.rm=T))
coyoteS$density=logGr(coyoteS$density)
names(coyoteS)[3]="logGr"
coyoteS=coyoteS[-which(is.na(coyoteS$logGr)),]

#Rabbits------------------------------------------------------------
rabbit$species="rabbit"
rabbitS=summarise(group_by(rabbit,year,species),density=mean(density_estimate,na.rm=T))
rabbitS$density=logGr(rabbitS$density)
names(rabbitS)[3]="logGr"
rabbitS=rabbitS[-which(is.na(rabbitS$logGr)),]

#Rabbits2--------------------------------------------------------------------------------
rabbits2S=summarise(group_by(rabbits,Year,month,leg,species),count=n())
rabbits2S$species[rabbits2S$species==""]=NA
rabbSpp=unique(rabbits2S$species)[-3]
rabbGr=list()
for(i in 1:length(rabbSpp)){
  
  tmp=dcast(subset(rabbits2S,species==rabbSpp[i]),Year~ leg,sum,value.var="count")
  for(l in 2:ncol(tmp)){ tmp[,l]=logGr(tmp[,l]) }
  tmp$species=rabbSpp[i]
  rabbGr[[i]]=tmp
  rabbGr[[i]]=melt(rabbGr[[i]],id.vars=c("Year","species"))
}
rabbGr=Reduce(function(...) rbind(...),rabbGr)
rabbGr=rabbGr[-which(is.na(rabbGr$value)),]
rabbGr=select(rabbGr,Year,variable,species,value)
names(rabbGr)=c("year","site","species","logGr")


#Ground arthropods----------------------------------------------------------------------------
arthrop$Count[arthrop$Count==-888]=NA
arthrop=arthrop[-which(arthrop$Genus==-888),] #Exclude NAs on 

##Spp replication withint Site
arthropS=summarise(group_by(arthrop,Site,Year,Genus,Species),count=sum(Count,na.rm=T))
arthropS$species=paste(arthropS$Genus,arthropS$Species,sep="_")
tmp=as.data.frame(arthropS)[,c("Year","species")]
arthList=colnames(table(tmp))[apply(table(tmp),2,sum)>20] 

#Population growth data
arthGr=list()
for(i in 1:length(arthList)){
  
  sppData=subset(arthropS,species==arthList[i])
  sppWide=dcast(sppData,Year ~ Site, sum,value.var="count")
  if(contiguousY(sppWide$Year)==F) { 
    time=data.frame(Year=seq(min(sppWide$Year),max(sppWide$Year),1))
    sppWide=merge(time,sppWide,all.x=T)
  }
  for(k in 2:ncol(sppWide)) { sppWide[,k]=logGr(sppWide[,k]) }
  sppWide$species=arthList[i]
  arthGr[[i]]=melt(sppWide,id.var=c("Year","species"))
  names(arthGr[[i]])[3:4]=c("site","logGr")
  
}
arthGr=Reduce(function(...) rbind(...),arthGr)
arthGr=arthGr[-which(is.na(arthGr$logGr)),]
arthGr=arthGr[,c("Year","site","species","logGr")]
names(arthGr)[1]="year"

#Small mammals-------------------------------------------------------------------------
mammS=summarise(group_by(mamm,year,location,species),count=n())
mammList=colnames(table(mammS[,c("year","species")]))
mammList=mammList[which(apply(table(mammS[,c("year","species")])>0,2,sum)>18)]

mammGr=list()
#Population growth data
for(i in 1:length(mammList)){
  
  tmp=subset(mammS,species==mammList[i])
  mammWide=dcast(tmp,year ~ location,sum,value.var="count")
  if(contiguousY(mammWide$year)==F) { 
    time=data.frame(year=seq(min(mammWide$year),max(mammWide$year),1))
    mammWide=merge(time,mammWide,all.x=T)
  }
  for(k in 2:ncol(mammWide)) { mammWide[,k]=logGr(mammWide[,k]) }
  mammWide$species=mammList[i]
  mammGr[[i]]=melt(mammWide,id.var=c("year","species"))
  names(mammGr[[i]])[3:4]=c("site","logGr")
  
}
mammGr=Reduce(function(...) rbind(...),mammGr)
mammGr=mammGr[-which(is.na(mammGr$logGr)),]



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

par(mfrow=c(3,3),mar=c(2.5,2.4,1,0.1),mgp=c(1.3,0.5,0),lwd=1)
resHoppers=analysisRandom(pop=hoppersD,speciesL=sppList,
                     meteoVars=names(pjMeteoData)[2:6],measure="count",
                     organism="arthropod",LTERsite="Sevilleta")
dev.off()


#coyotes--------------------------------------------
coyoteD=merge(coyoteS,pjMeteoData)

tiff("Results1/sevilleta_Mammals.tiff",unit="in",height=8,width=6.3,res=500,compression="lzw")
par(mfrow=c(4,4),mar=c(2.5,2.4,1,0.1),mgp=c(1.3,0.5,0),lwd=1)

#Small rodents---------------------------------------------------------------
mammD=merge(mammGr,pjMeteoData)

resMamm=analysisRandom(pop=mammD,speciesL=mammList,
                     meteoVars=names(pjMeteoData)[2:7],measure="count",
                     organism="mammal",LTERsite="Sevilleta")

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
rabbit2D=merge(rabbGr,pjMeteoData)

resRabb2=analysisRandom(pop=rabbit2D,speciesL=unique(rabbGr$species),
                     meteoVars=names(pjMeteoData)[2:7],measure="count",
                     organism="mammal",LTERsite="Sevilleta")
dev.off()


#Ground arthropods----------------------------------------------------------------
arthD=merge(arthGr,pjMeteoData)

tiff("Results1/sevilleta_Arthr.tiff",unit="in",height=8,width=6.3,res=500,compression="lzw")
par(mfrow=c(5,4),mar=c(2.5,2.4,1,0.1),mgp=c(1.3,0.5,0),lwd=1)
resArth=analysisRandom(pop=arthD,speciesL=arthList,
                      meteoVars=names(pjMeteoData)[2:7],measure="count",
                      organism="arthropod",LTERsite="Sevilleta")
dev.off()


#Store summary results
results=rbind(resHoppers,resRabb2,resArth,resMamm)
write.csv(results,"Results1/resultsSevilleta.csv",row.names=F)
