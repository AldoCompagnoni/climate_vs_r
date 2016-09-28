#Data on Palmer's zooplankton data
setwd("C:/Users/ac79/MEGA/Projects/LTER/")
#Functions (AICTable, formatMeteoDate, logGr, contiguousY)
source("C:/Users/ac79/Documents/CODE/climate_vs_r/functions.R")

#DATA#####################################################################

#Meteorological data------------------------------------------------------
meteo=read.csv("Data/Palmer/PalmerStationWeather_Daily.csv")

#Zooplankton data---------------------------------------------------------
zoo1=read.csv("Data/Palmer/Zooplankton Density_Historical.csv")
zoo2=read.csv("Data/Palmer/Zooplankton Density_Current.csv")

#Penguin data--------------------------------------------------------------
peng=read.csv("Data/Palmer/AdeliePenguinCensus.csv")

#Bird Census data----------------------------------------------------------
birdStat=read.csv("Data/Palmer/Bird Census Stationary - Summer.csv")


#FORMAT####################################################################

#Zooplankton----------------------------------------------------------------------------------------
zoo1$year=formatMeteoDate(zoo1$Date.GMT,c("year","month","day"),"-")[,1]
zoo1tmp=zoo1[,c(6,27:ncol(zoo1))]
zooLong=melt(zoo1tmp,id.vars=c("year","GridLine"))
names(zooLong)=c("year","site","species","abund")

zoo2$year=formatMeteoDate(zoo2$Date.GMT,c("year","month","day"),"-")[,1]
zoo2tmp=zoo2[,c(6,28:ncol(zoo2))]
zooLong2=melt(zoo2tmp,id.vars=c("year","GridLine"))
names(zooLong2)=c("year","site","species","abund")

zooAll=rbind(zooLong,zooLong2)
zooAll$species=as.character(zooAll$species)
zooGr=formatDemographicSite(zooAll)


#Penguin data----------------------------------------------------------------------------------------
peng=cbind(formatMeteoDate(peng$Date.GMT,format=c("year","month","day"),"-"),peng)
pengTmp=peng[,c("year","Island","Breeding.Pairs")]
pengTmp$species="AdeliePenguin"
names(pengTmp)=c("year","site","abund","species") 
pengTmp=pengTmp[,c("year","site","species","abund")]
pengGr=formatDemographicSite(pengTmp)

#Birds--------------------------------------------------------------------------

birdStat$year=NA
#identify years sing substr()
fourID=which(nchar(birdStat$Year...Month)==4)
threeID=which(nchar(birdStat$Year...Month)==3)
oneID=which(nchar(birdStat$Year...Month)==1)
birdStat$year[fourID]=substr(birdStat$Year...Month,1,2)[fourID]
birdStat$year[threeID]=paste0("200",substr(birdStat$Year...Month,1,1)[threeID])
birdStat$year[oneID]="2000"
#Finish creating "year" column
birdStat$year=as.numeric(birdStat$year)
nietiesID=which(birdStat$year<100 & birdStat$year>50)  
twoKID=which(birdStat$year>0 & birdStat$year<20)  
birdStat$year=as.character(birdStat$year)
birdStat$year[nietiesID]=paste0("19",birdStat$year[nietiesID])
birdStat$year[twoKID]=paste0("20",birdStat$year[twoKID])

birdTmp=birdStat[,c("year","Station","Species","Number")]
names(birdTmp)=c("year","site","species","abund")
birdTmp$year=as.numeric(birdTmp$year)
birdGr=formatDemographicSite(birdTmp)


#Meteorological data------------------------------------------------------
meteo=cbind(formatMeteoDate(meteo$Date,format=c("year","month","day"),"-"),meteo)
meteo=meteo[,c("year","month","day","Sea.Surface.Temperature..C.","Temperature..avg...C.",
         "Precipitation..melted...mm.")]
names(meteo)=c("year","month","day","SeaSurfaceTempC","TempAvgC","Precip_mm")

#There are no NA (-999) or TRACE (-99) values
sum(meteo$SeaSurfaceTempC==-9999,na.rm=T)
sum(meteo$SeaSurfaceTempC==-999,na.rm=T)
sum(meteo$SeaSurfaceTempC==-99,na.rm=T)
sum(meteo$TempAvgC==-9999,na.rm=T)
sum(meteo$TempAvgC==-999,na.rm=T)
sum(meteo$TempAvgC==-999,na.rm=T)
sum(meteo$TempAvgC==-99,na.rm=T)
sum(meteo$Precip_mm==-998,na.rm=T)
sum(meteo$Precip_mm==-99,na.rm=T)

##Averages (year +1)
tmpMet=subset(subset(meteo,year>1992))
tmpMet$year=tmpMet$year+1 #This associates temp of year t with growth r in year t+1  
yrSea=as.data.frame(summarise(group_by(tmpMet, year),yrSea=mean(SeaSurfaceTempC,na.rm=T)))
yrTemp=as.data.frame(summarise(group_by(tmpMet, year),yrTemp=mean(TempAvgC,na.rm=T)))
#I call "growing season" December of previous year, plus january of the next
#Clunky code because "year" does not align with growing season
tmpMet=subset(subset(meteo,year>1990), month==12)
tmpMet$year=tmpMet$year+1 #Associate this with "next year"
gsSea1=as.data.frame(summarise(group_by(tmpMet, year),gsSea1=mean(SeaSurfaceTempC,na.rm=T)))
gsTemp1=as.data.frame(summarise(group_by(tmpMet, year),gsTemp1=mean(TempAvgC,na.rm=T)))
tmpMet=subset(subset(meteo,year>1990), month==1)
gsSea2=as.data.frame(summarise(group_by(tmpMet, year),gsSea2=mean(SeaSurfaceTempC,na.rm=T)))
gsTemp2=as.data.frame(summarise(group_by(tmpMet, year),gsTemp2=mean(TempAvgC,na.rm=T)))
gsSea=merge(gsSea1,gsSea2)    ; gsSea$gsSea=apply(gsSea[,2:3],1,mean,na.rm=T)
gsTemp=merge(gsTemp1,gsTemp2) ; gsTemp$gsTemp=apply(gsTemp[,2:3],1,mean,na.rm=T)
gsSea=gsSea[,c("year","gsSea")]
gsTemp=gsTemp[,c("year","gsTemp")]

#Put it all together
meteoList=list(yrSea,yrTemp,gsSea,gsTemp)
metTemp=Reduce(function(...) merge(...),meteoList)
#Cubit terms
meteoData1=cbind(metTemp,metTemp[,-1]^2)
names(meteoData1)[6:9]=paste0(names(metTemp)[2:5],"_2")
meteoData=meteoData1 #This is the "standard" for the zooplankton data


##Averages, year = year (dataset only for Penguins)
tmpMet=subset(subset(meteo,year>1992))
yrSea=as.data.frame(summarise(group_by(tmpMet, year),yrSeaTm1=mean(SeaSurfaceTempC,na.rm=T)))
yrTemp=as.data.frame(summarise(group_by(tmpMet, year),yrTempTm1=mean(TempAvgC,na.rm=T)))
#I call "growing season" December of previous year, plus january of the next
#Clunky code because "year" does not align with growing season
tmpMet=subset(subset(meteo,year>1990), month==12)
gsSea1=as.data.frame(summarise(group_by(tmpMet, year),gsSea1=mean(SeaSurfaceTempC,na.rm=T)))
gsTemp1=as.data.frame(summarise(group_by(tmpMet, year),gsTemp1=mean(TempAvgC,na.rm=T)))
tmpMet=subset(subset(meteo,year>1990), month==1)
gsSea2=as.data.frame(summarise(group_by(tmpMet, year),gsSea2=mean(SeaSurfaceTempC,na.rm=T)))
gsTemp2=as.data.frame(summarise(group_by(tmpMet, year),gsTemp2=mean(TempAvgC,na.rm=T)))
gsSea=merge(gsSea1,gsSea2)    ; gsSea$gsSeaTm1=apply(gsSea[,2:3],1,mean,na.rm=T)
gsTemp=merge(gsTemp1,gsTemp2) ; gsTemp$gsTempTm1=apply(gsTemp[,2:3],1,mean,na.rm=T)
gsSea=gsSea[,c("year","gsSeaTm1")]
gsTemp=gsTemp[,c("year","gsTempTm1")]

#Put it all together
meteoList=list(yrSea,yrTemp,gsSea,gsTemp)
metTemp=Reduce(function(...) merge(...),meteoList)
#Cubit terms
meteoData2=cbind(metTemp,metTemp[,-1]^2)
names(meteoData2)[6:9]=paste0(names(metTemp)[2:5],"_2")



#Zooplankton data---------------------------------------------------------
#names(zoo1)[c(7,10,16)]=c("GridStation","TimeStart.CLST","TowDuration")
#zoo1$ZAID=NA

#zoo1=zoo1[,intersect(names(zoo1),names(zoo2))]
#zoo2=zoo2[,intersect(names(zoo1),names(zoo2))]
#zoo=rbind(zoo1,zoo2)
#zoo=cbind(formatMeteoDate(zoo$Date.GMT,format=c("year","month","day"),"-"),zoo)
#Remove "messy" symbols
#names(zoo)=gsub("..num.1000m³.","",names(zoo))
#names(zoo)=gsub("..ml.1000m3.","",names(zoo))


#Select only GridLines with multiple replication
#gridLines=c(300,400,500,600)
#zooA=subset(zoo,GridLine == 300 | GridLine == 400 | GridLine == 500 | GridLine == 600)


#Retain only estimates of "NUMBERS"
#zooA=zooA[,c(1:30,grep("Num",names(zooA)))]

#remove "trace" values (-1) with arbitrarily small value (0.0001)
#remove NAs, put a 0 instead.
#for(col in 30:54){
#  tmp=zooA[,col]
#  tmp[tmp==-1]=0.0001
#  tmp[is.na(tmp)]=0
#  zooA[,col]=tmp
#}

#SPATIALLY REPLICATED DATA
#eval(parse(n=1,text=paste0("zooDat=aggregate(cbind(",paste0(names(zooA[,c(30:54)]),collapse=","),
#                           ") ~ GridLine + year,FUN=sum,na.rm=T,data=zooA)")))
#only four species abundant enough to be present in every 
#apply(zooDat[,-c(1,2)]>1,2,sum) #ThysanoeNum has ALL the DATA NEEDED!

#NOT-SPATIALLY REPLICATED FILE
#eval(parse(n=1,text=paste0("zooDatAll=aggregate(cbind(",paste0(names(zooA[,c(30:54)]),collapse=","),
#                           ") ~ year,FUN=sum,na.rm=T,data=zooA)")))
#apply(zooDatAll[,-c(1,2)]>1,2,sum)
#24 years:LimacinaNum,PseudosaNum,EsuperbaNum,EcrystalNum,ThysanoeNum
#23 years:TomopterNum,ThemistoNum

###Format zooplankton data------------------------------------------------------------------------
#zooReplic=zooDat[,c("GridLine","year","ThysanoeNum")]
#Thysanoe=zooReplic[1:(nrow(zooReplic)-1),]
#Thysanoe=dcast(zooReplic,year ~ GridLine,value.var = "ThysanoeNum")
#Thysanoe=cbind(Thysanoe,matrix(NA,nrow(Thysanoe),4))
#calculate growth rates
#for(i in 2:5){
#  Thysanoe[,i+4]=logGr(Thysanoe[,i])
#}
#Thysanoe=Thysanoe[,c(1,6:9)]
#Thysanoe=melt(Thysanoe,id.var="year")
#Thysanoe=merge(Thysanoe,meteoData)
#names(Thysanoe)[2:3]=c("species","logGr")
#Thysanoe$species="ThysanoeNum"

#nonreplicated species
#zooSppList=c("LimacinaNum","PseudosaNum","EsuperbaNum","EcrystalNum")
#zooNonRep=zooDatAll[,c("year",zooSppList)]
#zooNonRep=cbind(zooNonRep,matrix(NA,nrow(zooNonRep),length(zooSppList)))
#for(i in 2:5){
#  zooNonRep[,i+4]=logGr(zooNonRep[,i])
#}
#zooNonRep=zooNonRep[-1,c(1,6:9)]
#names(zooNonRep)[2:5]=zooSppList
#zooNonRep=melt(zooNonRep,id.var="year")
#zooNonRep=merge(zooNonRep,meteoData)
#names(zooNonRep)[2:3]=c("species","logGr")
#Little test
#ncol(zooNonRep)==ncol(Thysanoe)
#length(intersect(names(Thysanoe),names(zooNonRep)))==ncol(zooNonRep)
#zooD=rbind(Thysanoe,zooNonRep) #Full zooplankton Data set


#Penguin data----------------------------------------------------------------------------------------
#peng=read.csv("Data/Palmer/AdeliePenguinCensus.csv")
#peng=cbind(formatMeteoDate(peng$Date.GMT,format=c("year","month","day"),"-"),peng)
#pengTmp=dcast(peng,year ~ Island,sum,value.var="Breeding.Pairs")
#remove LIT, went extinct
#pengData=cbind(pengTmp[,-5],matrix(NA,nrow(pengTmp),4))
#calculate growth rates
#for(i in 2:5){
#  pengData[,i+4]=logGr(pengData[,i])
#}
#pengData=pengData[-1,c(1,6:9)]
#metPeng=merge(meteoData1,meteoData2) #Meteorological data for Adelie penguins
#pengA=merge(melt(pengData,id.var="year")[,-2],metPeng)
#names(pengA)[2]="logGr"


#Bird census data-------------------------------------------------------------------

#Control how many species have a continuous record of 21 years
#sppNames=names(table(birdStat$Species)[which(table(birdStat$Species)>100)])
#keepNames=NULL
#for(i in 1:length(sppNames)){
#  tmp=aggregate(Number ~ year,sum,data=subset(birdStat,Species==sppNames[i]))
#  if(nrow(tmp)==21) { keepNames=c(keepNames,sppNames[i])}
#}
#sppList=keepNames


###############################################################################################
#ANALYZE##########################
###############################################################################################


#Analyses and graphs---------------------------------------------------------------------------
zooD=merge(zooGr,meteoData1)
for(i in 1:length(unique(zooD$species))){
  newName=unique(zooD$species)[i]
  newName=gsub("..num.1000m³.","",newName)
  newName=gsub("..ml.1000m3.","",newName)
  zooD$species[zooD$species %in% unique(zooD$species)[i]]=newName
}

#Do it by hand
orig=c("LimacinaNum",   "Amphipod_oNum", "TeleosteNum",   "PseudosaNum", 
"EsuperbaNum",   "ThysanoeVol",   "SalpAggVol",    "ThysanoeNum",  
"TeleosteVol",   "EcrystalNum",   "Chaetogn_oNum", "PseudosaVol",  
"PantarctVol",   "EsuperbaVol",   "SalpAggNum",    "Chaetogn_oVol",
"PantarctNum",   "Polychae_oVol", "Amphipod_oVol", "TomopterVol",  
"TomopterNum",   "SiphonVol",     "ThemistoNum",   "CtenophNum",   
"ThemistoVol")
newNames=rep(NA,25)
for(i in 1:25) newNames[i]=gsub("_oNum","",orig[i])
for(i in 1:25) newNames[i]=gsub("_oVol","",newNames[i])
for(i in 1:25) newNames[i]=gsub("Num","",newNames[i])
for(i in 1:25) newNames[i]=gsub("Vol","",newNames[i])
for(i in 1:25) zooD$species[zooD$species %in% orig[i]]=newNames[i]

zooD <- zooD[-which(zooD$Nt0 < 0),]
zooD <- zooD[-which(is.na(zooD$logNt0)),]

tiff("Results1/plamer_zooplankton.tiff",unit="in",height=5,width=6.3,res=500,compression="lzw")

par(mfrow=c(6,5),mar=c(2.3,2.8,1,0.1),mgp=c(1.3,0.5,0))
resZoo=analysisRandom(pop=zooD,
                      #speciesL=unique(zooD$species)[1:25],
                      measure="density",
                      meteoVars=names(meteoData1)[2:5],
                      organism="zooplankton",
                      LTERsite="Palmer")
dev.off()


#Analyses and graphs-----------------------------------------------------------------------------
birdD=merge(birdGr,meteoData1)
  
tiff("Results1/plamer_birds.tiff",unit="in",height=8,width=6.3,res=500,compression="lzw")
par(mfrow=c(4,2),mar=c(3,2.8,1.2,0.1),mgp=c(1.5,0.6,0))
resBirds=analysisRandom(pop=birdD,
                     speciesL=unique(birdD$species),measure="count",
                     meteoVars=c("yrSea","yrTemp","gsSea","gsTemp"),
                     organism="bird",LTERsite="Palmer")

#PENGUIN---------------------------------------------------------------------------------------
#Set up data
pengD=merge(pengGr,meteoData2)

resPeng=analysisRandom(pop=pengD,speciesL="AdeliePenguin",
                       meteoVars=names(meteoData2)[2:5],measure="count",
                       organism="bird",LTERsite="Palmer")

dev.off()

#Writeup results
palmerRes=rbind(resZoo,resBirds,resPeng)
write.csv(palmerRes,"Results1/resultsPalmer.csv",row.names=F)
