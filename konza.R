setwd("C:/Users/ac79/MEGA/Projects/LTER/")
source("C:/Users/ac79/Documents/CODE/climate_vs_r/functions.R")

#Data
#climatic
climK=read.table("Data/Konza/dailyMeteo_AWE012.txt",header=T,sep=",")

#Population
#vegetation
#biomass=read.table("Data/Konza/veg_PAB011_AlanKnapp.txt",header=T,sep=",")
#Cover
cov=read.table("Data/Konza/vegCover_PVC021_David C. Hartnett.txt",header=T,sep=",")

#grasshoppers
hoppers=read.table("Data/Konza/GrassHopp_CGR022_indiv_Anthony Joern.txt",header=T,sep=",")

#Fishes
fish=read.table("Data/Konza/Fish_CFP012_KeithGido.txt",header=T,sep=",")

#Small mammals
mamm=read.table("Data/Konza/smallMamm_CSM011_DonaldKaufman.txt",header=T,sep=",")


###########################################################################
#FORMAT
###########################################################################

#Vegetation cover----------------------------------------------------------
cov$Cover=factor(cov$Cover,labels=paste0(c(0.5,3,15,37.5,62.5,85,97.5)))
cov$Cover=as.numeric(as.character(cov$Cover))
#Find watersheds with CONTINUOUS DATA
tempRepW=table(select(cov,RECYEAR,WATERSHED))
watershedKeep=colnames(tempRepW)[apply(tempRepW==0,2,sum)==0]
wI=which(cov$WATERSHED %in% watershedKeep)
cov=cov[wI,]
cov$species=paste(cov$AB_GENUS,cov$AB_SPECIES,sep="_")
names(cov)[c(3,6,14)]=c("year","site","abund")
covGr=formatDemographicSite(cov[,c("year","site","species","abund")])

#species rank
sppRank=summarise(group_by(cov[,c("year","site","species","abund")],
                       species),abundance=sum(abund))
sppRank$species=sppRank$species[order(sppRank$abundance,decreasing=T)]

#GrassHoppers-------------------------------------------------------------------
#mistake
hoppers[which(hoppers$TOTAL=="1 01"),]$TOTAL="1"
hoppers$TOTAL[hoppers$TOTAL==""]=NA
hoppers$TOTAL=as.numeric(hoppers$TOTAL)
#WHAT IS REPSITE?!?
hopGr=summarise(group_by(hoppers,RECYEAR,WATERSHED,SPECIES),
                count=sum(TOTAL))
names(hopGr)=c("year","site","species","abund")
hopGr=formatDemographicSite(hopGr)


#fishes-------------------------------------------------------------------------------------
tmpMat=matrix(unlist(strsplit(fish$recdate,"/")),length(fish$recdate),3,byrow=T)
tmpMat[,3]=matrix(unlist(strsplit(tmpMat[,3]," ")),nrow(tmpMat),3,byrow=T)[,1]
colnames(tmpMat)=c("month","day","year")
tmpMat=as.data.frame(tmpMat)
for(i in 1:3){ tmpMat[,i]=as.numeric(tmpMat[,i]) }
fish=cbind(tmpMat,fish)

fishGr=melt(fish[,c(3,7,10:28)],id.vars=c("year","watershed"))
names(fishGr)=c("year","site","species","abund")
fishGr$species = as.character(fishGr$species)
fishGr=formatDemographicSite(fishGr)


#small mammals-----------------------------------------------------------------
mammGr=melt(mamm[,c(4,6,7:21)],id.vars=c("RECYEAR","WATERSHED.LINE"))
names(mammGr)=c("year","site","species","abund")
mammGr$species = as.character(mammGr$species)
mammGr=formatDemographicSite(mammGr)


#Format climate data---------------------------------------------------------------------------

#introduce NAs and convert to numeric
#maximum temperature
climK$TMAX[climK$TMAX==""]=NA
climK$TMAX[climK$TMAX=="."]=NA
#minimum temperature
climK$TMIN[climK$TMIN==""]=NA
climK$TMIN[climK$TMIN=="."]=NA
#precipitation
climK$DPPT[climK$DPPT=="."]=NA
climK$DPPT[climK$DPPT==""]=NA
#convert to numeric
climK$DPPT=as.numeric(climK$DPPT)
climK$TMIN=as.numeric(climK$TMIN)
climK$TMAX=as.numeric(climK$TMAX)
climK[which(climK$TMIN< -50),][,c("TMAX","TMIN")]=NA
#mean Temperature
climK$TMEAN=(climK$TMIN+climK$TMAX)/2

#Summary statistics
#Year data
yrClim=summarise(group_by(climK,RECYEAR),yrTemp=mean(TMEAN,na.rm=T),yrPrec=sum(DPPT,na.rm=T))
yrClim$RECYEAR=yrClim$RECYEAR + 1
#growing seasons, two possibilities (April-September,April-July)
gs1D=subset(climK,RECMONTH>3 & RECMONTH < 9)
gs2D=subset(climK,RECMONTH>3 & RECMONTH < 8)
gs1=summarise(group_by(gs1D,RECYEAR),gs1Temp=mean(TMEAN,na.rm=T),gs1Prec=sum(DPPT,na.rm=T))
gs2=summarise(group_by(gs2D,RECYEAR),gs2Temp=mean(TMEAN,na.rm=T),gs2Prec=sum(DPPT,na.rm=T))

meteoList=list(yrClim,gs1,gs2)
meteoData1=Reduce(function(...) merge(...),meteoList)
meteoData2=meteoData1[,-1]^2
colnames(meteoData2)=paste0(names(meteoData1[,-1]),"_2")
meteoData=cbind(meteoData1,meteoData2)
names(meteoData)[1]="year"

#############################################################################
#ANALYSIS###############
#############################################################################

#Grasshoppers----------------------------------------------------------------
hopD=merge(hopGr,meteoData)
#names(hopKonzaD)[1:2]=c("year","species")

tiff("Results1/konza_grasshopper.tiff",unit="in",height=5,width=6.3,res=500,compression="lzw")

par(mfrow=c(5,3),mar=c(2.3,2.5,1.2,0.1),mgp=c(1.3,0.5,0),oma=c(0,0,0,0.5))
resHopper=analysisRandom(pop=hopD,speciesL=unique(hopD$species),
                       meteoVars=names(meteoData)[2:7],measure="count",
                       organism="arthropod",LTERsite="Konza")
dev.off()

#Vegetation cover data---------------------------------------------------------------------------------
covD=merge(covGr,meteoData)

tiff("Results1/konza_vegetation1.tiff",unit="in",height=5,width=6.3,res=500,compression="lzw")

par(mfrow=c(5,5),mar=c(2.3,2.5,1.2,0.2),mgp=c(1.3,0.5,0),oma=c(0,0,0,0.5))
resVegCover1=analysisRandom(pop=covD,speciesL=sppRank$species[1:25],
                            meteoVars=names(meteoData)[2:7],measure="cover",
                            organism="plant",LTERsite="Konza")
dev.off()

# all results
resVegCover=analysisRandom(pop=covD,
                           meteoVars=names(meteoData)[2:7],measure="cover",
                           organism="plant",LTERsite="Konza")

#tiff("Results1/konza_vegetation2.tiff",unit="in",height=8,width=6.3,res=500,compression="lzw")

#par(mfrow=c(5,5),mar=c(2.3,2.5,1.2,0.1),mgp=c(1.3,0.5,0))
#resVegCover2=analysisRandom(pop=covD,speciesL=sppRank$species[c(26:28,30:50)],
#                            meteoVars=names(meteoData)[2:7],measure="cover",
#                            organism="plant",LTERsite="Konza")
#dev.off()
#tiff("Results1/konza_vegetation3.tiff",unit="in",height=8,width=6.3,res=500,compression="lzw")

#par(mfrow=c(5,5),mar=c(2.3,2.5,1.2,0.1),mgp=c(1.3,0.5,0))
#resVegCover3=analysisRandom(pop=covD,speciesL=sppRank$species[c(51:54,56:65,67:75)],
#                            meteoVars=names(meteoData)[2:7],measure="cover",
#                            organism="plant",LTERsite="Konza")
#dev.off()
#tiff("Results1/konza_vegetation4.tiff",unit="in",height=8,width=6.3,res=500,compression="lzw")#

#par(mfrow=c(5,5),mar=c(2.3,2.5,1.2,0.1),mgp=c(1.3,0.5,0))
#resVegCover4=analysisRandom(pop=covD,speciesL=sppRank$species[c(76:98,100)],#99 does not converge!
#                            meteoVars=names(meteoData)[2:7],measure="cover",
#                            organism="plant",LTERsite="Konza")
#dev.off()
#tiff("Results1/konza_vegetation5.tiff",unit="in",height=8,width=6.3,res=500,compression="lzw")

#par(mfrow=c(5,5),mar=c(2.3,2.5,1.2,0.1),mgp=c(1.3,0.5,0))
#resVegCover5=analysisRandom(pop=covD,speciesL=sppRank$species[c(101:132)],
#                            meteoVars=names(meteoData)[2:7],measure="cover",
#                            organism="plant",LTERsite="Konza")
#dev.off()


# Fish data---------------------------------------------------------------------------------
fishD=merge(fishGr,meteoData)

tiff("Results1/konza_fish.tiff",unit="in",height=3.15,width=6.3,res=500,compression="lzw")

par(mfrow=c(1,2),mar=c(2.3,2.5,1.2,0.1),mgp=c(1.3,0.5,0))
resFish=analysisRandom(pop=fishD,speciesL=unique(fishD$species),
                     meteoVars=names(meteoData)[2:7],measure="count",
                     organism="fish",LTERsite="Konza")
dev.off()


#Small mammals data---------------------------------------------------------------------------------
mammD=merge(mammGr,meteoData)

tiff("Results1/konza_smallMammals.tiff",unit="in",height=8,width=6.3,res=500,compression="lzw")

par(mfrow=c(4,2),mar=c(2.3,2.5,1.2,0.1),mgp=c(1.3,0.5,0))
resMamm=analysisRandom(pop=mammD,speciesL=unique(mammD$species),
                     meteoVars=names(meteoData)[2:7],measure="count",
                     organism="mammal",LTERsite="Konza")
dev.off()


#write out all results
results=rbind(resHopper,resVegCover,resFish,resMamm)
write.csv(results,"Results1/resultsKonza.csv",row.names=F)








#Vegetation cover----------------------------------------------------------
#Convert cover values to cover abundances
#categories are 1=0.5,2=3,3=15,4=37.5,5=62.5,6=85,7=97.5
#cov$Cover=factor(cov$Cover,labels=paste0(c(0.5,3,15,37.5,62.5,85,97.5)))
#cov$Cover=as.numeric(as.character(cov$Cover))

#Find watersheds with CONTINUOUS DATA
#tempRepW=table(select(cov,RECYEAR,WATERSHED))
#watershedKeep=colnames(tempRepW)[apply(tempRepW==0,2,sum)==0]
#wI=which(cov$WATERSHED %in% watershedKeep)
#cov=cov[wI,]

#Mean covers across watershed, every year, per species
#meanCov=summarise(group_by(cov,SPECODE,RECYEAR,WATERSHED),Cover=mean(Cover,na.rm=T))
#meanCov$notNa=as.numeric(!is.na(meanCov$Cover))
#counts=summarise(group_by(meanCov,SPECODE),count=sum(notNa))
#What species have replication 4 across all 33 years? 
#sppCodList=counts$SPECODE[counts$count==132]


#Population data
#popGr=list()
#for(i in 1:length(sppCodList)){

#  tmp=subset(meanCov,SPECODE==sppCodList[i])
#  wideForm=dcast(tmp,RECYEAR + SPECODE ~ WATERSHED,value.var = "Cover")
#  logGrowthRate=NULL
#  for(w in 1:length(watershedKeep)){
#    logGrowthRate=cbind(logGrowthRate,logGr(wideForm[,watershedKeep[w]]))
#  }
#  colnames(logGrowthRate)=watershedKeep
#  popGr[[i]]=cbind(wideForm$RECYEAR,wideForm$SPECODE,logGrowthRate)

#}
#grDf=as.data.frame(Reduce(function(...) rbind(...),popGr))
#names(grDf)[1:2]=c("year","species")
#grLong=melt(grDf,id=c("year","species"))[,-3]
#names(grLong)[3]="logGr"
#Remove NA values
#grLong=grLong[-which(is.na(grLong$logGr)),]


#GrassHoppers-------------------------------------------------------------------
#mistake
#hoppers[which(hoppers$TOTAL=="1 01"),]$TOTAL="1"
#hoppers$TOTAL[hoppers$TOTAL==""]=NA
#hoppers$TOTAL=as.numeric(hoppers$TOTAL)
#WHAT IS REPSITE?!?
#hoppers=summarise(group_by(hoppers,RECYEAR,WATERSHED,SPECIES),
#                   count=sum(TOTAL))
#hopKonza=summarise(group_by(hoppers,RECYEAR,SPECIES),
#                   count=sum(TOTAL))

#spp with at least 24 years of data
#presence=apply(table(hopKonza[,1:3]),2,sum)
#spp24=names(presence)[presence==24]
#spp23=names(presence)[presence==23]

#Two continuous periods for spp24: 1982:1991 ; 1996:2008
#Three continuous periods for spp23: 1982:1988, 1990:1991 ; 1996:2008
#sppSpecD=as.data.frame(subset(hopKonza,SPECIES==spp24[1]))
#period24=list(first=c(1982:1991),second=c(1996:2008))
#period23=list(first=c(1982:1988),second=c(1990:1991),third=c(1996:2008))

#species with 24 years
#gs24=NULL
#for(i in 1:length(spp24)){
#  sppTmp=as.data.frame(subset(hopKonza,SPECIES==spp24[i]))
#  perGr=list(NULL)
#  for(p in 1:length(period24)){
#    perGr[[p]]=sppTmp[which(sppTmp$RECYEAR %in% period24[[p]]),]
#    perGr[[p]]$logGr=logGr(perGr[[p]]$count)
#  } 
#  gs24=rbind(gs24,Reduce(function(...) rbind(...),perGr))
#} 

#species with 23 years
#gs23=NULL
#for(i in 1:length(spp23)){
#  sppTmp=as.data.frame(subset(hopKonza,SPECIES==spp23[i]))
#  perGr=list(NULL)
#  for(p in 1:length(period23)){
#    perGr[[p]]=sppTmp[which(sppTmp$RECYEAR %in% period23[[p]]),]
#    perGr[[p]]$logGr=logGr(perGr[[p]]$count)
#  } 
#  gs23=rbind(gs23,Reduce(function(...) rbind(...),perGr))
#} 
#hopKonzaGr=rbind(gs24,gs23)[,c("RECYEAR","SPECIES","logGr")]


#fish--------------------------------------------------------------------------------
#tmpMat=matrix(unlist(strsplit(fish$recdate,"/")),length(fish$recdate),3,byrow=T)
#tmpMat[,3]=matrix(unlist(strsplit(tmpMat[,3]," ")),nrow(tmpMat),3,byrow=T)[,1]
#colnames(tmpMat)=c("month","day","year")
#tmpMat=as.data.frame(tmpMat)
#for(i in 1:3){ tmpMat[,i]=as.numeric(tmpMat[,i]) }
#fish=cbind(tmpMat,fish)

#Means
#fishM=as.data.frame(summarise(group_by(fish,year),CAMANO=sum(CAMANO,na.rm=T),
#                CATCOM=sum(CATCOM,na.rm=T),CYPLUT=sum(CYPLUT,na.rm=T),
#                ETHNIG=sum(ETHNIG,na.rm=T),ETHSPE=sum(ETHSPE,na.rm=T),
#                GAMAFF=sum(GAMAFF,na.rm=T),LEPCYA=sum(LEPCYA,na.rm=T),
#                LEPHUM=sum(LEPHUM,na.rm=T),LEPMAC=sum(LEPMAC,na.rm=T),
#                LEPMEG=sum(LEPMEG,na.rm=T),LUXCOR=sum(LUXCOR,na.rm=T),
#                NOTEXI=sum(NOTEXI,na.rm=T),NOTSTR=sum(NOTSTR,na.rm=T),
#                PHEMIR=sum(PHEMIR,na.rm=T),PHOERY=sum(PHOERY,na.rm=T),
#                PIMNOT=sum(PIMNOT,na.rm=T),PIMPRO=sum(PIMPRO,na.rm=T),
#                SEMATR=sum(SEMATR,na.rm=T),MOXMIC=sum(MOXMIC,na.rm=T)))

#Select fishes
#fishList=names(fishM)[apply(fishM==0,2,sum)==0][-1]

#fishGr=list()
#for(i in 1:length(fishList)){

#  tmp=logGr(fishM[,c(fishList[i])])
#  fishGr[[i]]=as.data.frame(cbind(fishM$year,tmp))
#  names(fishGr[[i]])=c("year",fishList[i])

#}
#fishGr=Reduce(function(...) merge(...),fishGr)
#fishGr=fishGr[-1,]
#fishGr=melt(fishGr,id.var="year")  
#names(fishGr)[2:3]=c("species","logGr")



#small mammals-----------------------------------------------------------------

#yearly sums 
#mammTmp=mamm[,c(4,7:21)] 
#mammS=as.data.frame(mammTmp %>% group_by(RECYEAR) %>% summarise_each(funs(sum)))
#names(mammS)[1]="year"

#Select species, and prepare pop. growth rates
#mammList=names(mammS[,-1])[which(apply(mammS[,-1]==0,2,sum)<8)]

#Population growth rates
#mammGr=list()
#for(i in 1:length(mammList)){

#  tmp=logGr(mammS[,c(mammList[i])])
#clean up "bad data"
#  tmp[tmp==Inf]=NA
#  tmp[tmp==-Inf]=NA
#  tmp[is.nan(tmp)]=NA
#Finish up
#  mammGr[[i]]=as.data.frame(cbind(mammS$year,tmp))
#  names(mammGr[[i]])=c("year",mammList[i])

#}
#mammGr=Reduce(function(...) merge(...),mammGr)
#mammGr=melt(mammGr,id.var="year")  
#names(mammGr)[2:3]=c("species","logGr")
#mammGr=mammGr[which(!is.na(mammGr$logGr)),] #Remove NAs
