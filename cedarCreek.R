#THE PLAN:
setwd("C:/Users/ac79/MEGA/Projects/RICE/LTER/")
source("Analysis1/functions.R")

#http://www.cedarcreek.umn.edu/research/data?EXPERIMENT_SEARCH=%25&search_words=&CORE_AREA=%25&RESEARCHERS_LNAME=&SEARCH=Search
#Climate Data
climateCC=read.csv("Data/Cedar Creek/DailyClimateSummary.csv")

#Demographic data
plantBioCC=read.csv("Data/Cedar Creek/CedarCreekBiomass.csv")

#Arthropods (this does not have 10 years of continuous data, unfortunately)
#insectCC=read.table("Data/Cedar Creek/e120_Main Plots All Arthropod Insect Sweepnet Sampling_DavidTilman.txt",header=T,sep="\t")

#Grasshopper1
hopper=read.table("Data/Cedar Creek/e014_Core Old Field Grasshopper Sampling .txt",header=T,sep="\t")


#####################################################################################################
#Format plant biomass data##########################################################################
#####################################################################################################

#Biomass data-----------------------------------------------------------------------------------------
#Remove nitrogen addition plots
ccDat=subset(plantBioCC,NitrAdd==0)
names(ccDat)[c(2,3,9,10)]=c("year","site","species","abund")
bioGr=formatDemographicSite(ccDat[,c(2,3,9,10)])

plantList=summarise(group_by(ccDat,species),abundance=sum(abund))
plantList=plantList$species[order(plantList$abundance,decreasing=T)]

#Grasshoppers-------------------------------------------------------------------
hopper$species=paste(hopper$Genus,hopper$Specific.epithet,sep="_")
names(hopper)[c(3,4,11,12)]=c("site","year","abund","species")
hopperGr=formatDemographicSite(hopper[,c(4,3,12,11)])


#####################################################################################################
#FORMAT CLIMATE DATA#################################################################################
#####################################################################################################
climateCC=select(climateCC,Date,MaxTemp.degF.,MinTemp.degF.,Precip.inches.)
climateCC=cbind(formatMeteoDate(climateCC$Date,format=c("month","day","Year"),"/"),climateCC)
climateCC$meanTemp=(climateCC$MaxTemp.degF. + climateCC$MinTemp.degF.)/2

#LAG year summaries
yrTemp=summarise(group_by(climateCC,Year),yrTemp=mean(meanTemp))
yrPrec=summarise(group_by(climateCC,Year),yrPrec=sum(Precip.inches.))
#this makes it "the previous year"
yrTemp$Year=yrTemp$Year-1 ; yrPrec$Year=yrPrec$Year-1

#Two options for CURRENT growing season summaries (May-August,April-July)
gs1=subset(climateCC,month > 4 & month < 9)
gs2=subset(climateCC,month > 3 & month < 8)
#temperatures
gs1Temp=summarise(group_by(gs1,Year),gs1Temp=mean(meanTemp))
gs2Temp=summarise(group_by(gs2,Year),gs2Temp=mean(meanTemp))
#Precipitations
gs1Prec=summarise(group_by(gs1,Year),gs1Prec=sum(Precip.inches.))
gs2Prec=summarise(group_by(gs2,Year),gs2Prec=sum(Precip.inches.))

#final list
meteoList=list(yrTemp,yrPrec,gs1Temp,gs2Temp,gs1Prec,gs2Prec)
meteoData=Reduce(function(...) merge(...,all=T),meteoList)
#quadratic terms
quadratic=meteoData[,-1]*meteoData[,-1]
names(quadratic)=paste0(names(quadratic),"_2")
meteoData=cbind(meteoData,quadratic)
names(meteoData)[1]="year"

#####################################################################################################
#ANALYSIS############################################################################
#####################################################################################################

#plant biomass-------------------------------------------------------------------------
ccD=merge(bioGr,meteoData)
ccD=subset(ccD,species!="Miscellaneous litter")


#grasshopper---------------------------------------------------------------------------
hopperD=merge(hopperGr,meteoData)

tiff("Results1/cedarCreek_hopper.tiff",unit="in",height=5,width=6.3,res=500,compression="lzw")

par(mfrow=c(5,5),mar=c(2.3,2.5,1.2,0.1),mgp=c(1.3,0.5,0),oma=c(0,0,0,0.5))
resHopp=analysisRandom(pop=hopperD,speciesL=unique(hopperD$species),
                   meteoVars=names(meteoData)[2:7],
                   organism="arhtopod",measure="count",
                   LTERsite="CedarCreek")
dev.off()

#EIGHTEEN SPECIES WITH only TEMPORAL REPLICATION
ccD=merge(ccD,meteoData)
tiff("Results1/cedarCreek_plant.tiff",unit="in",height=8,width=6.3,res=500,compression="lzw")
par(mfrow=c(5,4),mar=c(2.3,2.5,1.2,0.1),mgp=c(1.3,0.5,0))
resCC=analysisRandom(pop=ccD,speciesL=plantList[c(2:3,5:length(plantList))],
                     meteoVars=names(meteoData)[2:7],measure="biomass",
                     organism="plant",LTERsite="CedarCreek")
dev.off()

tiff("Results1/PoaPratensis.tiff",unit="in",height=6.3,width=6.3,res=500,compression="lzw")
par(mfrow=c(1,1),mar=c(2.3,2.5,1.2,0.1),mgp=c(1.3,0.5,0))
resCC=analysisRandom(pop=ccD,speciesL=plantList[3],
                     meteoVars=names(meteoData)[2:7],measure="biomass",
                     organism="plant",LTERsite="CedarCreek")
dev.off()


tiff("Results1/cedarCreek_plant.tiff",unit="in",height=8,width=6.3,res=500,compression="lzw")
par(mfrow=c(5,4),mar=c(2.3,2.5,1.2,0.1),mgp=c(1.3,0.5,0))
resCC=analysisRandom(pop=ccD,speciesL=plantList[c(2:3,5:length(plantList))],
                     meteoVars=names(meteoData)[2:7],measure="biomass",
                     organism="plant",LTERsite="CedarCreek")


resCC=analysisRandom(pop=ccD,speciesL=plantList[2],
                     meteoVars=names(meteoData)[2:7],measure="biomass",
                     organism="plant",LTERsite="CedarCreek")


resCC=analysisRandom(pop=ccD,speciesL=plantList[2:3],
                     meteoVars=names(meteoData)[2:7],measure="biomass",
                     organism="plant",LTERsite="CedarCreek")


#store results
results=rbind(resCC,resHopp)
write.csv(resCC,"Results1/results_CedarCreek.csv",row.names=F)






#Find species with 30 years of data available---------------------------------------------------------
#yrSppBio=aggregate(Biomass ~ Species + Year,sum,data=ccDat)
#sppList=unique(yrSppBio$Species)
#allSppKeep=NULL
#for(i in 1:length(sppList)){
#  tmp=subset(yrSppBio,Species==sppList[i])
#  if( length(unique(tmp$Year)) > 29)  allSppKeep=c(allSppKeep,sppList[i])
#}
#18 (eighteen!!!) species have at least 30 years of data
#sppKeep=setdiff(allSppKeep,c("Poa pratensis","Schizachyrium scoparium"))
#bioCC=aggregate(Biomass ~ Species + Year,sum,data=ccDat)
#calculate growth rates
#nonRepSppD=list()
#for(i in 1:length(sppKeep)){
#  tmp=subset(bioCC,Species==sppKeep[i])
#  tmp$logGr=logGr(tmp$Biomass)
#  nonRepSppD[[i]]=tmp  
#}
#nonRepD=Reduce(function(...) rbind(...),nonRepSppD)


#Find species with multiple fields AND plots available-------------------------------------------
#repLength=subset(aggregate(Biomass ~ Species + Field + Plot,length,data=ccDat),Biomass==30)
#only two species which have 
#sppList=c("Poa pratensis","Schizachyrium scoparium")
#Test: both species just miss year 2013
#for(i in 1:length(sppList)){
#  tmp=subset(repLength,Species==sppList[i])
#  print(sppList[i])
#  for(i in 1:nrow(tmp)){
#    yearTmp=subset(ccDat,Species==tmp$Species[i] & Field ==tmp$Field[i] & Plot==tmp$Plot[i])
#    print(sum(unique(yearTmp$Year)==2013))
#  }
#}


#POPR file (replication==8)
#select only the spatially replicated plots
#fpPopr=filter(repLength,Species=="Poa pratensis") 
#poprList=list()
#for(i in 1:nrow(fpPopr)){ #Calculate growth rates for each plot
#  poprList[[i]]=subset(ccDat,Species==fpPopr$Species[i] & Field==fpPopr$Field[i] & Plot==fpPopr$Plot[i])
#  poprList[[i]]=select(poprList[[i]],Year,Field,Plot,Species,Biomass)
#  poprList[[i]]$logGr=logGr(poprList[[i]]$Biomass)
#}
#popr=Reduce(function(...) rbind(...),poprList)


#SCSC file (replication==7)
#select only the spatially replicated plots
#fpScsc=filter(repLength,Species=="Schizachyrium scoparium")
#scscList=list()
#for(i in 1:nrow(fpScsc)){ #Calculate growth rates for each plot
#  scscList[[i]]=subset(ccDat,Species==fpScsc$Species[i] & Field==fpScsc$Field[i] & Plot==fpScsc$Plot[i])
#  scscList[[i]]=select(scscList[[i]],Year,Field,Plot,Species,Biomass)
#  scscList[[i]]$logGr=logGr(scscList[[i]]$Biomass)
#}
#scsc=Reduce(function(...) rbind(...),scscList)

#Spat replication spp
#repD=rbind(popr,scsc)
#repD=repD[,c("Year","Species","Biomass","logGr")]

#All Population data
#popD=rbind(repD,nonRepD)[,c("Year","Species","logGr")]
#names(popD)=c("year","species","logGr")
#popD=popD[!is.na(popD$logGr),] #Remove NAs

