#This script runs preliminary analyses on LTER data
setwd("C:/Users/ac79/MEGA/Projects/LTER/")
source("C:/Users/ac79/Documents/CODE/climate_vs_r/functions.R")

# PROBLEM HERE: meteorological data starts from 2002 only.


#DATA
#Climate
SBClimate=read.csv("Data/Santa Barabara/bottom_temp_all_years_20130410.csv")

#Demography
#Fishes  http://sbc.lternet.edu/cgi-bin/showDataset.cgi?docid=knb-lter-sbc.17
fish=read.csv("Data/Santa Barabara/all_fish_all_years_20140903.csv")
#Fish 2 (Surfperch and Garibaldi) http://sbc.lternet.edu/cgi-bin/showDataset.cgi?docid=knb-lter-sbc.39
fish2=read.csv("Data/Santa Barabara/SCI_Fish_All_Years_SchmittRussell.csv")

#Invertebrates  http://sbc.lternet.edu/cgi-bin/showDataset.cgi?docid=knb-lter-sbc.19
inv=read.csv("Data/Santa Barabara/quad_swath_all_years_20140908_ReedDaniel.csv")
#Sessile organisms http://sbc.lternet.edu/cgi-bin/showDataset.cgi?docid=knb-lter-sbc.15
ses=read.csv("Data/Santa Barabara/cover_all_years_20140902_ReedDaniel.csv")
#Benthic organisms weight http://sbc.lternet.edu/cgi-bin/showDataset.cgi?docid=knb-lter-sbc.47
ben_w=read.csv("Data/Santa Barabara/SCI_Algal_Weights_All_Years_SchmittRussell.csv")
#Benthic organisms cover http://sbc.lternet.edu/cgi-bin/showDataset.cgi?docid=knb-lter-sbc.38
ben_cov=read.csv("Data/Santa Barabara/SCI_Cover_All_Years_SchmittRussell.csv")



##################################################################################################
#FORMAT
##################################################################################################

#Fishes------------------------------------------------------------------------------
fish1=fish[-which(fish$COUNT==-99999),] #Remove NAs
names(fish1)[c(5,6,7,8,10)]=c("sp_1","sp_2","sp_3","species","abund")
fish1=fish1[,c("year","site","sp_1","sp_2","sp_3","species","abund")]
fish1=subset(fish1,species != "NF") # Remove "not found" species
fish1Gr=formatDemographicSite(fish1)

#Fishes 2------------------------------------------------------------------------------
fish2=fish2[-which(fish2$COUNT==-99999),]  #Remove NAs
#Format names for function
names(fish2)[c(1,4:8)]=c("year","site","sp_1","sp_2","species","abund")
fish2Gr=formatDemographicSite(fish2)
#Put both fish datasets together!
fishGr=rbind(fish1Gr,fish2Gr)

#Invertebrates----------------------------------------------------------------------------------
inv=inv[-which(inv$COUNT==-99999),] #Remove NAs
names(inv)[c(1,4,5,6,7,8,9)]=c("year","site","sp_1","sp_2","sp_3","species","abund") #format columns
invGr=formatDemographicSite(inv[,c("year","site","sp_1","sp_2","sp_3","species","abund")])

#Sessile organisms-------------------------------------------------------------------------------
ses=ses[-which(ses$PERCENT_COVER==-99999.00),] #Remove NAs
names(ses)[c(1,4,5,6,7)]=c("year","site","sp_1","species","abund") #format columns
sesGr=formatDemographicSite(ses[,c("year","site","sp_1","species","abund")])

#Benthic organisms weight--------------------------------------------
names(ben_w)[c(1,4,6,7,8)]=c("year","site","sp_1","species","abund")
ben_wGr=formatDemographicSite(ben_w[,c("year","site","sp_1","species","abund")])

#Benthic organisms cover--------------------------------------------
#Need to summarize counts here
ben_cov=ben_cov[-which(ben_cov$SP_CODE==""),] #remove missing species
ben_covS=summarise(group_by(ben_cov,YEAR,SITE,SP_CODE),abund=n())
names(ben_covS)=c("year","site","species","abund")
ben_covGr=formatDemographicSite(ben_covS)


#Climate data---------------------------------------------------------------------------------- 
#format date
SBClimate$Date[is.na(SBClimate$Date)]="NA-NA-NA"
SBClimate$Date[SBClimate$Date==""]="NA-NA-NA"
SBClimate=cbind(formatMeteoDate(SBClimate$Date,c("year","month","day"),"-"),SBClimate)

SBClimate$temp_c[SBClimate$temp_c==-99999]=NA

#Make yearly means
climD=SBClimate[,c("temp_c","site","year","month")]
yrTemp=summarise(group_by(climD,year), yrTemp=mean(temp_c,na.rm=T))
gsTemp=summarise(group_by(subset(climD,month>4 & month<9),year), gsTemp=mean(temp_c,na.rm=T))
meteoData=merge(yrTemp,gsTemp)
meteoData=cbind(meteoData,meteoData[,2:3]^2)
names(meteoData)[4:5]=paste0(names(meteoData)[2:3],"_2")



########################################################################################################################################################
#ANALYSIS############################################################################
########################################################################################################################################################


#Fishes----------------------------------------------------------------------------------------
fishD=merge(fishGr,meteoData)

tiff("Results1/SB_Fishes.tiff",unit="in",width=6.3,height=8,res=500,compression="lzw")
par(mfrow=c(6,5),mar=c(2,2,1,0.2),mgp=c(1,0.5,0))
resFish=analysisRandom(pop=fishD,
                       meteoVars=names(meteoData)[2:3],
                       measure="count",
                       organism="fish",
                       LTERsite="SantaBarbara")
dev.off()

#Invertebrates---------------------------------------------------------------------------------
invD=merge(invGr,meteoData)
inv[inv$species %in% unique(invD$species)[21],]$TAXON_PHYLUM="Mollusca"
phylum=NULL ; for(t in 1:length(unique(invD$species))){ phylum=c(phylum,unique(subset(inv,species==unique(invD$species)[t])$inv[inv$species %in% unique(invD$species)[21],]$TAXON_PHYLUM)) }

# remove a species whose models do not converge
invD <- subset(invD,species != "PAST")

tiff("Results1/SB_Invertebrates.tiff",unit="in",width=6.3,height=9,res=500,compression="lzw")
par(mfrow=c(9,5),mar=c(2,2,1,0.2),mgp=c(1,0.5,0))
resInv=analysisRandom(pop=invD,
                      meteoVars=names(meteoData)[2:3],
                      measure="count",
                      organism=phylum,
                      LTERsite="SantaBarbara")
dev.off()

#sessile organisms--------------------------------------------------------------------------
sesD=merge(sesGr,meteoData)
phylum=NULL ; for(t in 1:length(unique(sesD$species))){ phylum=c(phylum,unique(subset(ses,species==unique(sesD$species)[t])$TAXON_PHYLUM)) }

tiff("Results1/SB_Sessiles.tiff",unit="in",width=6.3,height=9,res=500,compression="lzw")
par(mfrow=c(10,5),mar=c(2,2,1,0.2),mgp=c(1,0.5,0))
resSes=analysisRandom(pop=sesD,meteoVars=names(meteoData)[2:3],measure="cover",
                      organism=phylum,LTERsite="SantaBarbara")
dev.off()

#Benthic weight------------------------------------------------------------------------------
ben_wD=merge(ben_wGr,meteoData)

#Do NOT REPORT GTHESE RESULTS (re-check later)
#par(mfrow=c(3,2),mar=c(2,2,1,0.2),mgp=c(1,0.5,0))
#resBen_W=analysisRandom(pop=ben_wD,
#                      meteoVars=names(meteoData)[2:3],measure="weight",
#                      organism="algae",LTERsite="SantaBarbara")

#Benthic cover------------------------------------------------------------------------------
ben_covD=merge(ben_covGr,meteoData)
phylum=NULL ; for(t in 1:length(unique(ben_covD$species))){ phylum=c(phylum,unique(subset(ben_cov,SP_CODE==unique(ben_covD$species)[t])$TAXON_PHYLUM)) }
phylum=phylum[-26]
# Remove species whose model does not converge
ben_covD=subset(ben_covD,species != "Lf")


tiff("Results1/SB_Algae.tiff",unit="in",width=6.3,height=9,res=500,compression="lzw")
par(mfrow=c(7,4),mar=c(2,2,1,0.2),mgp=c(1,0.5,0))
resAlgae=analysisRandom(pop=ben_covD, #26 and 27 does not converge
                        meteoVars=names(meteoData)[2:3],
                        measure="cover",
                        organism="algae",LTERsite="SantaBarbara")
dev.off()


#write out
results=rbind(resFish,resInv,resSes,resAlgae) #resBen_W
write.csv(results,"Results1/resultsSantaBarbara.csv",row.names=F)
  