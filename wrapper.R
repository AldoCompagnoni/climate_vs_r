#Wrapper that runs all analyses at once.
setwd("C:/Users/ac79/MEGA/Projects/LTER/")
dir="C:/Users/ac79/Documents/CODE/climate_vs_r/"


#Analyses by LTER Site 
source(paste0(dir,"santaBarbara.R"))
source(paste0(dir,"northTemperateL.R"))
source(paste0(dir,"cedarCreek.R"))
source(paste0(dir,"sevilleta.R"))
source(paste0(dir,"palmer.R"))
source(paste0(dir,"konza.R"))
#source(paste0(dir,"bonanzaCreek.R"))

#Meta analysis of results
