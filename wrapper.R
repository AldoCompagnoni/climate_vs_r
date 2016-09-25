#Wrapper that runs all analyses at once.
setwd("C:/Users/ac79/MEGA/Projects/RICE/LTER/")

#Analyses by LTER Site 
source("Analysis1/santaBarbara.R")
source("Analysis1/northTemperateL.R")
source("Analysis1/cedarCreek.R")
source("Analysis1/sevilleta.R")
source("Analysis1/palmer.R")
source("Analysis1/konza.R")
source("Analysis1/bonanzaCreek.R")

#Meta analysis of results
