# Bayesian metapopulation model with grasshopper data 
setwd("C:/Users/ac79/MEGA/Projects/LTER/")
source("C:/Users/ac79/Documents/CODE/climate_vs_r/functions.R")

#Read data ---------------------------------------------------------------------------
hoppers <- read.table("Data/Sevilleta/sev106_hopperdynamics_20150826.txt",header=T,sep=",")


# Format =============================================================================

# Species count data------------------------------------------------------------------
hoppers <- mutate(hoppers, YEAR = as.numeric(substr(hoppers$PER,1,4)),
                  SEASON = substr(hoppers$PER,5,5))
hoppers <- rename(hoppers,year=YEAR,site=SITE,species=SPECIES,abund=CNT)

# experimental design
exp_des <-  distinct(select(hoppers,year,site,WEB,TRN))

# aggregate counts regardless of individual characteristics
# IGNORE SUBSTRATE for now.
hoppers <-  hoppers %>% 
            group_by(year,site,species,WEB,TRN) %>% 
            summarise(abund = sum(abund))

# Calculate species rank (by total abundance)
spp_abund <- hoppers %>% 
             group_by(species) %>% 
             summarise(total = sum(abund))
spp_rank  <- arrange(spp_abund,desc(total))

# Focus abundant species, one site only (simplicity)
trpa <- filter(hoppers, species == "TRPA" & site == "LATR" & 
                 year == 1992)

# Account for zeroes: merge with 
trpa <- merge(filter(exp_des,site == "LATR" & year == 1992),
                     trpa,all=T)
trpa$abund[is.na(trpa$abund)] <- 0  #look up "replace"


# Model abundance --------------------------------------------------------------

# matrix of observations
y <- select(trpa,WEB,TRN,abund) %>% 
        spread(TRN,abund) %>% 
          select(-WEB) %>% as.matrix()

#Format data for WinBUGS
sev_dat <- list(n      = 5, 
                J      = 6,
                y      =  y) 


# WinBUGS model
sink("Analysis1/SEV_hop_9.28.txt")
cat("
  model {

    #Priors
    beta ~ dnorm(0, 1.0E-4)
    lambda <- exp(beta)
    p ~ dunif(0,1)

    #likelihood
    for (i in 1:n) {
      N[i] ~ dpois(lambda)
      for (j in 1:J) {
        y[i,j] ~ dbin(p, N[i])
      }
    }
  }",fill = TRUE)
sink()


# Initial values
inits <- function(){ list(beta = rnorm(1, 0, 1), p = runif(1,0,1), 
                          N = rpois(5, 10)) }

# Parameters monitored
params <- c("beta","p","lambda","N")

# MCMC settings
ni <- 100000
nt <- 50
nb <- 10000
nc <- 3

# Call WinBUGS from R (BRT 215 min)
library(R2WinBUGS)
out2 <- bugs(sev_dat, inits, params, "C:/Users/ac79/Mega/Projects/LTER/Analysis1/SEV_hop_9.28.txt", 
             n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE, working.directory = getwd())
