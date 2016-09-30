# Bayesian metapopulation model with grasshopper data 
setwd("C:/Users/ac79/MEGA/Projects/LTER/")
gitDir="C:/Users/ac79/Documents/CODE/climate_vs_r/"
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


# seasonsal trends in TRPA
hoppers_season <-  hoppers %>% 
        group_by(year,SEASON,site,species) %>% 
              summarise(abund = sum(abund))

tmp <- filter(hoppers_season, species == "TRPA" & site == "LATR")
tmp <- tmp %>% select(year,SEASON,abund) %>% spread(SEASON,abund,fill=0) %>% data.frame() 
tmp <- tmp  %>% select(year,E,L) %>% as.matrix()
matplot(tmp[,1],tmp[,-1],type = "l")
# low correlations between seasons
cor(tmp[,-1])

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
sink(paste0(gitDir,"SEV_hop_9.28.txt"))
cat("
  model {

    #Priors
    #beta ~ dnorm(0, 1.0E-4)
    #lambda <- exp(beta)
    p ~ dunif(0,1)

    # WEB level means 
    for (k in 1:n) { beta0[k] ~ dnorm(0,.01) }

    #likelihood
    for (i in 1:n) {

      N[i] ~ dpois(lambda[i])
      log(lambda[i]) <- beta0[i]

      for (j in 1:J) {
        y[i,j] ~ dbin(p, N[i])
      }
    }
  }",fill = TRUE)
sink()


# Initial values
inits <- function(){ list(beta0 = rnorm(1, 0, 1), p = runif(1,0,1), 
                          N = rpois(5, 10)) }

# Parameters monitored
params <- c("beta0","p","lambda","N")

# MCMC settings
ni <- 100000
nt <- 50
nb <- 10000
nc <- 3

# Call WinBUGS from R (BRT 215 min)
library(R2WinBUGS)
out2 <- bugs(sev_dat, inits, params, paste0(gitDir,"SEV_hop_9.28.txt"), 
             n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE, working.directory = getwd())


# Model abundance by year -------------------------------------------------------

# Data from all years
trpa <- filter(hoppers, species == "TRPA" & site == "LATR")

# Account for zeroes: merge with 
trpa <- merge(filter(exp_des,site == "LATR"),trpa,all=T)
trpa$abund[is.na(trpa$abund)] <- 0  #look up "replace"


# matrix of observations
obs <- select(trpa,year,WEB,TRN,abund)
#obs <- mutate(obs, year_web = paste(year,WEB,sep="_") )

# Site index
site_i <- obs$WEB       
# Number of sites
n_s    <- 5 
# Year index
year_i <- obs$year - 1991
# Number of years
n_y    <- length(unique(year_i))
# Abundance data
y      <- obs$abund

# Random
re_i   <- 1:110

#Format data for WinBUGS
sev_dat <- list(n      = length(y), 
                site_i = site_i, n_s = n_s,
                year_i = year_i, n_y = n_y,
                y      =  y) 


# WinBUGS model 
sink(paste0(gitDir,"SEV_hop_9.28.txt"))
cat("
  model {

    # Hyperpriors - year random effects

    # Priors
    # WEB level means 
    for (k in 1:n_s) { beta0[k] ~ dnorm(0,.01) }
    
    # Year random effects
    for (i in 1:n_y){ beta1[i] ~ dnorm(0, 0.001) }

    # Observation error
    p ~ dunif(0,1)

    # Models ----------------------------------------------
    # Abundance model
    for (i in 1:n) {

      N[i] ~ dpois(lambda[i])
      log(lambda[i]) <- beta0[site_i[i]] + beta1[year_i[i]]  
    
      # Observation model
      y[i] ~ dbin(p, N[i])

    }    

  }",fill = TRUE)
sink()



# Initial values
inits <- function(){ list(beta1 = rnorm(n_y, 0, 1), beta0 = rnorm(n_s, 0, 1),
                          p = runif(1,0,1), N = rpois(sev_dat$n, 30))
                          } 

# Parameters monitored
params <- c("beta1","beta0","p","N")

# MCMC settings
ni <- 100000
nt <- 500
nb <- 10000
nc <- 3

# Call WinBUGS from R (BRT 215 min)
library(R2WinBUGS)
out2 <- bugs(sev_dat, inits, params, paste0(gitDir,"SEV_hop_9.28.txt"), 
             n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, debug = T, working.directory = getwd())

