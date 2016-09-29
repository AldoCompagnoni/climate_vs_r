#THE PLAN:
#Climate data: http://sev.lternet.edu/data/sev-1
#Demographic data: http://sev.lternet.edu/data/sev-106
setwd("C:/Users/ac79/Mega/Projects/RICE/LTER/")
options(stringsAsFactors = F)
library(lubridate)
library(bbmle)
source("Analysis/AICtable.R")

#Data
GrasshopperSEV=read.table("Data/Sevilleta/sev106_hopperdynamics_20150826.txt",header=T,sep=",")
climateSEV=read.csv("Data/Sevilleta/sev1_meteorology_20140111.csv")


#Format
#Species count data
GrasshopperSEV$YEAR=as.numeric(substr(GrasshopperSEV$PER,1,4))
GrasshopperSEV$SEASON=substr(GrasshopperSEV$PER,5,5)


#Climate data
climateSEV$Temp_C[climateSEV$Temp_C== -999]=NA
climateSEV$Temp_C[climateSEV$Temp_C < -39.9]=NA
climateSEV$Precip[climateSEV$Precip== -999]=NA
climateSEV$Precip[climateSEV$Precip<0]=NA
climateSEV=subset(climateSEV,Station_ID == 42 | Station_ID == 49 | Station_ID == 50)#exclude sites we don't care about!

#Daily means/sums
dailyTemp=aggregate(Temp_C ~ Jul_Day + Year + Station_ID,mean,data=climateSEV,na.rm=T) #Average across all sites
dailyPrecip=aggregate(Precip ~ Jul_Day + Year + Station_ID,sum,data=climateSEV,na.rm=T)
dailyPrecip=aggregate(Precip ~ Jul_Day + Year + Station_ID,mean,data=dailyPrecip,na.rm=T) #Average across all summed up precipitation
climateData=merge(dailyTemp,dailyPrecip)

#Spring and summer temp./precip.
sprTemp=aggregate(Temp_C ~ Year + Station_ID,mean,data=subset(climateData, Jul_Day<152)) ; names(sprTemp)[3]="sprTemp"
sprPrecip=aggregate(Precip ~ Year + Station_ID,sum,data=subset(climateData, Jul_Day<152)) ; names(sprPrecip)[3]="sprPrecip"
summTemp=aggregate(Temp_C ~ Year + Station_ID,mean,data=subset(climateData, Jul_Day<274)) ; names(summTemp)[3]="summTemp"
summPrecip=aggregate(Precip ~ Year + Station_ID,sum,data=subset(climateData, Jul_Day<274)) ; names(summPrecip)[3]="summPrecip"
tempTm1=aggregate(Temp_C ~ Year + Station_ID,mean,data=climateData) ; names(tempTm1)[3]="tempTm1"
precipTm1=aggregate(Precip ~ Year + Station_ID,sum,data=climateData) ; names(precipTm1)[3]="precipTm1"
tempTm1$Year=tempTm1$Year+1      #Associate this data with NEXT year!
precipTm1$Year=precipTm1$Year+1  #Associate this data with NEXT year!

climateList=list(sprTemp,sprPrecip,summTemp,summPrecip,tempTm1,precipTm1)
annualClimate=Reduce(function(...) merge(...),climateList)
#Let's just use the PINION-JUNIPER data (Cerro Montoso)
pjClimate=subset(annualClimate, Station_ID==42)[,-2]
names(pjClimate)[1]="YEAR"

#Create squared values
pjClimate=cbind(pjClimate,pjClimate[,-1]*pjClimate[,-1])
names(pjClimate)[8:13]=paste(names(pjClimate)[2:7],"_2",sep="")



#EXAMPLE OF BAYESIAN ANALYSIS USING SPECIES "trpa" at site "BOER"
trpaBoer=subset(GrasshopperSEV,SPECIES=="TRPA" & SITE=="BOER")
yearWebTrnGrid=expand.grid(YEAR=c(1992:2013),WEB=c(1:5),TRN=c(12,36,60,84,108,132))
bData=merge(yearWebTrnGrid,trpaBoer,all.x=T)[,c("YEAR","WEB","TRN","CNT")]
bData$CNT[is.na(bData$CNT)]=0

### Four-way unscaled Wishart
sev_dat <- list(n = nrow(bData), 
                nYears = as.integer(length(unique(bData$YEAR))), 
                nSites = as.integer(length(unique(bData$WEB))),
                years = bData$YEAR - 1991, sites = bData$WEB,
                y = bData$CNT) 


# 12.3.3. Binomial-mixture model with overdispersion in both abundance and detection
sink("Analysis/SevBUGS_simple.txt")
cat("
    model{
      
      # Priors
      for(ss in 1:nSites){
        bSite[ss]  ~ dnorm(0, tau.site)
      }
      for(yy in 1:nYears){
        bYear[yy] ~ dnorm(0, tau.year)
      }
      # Abundance site and detection site-by-day random effects
      tau.site <- 1 / (sd.site * sd.site)
      sd.site ~ dunif(0, 3)
      tau.year <- 1 / (sd.year * sd.year)
      sd.year ~ dunif(0, 3)
     
      p <- 1 / (1 + exp(-lp)) 
      lp ~ dnorm(0, 0.1) # random delta defined implicitly
      
      # Likelihood
      # Ecological model for true abundance
      for (i in 1:n){       # Loop over R sites (95)
      
        log(lambda[i]) <- bYear[years[i]] + bSite[sites[i]]
        N[i] ~ dpois(lambda[i])               # Abundance
        y[i] ~ dbin(p, N[i])      # Detection
      
      }
      
    }",fill = TRUE)
sink()

# Initial values
N = sev_dat$y
N[N==0]=1
inits <- function(){list(N = N, lp = runif(1, -4, 4),
                         bYear = runif(22, -3, 3),  bSite = runif(5, -3, 3), 
                         sd.year = runif(1, 0, 1), sd.site = runif(1, 0, 1))}

# Parameters monitored
params <- c("bYear","bSite","sd.site","sd.year","p","N")

# MCMC settings
ni <- 100000
nt <- 50
nb <- 10000
nc <- 3

# Call WinBUGS from R (BRT 215 min)
library(R2WinBUGS)
out2 <- bugs(sev_dat, inits, params, "C:/Users/ac79/Mega/Projects/RICE/LTER/Analysis/SevBUGS_simple.txt", 
             n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE, working.directory = getwd())



#EXTRACT SPECIES DATA############################################################################
#only first 9 species in terms of abundance
sppList=names(table(GrasshopperSEV$SPECIES)[order(table(GrasshopperSEV$SPECIES),decreasing=T)])[1:9]

tiff("Results/sevilleta_Grasshoppers.tiff",unit="in",height=6.3,width=6.3,res=500,compression="lzw")

par(mfrow=c(3,3),mar=c(2.8,2.3,2,0.1),mgp=c(1.5,0.5,0))
for(i in 1:length(sppList)){
  
  oneSpp=subset(GrasshopperSEV,SPECIES==sppList[i])
  oneSpp=aggregate(CNT ~ YEAR + SITE,sum,data=oneSpp)
  
  oneSpp1=subset(oneSpp,YEAR!=2013) ; names(oneSpp1)[3]="CNT1"
  oneSpp2=subset(oneSpp,YEAR!=1992) ; names(oneSpp2)[3]="CNT2"
  oneSpp1$YEAR=oneSpp1$YEAR+1
  
  gr=merge(oneSpp1,oneSpp2)
  gr$logGr=log(gr$CNT2 / gr$CNT1)
  
  #climate variables
  gr=merge(gr,pjClimate)
  
  #Model selection
  models=list()
  models[[1]]=lm(logGr ~ sprTemp + sprTemp_2,data=gr)
  models[[2]]=lm(logGr ~ summTemp + summTemp_2,data=gr)
  models[[3]]=lm(logGr ~ tempTm1 + tempTm1_2,data=gr)
  models[[4]]=lm(logGr ~ sprPrecip + sprPrecip_2,data=gr)
  models[[5]]=lm(logGr ~ summPrecip + summPrecip_2,data=gr)
  models[[6]]=lm(logGr ~ precipTm1 + precipTm1_2,data=gr)
  
  models[[7]]=lm(logGr ~ sprTemp,data=gr)
  models[[8]]=lm(logGr ~ summTemp,data=gr)
  models[[9]]=lm(logGr ~ tempTm1,data=gr)
  models[[10]]=lm(logGr ~ sprPrecip,data=gr)
  models[[11]]=lm(logGr ~ summPrecip,data=gr)
  models[[12]]=lm(logGr ~ precipTm1,data=gr)
  
  names(models)=c("sprTemp_2","summTemp_2","tempTm1_2",
                  "sprPrecip_2","summPrecip_2","precipTm1_2",
                  "sprTemp","summTemp","tempTm1",
                  "sprPrecip","summPrecip","precipTm1")
  
  bestVariab=min(grep(AICtable(models)$model[1],names(gr)))
  bestVariab=gsub("_2","",names(gr)[bestVariab])
  xI=min(grep(bestVariab,names(gr)))
  
  plot(gr$logGr ~ gr[,xI],pch=16,main=sppList[i],
       xlab=names(gr)[xI],ylab="log_growth rate")
  xSeq=seq(min(gr[,xI]),max(gr[,xI]),by=0.1)
  
  bestModN=AICtable(models)$model_n[1]
  betas=coef(models[[bestModN]])
  
  if(bestModN > 6){ #best model is linear
    yMean=betas[1] + betas[2]*xSeq
  } else {  #best model is quadratic
    yMean=betas[1] + betas[2]*xSeq + betas[3]*(xSeq^2)
  }
  lines(xSeq,yMean) 
  
}

dev.off()





##################################################################################################################################################################
###OTHER BAYESIAN CODE#################################################################################
##################################################################################################################################################################
sink("Analysis/SevBUGS.txt")
cat("
    model{
    
    # Priors
    for(ss in 1:nSites){
    bSite[ss]  ~ dnorm(0, tau.site)
    lpSite[ss] ~ dnorm(0, tau.lp)
    pSite[ss]  <- 1 / (1 + exp(-lpSite[ss])) 
    }
    for(yy in 1:nYears){
    bYear[yy] ~ dnorm(0, tau.year)
    }
    # Abundance site and detection site-by-day random effects
    tau.site <- 1 / (sd.site * sd.site)
    sd.site ~ dunif(0, 3)
    tau.year <- 1 / (sd.year * sd.year)
    sd.year ~ dunif(0, 3)
    tau.lp <- 1 / (sd.lp * sd.lp)
    sd.lp ~ dunif(0, 3)
    
    #p <- 1 / (1 + exp(-lp)) 
    #lp ~ dnorm(0, 0.1) # random delta defined implicitly
    
    # Likelihood
    # Ecological model for true abundance
    for (i in 1:n){       # Loop over R sites (95)
    
    log(lambda[i]) <- bYear[years[i]] + bSite[sites[i]]
    N[i] ~ dpois(lambda[i])               # Abundance
    y[i] ~ dbin(pSite[sites[i]], N[i])      # Detection
    
    }
    
    }",fill = TRUE)
sink()

# Initial values
N = sev_dat$y
N[N==0]=1
inits <- function(){list(N = N, 
                         bSite = runif(5, -3, 3), bYear = runif(22, -3, 3), lpSite = runif(5, -3, 3), 
                         sd.year = runif(1, 0, 1), sd.site = runif(1, 0, 1), sd.lp = runif(1, 0, 1))}

# Parameters monitored
params <- c("bYear","bSite","pSite","lpSite","sd.site","sd.year","sd.lp","N")

# MCMC settings
ni <- 35000
nt <- 30
nb <- 5000
nc <- 3

# Call WinBUGS from R (BRT 215 min)
library(R2WinBUGS)
out2 <- bugs(sev_dat, inits, params, "C:/Users/ac79/Mega/Projects/RICE/LTER/Analysis/SevBUGS.txt", n.chains = nc, 
             n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE, working.directory = getwd())


sink("Analysis/sevGrass.stan")
cat("
    data {
        
      int<lower=0> n;             
      int<lower=0> nYears;       // N. of years.
      int<lower=0> nSites;       // N. of sites.
  
      int<lower=0> years[n];     // Index for years.
      int<lower=0> sites[n];     // Index for years.
  
      int y[n];                  // log size at time t+1 
    
    }
    
    parameters {

      real uY;                // intercept mean
      real<lower=0> sdY;      // intercept standard deviation
      real uS;                // intercept mean
      real<lower=0> sdS;      // intercept standard deviation

      vector[nYears] bYear;          
      vector[nSites] bSite;             
      real<lower=0,upper=1> p;   // Scaling factor. Maximum potential probability

    }

    #transformed parameters {
    #  int<lower=0> N[n];  
    #}

    model {
  
      int indY;
      int indS;
      int N[n];
      vector[n] m; // place holder for predictor growth
  
      ## Hyperpriors
      uY ~ normal(0,100);
      sdY ~ inv_gamma(0.001, 0.001);
      uS ~ normal(0,100);
      sdS ~ inv_gamma(0.001, 0.001);
  
      ## Priors
      # Year effects
      for (nI in 1:nYears){ 
        bYear[nI] ~ normal(uY, sdY);
      }
      # Site effects
      for (nI in 1:nSites){ 
        bSite[nI] ~ normal(uS, sdS);
      }
      p ~ uniform(0,1);
      
      ## MODEL
      for(nI in 1:n){
        indY  <-years[nI];
        indS  <-sites[nI];
        m[nI] <- bYear[indY] + bSite[indS];
        N[nI] ~ poisson_log(m[nI]);
      }
      y ~ binomial_logit(N,p);

    }",fill=T)
sink()

fit <- stan(file = 'Analysis/sevGrass.stan', data = sev_dat, 
             iter = 5000, warmup = 1000, chains = 4)





