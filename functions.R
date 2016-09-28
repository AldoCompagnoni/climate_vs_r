library(lubridate) ; library(dplyr)
library(reshape2) ; library(moments)
library(nlme) ; library(tidyr)
options(stringsAsFactors = F)

#Make my own AIC table. (NO COMMENT on AICtab() function from bbmle)
AICtable=function(x){
  
  aicValues=lapply(x,AIC)
  
  out=data.frame(model=names(aicValues),model_n=c(1:length(x)),aic=unlist(aicValues),row.names=NULL)
  out=out[order(out$aic),]
  out$deltaAIC=out$aic-out$aic[1]
  out$relLik=exp(-0.5* out$deltaAIC)
  out$weights=out$relLik/sum(out$relLik)
  
  out$weights=round(out$weights,4)
  
  return(out)
  
}


#Calculate log growth rate out of a vector of abundances.
#NOTE: first element of the vector is an NA 
logGr=function(x){ 
  out=x[2:length(x)]/x[1:(length(x)-1)]
  out=c(NA,log(out))
  out[out==Inf]=NA
  out[out==-Inf]=NA
  out[is.nan(out)]=NA
  return(out) 
}


#formatMeteoDate 
formatMeteoDate=function(x,format,separator){
  
  out=as.data.frame(matrix(unlist(strsplit(as.character(x),separator)),length(x),3,byrow=T))
  names(out)=c(format)
  for(i in 1:3) {out[,i]=as.numeric(as.character(out[,i]))}
  return(out)
  
}


#Function tests whether series of ingeters is continuous
#I use it to understand if there is a gap in a series of years
contiguousY=function(x){
  minY=min(x)
  series=seq(minY,minY+length(x)-1,by=1)
  if(sum(series!=x)>0) return(F)
  if(sum(series==x)==length(x)) return(T)
}


# Select species with enough replication
enough_rep_sppList <- function(x){
  
  #Find species with at least 
  #i)   10 transitions (years)
  #ii)  2 sites at least
  #iii) 20 data points at least
  sppListRaw  <- unique(x$species)
  sppKeep     <- NULL
  for(i in 1:length(sppListRaw)){
    
    tmp       <- subset(x,species == sppListRaw[i])
    if(length(unique(tmp$site)) > 1 & length(unique(tmp$year)) > 9 & nrow(tmp) > 20) {
      sppKeep <- c(sppKeep,sppListRaw[i])
    } 
    
  }

  return(sppKeep)
  
}


# Formats demographic data
# Input data frame should have four columns:
# Columns: year, site, species, abund (the latter for abundance)
#'replication' is the minimum number of YEARS for which there is spp abundance data. 
formatDemographicSite=function(x,replication=10){
  
  #Summarize by biggest replicate
  summ    <- summarise(group_by(x,year,site,species),abund=sum(abund,na.rm=T))
  if(sum(summ$abund==0,na.rm=T)>0) summ=summ[-which(summ$abund==0),]
  # keep species with enough replication
  sppList <- enough_rep_sppList(summ)
  
  # Store growth rates in a new data frame 
  gr=NULL
  for(i in 1:length(sppList)){
    
    tmp=subset(summ,species==sppList[i])
    dataWide=spread(tmp,site,abund,fill=NA)
    # produce continuous sequence of years (if needed) 
    if(contiguousY(dataWide$year)==F) {
      time=data.frame(year=seq(min(dataWide$year),max(dataWide$year),1))
      dataWide=merge(time,dataWide,all.x=T)
    }
    
    # Calculate growth rates
    dataGr    <- dataWide
    for(k in 3:ncol(dataGr)) { 
      dataGr[,k]=logGr(data.frame(dataGr[,k])[,1])
    }
    sppGr <- gather(dataGr,site,logGr,-year,-species)
    
    # Store abundance at time Nt0
    dataAbund       <- dataWide
    dataAbund       <- gather(dataAbund,site,Nt0,-year,-species)
    dataAbund       <- mutate(dataAbund,logNt0 = log(Nt0))
    # Align years (sppGr$year refers to year t+1)
    dataAbund$year  <- dataAbund$year + 1
    
    # merge
    popGrowth       <- merge(dataAbund,sppGr)
    
    #include data set only if there are > 10 replicate years for growth rates
    grSampleIDs=which(!is.na(popGrowth$logGr))
    yrSample=unique(sppGr$year[grSampleIDs])
    if(length(yrSample)>replication) gr=rbind(gr,popGrowth)
    
  }
  if(!is.null(gr)) {gr=gr[-which(is.na(gr$logGr)),]}
  
  # re-check for species with enough replication
  sppList <- enough_rep_sppList(gr)
  out     <- gr[which(gr$species %in% sppList),]
  return(out) 
  
}




#performs analyses, plots results, and returns general results!!!
analysisRandom=function(pop,speciesL=NULL,meteoVars,
                            organism,measure,LTERsite){
  
  meteoVars2  <- paste0(meteoVars,"_2")
  pop$site    <- as.factor(pop$site)
  if(is.null(speciesL)) speciesL=unique(pop$species)
  
  #Summary files
  nSpp   =length(speciesL)
  results=data.frame(R2=rep(NA,nSpp),quadratic=rep(NA,nSpp),curvature=rep(NA,nSpp),
                     bestAICweight=rep(NA,nSpp),bestPredictor=rep(NA,nSpp),repT=rep(NA,nSpp),
                     repS=rep(1,nSpp),species=rep(NA,nSpp),reps=rep(NA,nSpp),LTER=rep(LTERsite,nSpp),
                     taxa=rep(NA,nSpp),skewPred=rep(NA,nSpp),measure=rep(measure,nSpp),
                     sens=rep(NA,nSpp),cov_Nt0Meteo=rep(NA,nSpp))
  if(length(organism) == 1) results$taxa=rep(organism,nSpp)
  if(length(organism) > 1) results$taxa=organism
  
  
  #Loop across species
  for(i in 1:length(speciesL)){
    
    grD=subset(pop,species == speciesL[i])
    meteo=grD[,c("year",meteoVars)]
    meteo2=grD[,c("year",meteoVars2)]
    
    mod=list()
    nMetVars=(ncol(meteo)-1)
    for(col in 1:nMetVars){ #linear models
      eval(parse(n=1,text=paste0("mod[[col]]=lme(logGr ~ ",
                                 meteoVars[col]," + logNt0,random = ~ 1 | site,data=grD)")))
    }
    for(col in 1:nMetVars){ #quadratic models
      eval(parse(n=1,text=paste0("mod[[col+nMetVars]]=lme(logGr ~ ",meteoVars[col]," + ",
                                 meteoVars2[col]," + logNt0, random = ~ 1 | site,data=grD)")))
    }
    # Only density dependence
    mod[[length(mod)+1]]=lme(logGr ~ logNt0, random = ~ 1 | year,data=grD)
    # NULL model
    mod[[length(mod)+1]]=lme(logGr ~ 1, random = ~ 1 | year,data=grD)
    
    #Model results
    names(mod)=c(names(meteo)[-1],names(meteo2)[-1],"DD","NULL")
    bestMod=AICtable(mod)$model_n[1]
    bestVariab=gsub("_2","",names(mod)[bestMod])
    xI=grep(paste0("\\b",bestVariab,"\\b"),names(meteo)) #"\\b" gets the exact match
    betas <- fixef(mod[[bestMod]])
    
    #Print if best model is just density dependent
    if(names(mod)[bestMod] == "DD"){
      plot(grD$logGr ~ grD$logNt0,pch=16,col="grey",main=speciesL[i],
           xlab=expression("log(N"[t]*")"),ylab=expression("log("*lambda*")"))
      xN    <- seq(min(grD$logNt0),max(grD$logNt0),length.out = 100)
      yPred <- betas[1] + betas[2] * xN
      lines(xN,yPred,lwd=2)
      sens=NA
    } else { #Otherwise, proceed as "usual"
      
      if(sum(xI)==0) { xI=2 }
      plot(grD$logGr ~ meteo[,xI],pch=16,main=speciesL[i],col=grD$site,
           xlab=names(meteo)[xI],ylab=expression("log("*lambda*")"))
      #Meteo sequence
      xSeq=seq(min(meteo[,xI]),max(meteo[,xI]),by=0.1)
      nSeq=seq(min(meteo[,xI]),max(meteo[,xI]),by=0.1)
      
      # Graphs (and sensitivity analysis)
      sens=NA
      mM=mean(meteo[,xI],na.rm=T)
      # best model is linear
      if(bestMod < (nMetVars+1)){ 
        yMean <- betas[1] + betas[2]*xSeq + betas[3]*mean(grD$logNt0)
        lines(xSeq,yMean,lwd=2) 
        #y1=betas[1] + betas[2]*(mM*1.01)
        #y2=betas[1] + betas[2]*(mM*0.99)
        #sens=y1-y2
      } 
      # best model is quadratic
      if(bestMod > nMetVars & bestMod < length(mod)){ 
        yMean <- betas[1] + betas[2]*xSeq + betas[3]*xSeq^2 + betas[4]*mean(grD$logNt0)
        lines(xSeq,yMean,lwd=2) 
        #y1=betas[1] + betas[2]*(mM*1.01) + betas[3]*(mM*1.01)^2
        #y2=betas[1] + betas[2]*(mM*0.99) + betas[3]*(mM*0.99)^2 
        #sens=y1-y2
      }
      
    }
    
    #Write up results
    #results$R2[i]=summary(mod[[bestMod]])$adj.r.squared
    if(bestMod == length(mod)) results$quadratic[i]="NULL"
    if(bestMod == (length(mod)-1) ) results$quadratic[i]="DD"
    if(bestMod < (nMetVars+1)) results$quadratic[i]="linear"
    if(bestMod > nMetVars & bestMod < (length(mod)-1) ) {
      results$quadratic[i]="quadratic" ; results$curvature[i]=fixef(mod[[bestMod]])[3]
    }
    results$species[i]=speciesL[i]
    results$bestAICweight[i]=AICtable(mod)$weights[1]
    results$bestPredictor[i]=bestVariab
    results$repT[i]=length(unique(grD$year))
    results$repS[i]=mean(table(grD$year))
    results$reps[i]=nrow(grD)
    results$sens[i]=sens
    #results$cov_Nt0Meteo[i]=
    if(results$bestPredictor[i]!="NULL" & results$bestPredictor[i]!="DD") {
      results$skewPred[i]=skewness(meteo[,results$bestPredictor[i]])
    }
    
  }
  
  return(results)
  
}

