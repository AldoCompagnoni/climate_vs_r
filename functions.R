library(lubridate) ; library(dplyr)
library(reshape2) ; library(moments)
library(nlme)
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

#Function returns if a series of ingeters is continuous
#I use it to understand if there is a gap in a series of years
contiguousY=function(x){
  minY=min(x)
  series=seq(minY,minY+length(x)-1,by=1)
  if(sum(series!=x)>0) return(F)
  if(sum(series==x)==length(x)) return(T)
}


#This function formats demographic data
#Format of the file should include:
#Columns: year, site, species, abund (the latter for abundance)
#'replication' is the minimum number of YEARS for which there is spp abundance data. 
formatDemographicSite=function(x,replication=10){
  
  #Summarize by biggest replicate
  summ=summarise(group_by(x,year,site,species),abund=sum(abund,na.rm=T))
  if(sum(summ$abund==0,na.rm=T)>0) summ=summ[-which(summ$abund==0),]
  
  #Find species with at least 10 transitions enough replication
  allSpp=colnames(table(summ[,c("year","species")])) #Total spp list
  sppList=allSpp[which(apply(table(summ[,c("year","species")])>0,2,sum)>10)]
  
  gr=NULL
  for(i in 1:length(sppList)){
    
    tmp=subset(summ,species==sppList[i])
    dataWide=dcast(tmp,year ~ site,mean,na.rm=T,value.var="abund")
    if(contiguousY(dataWide$year)==F) { #produce continuous sequence of years (if needed) 
      time=data.frame(year=seq(min(dataWide$year),max(dataWide$year),1))
      dataWide=merge(time,dataWide,all.x=T)
    }
    for(k in 2:ncol(dataWide)) { dataWide[,k]=logGr(dataWide[,k]) }
    dataWide$species=sppList[i]
    sppGr=melt(dataWide,id.var=c("year","species")) 
    names(sppGr)[3:4]=c("site","logGr")
    #include data set only if there are > 10 replicate years for growth rates
    grSampleIDs=which(!is.na(sppGr$logGr))
    yrSample=unique(sppGr$year[grSampleIDs])
    if(length(yrSample)>replication) gr=rbind(gr,sppGr)
    
  }
  if(!is.null(gr)) {gr=gr[-which(is.na(gr$logGr)),]}
  
  return(gr) 
  
}


#performs analyses, plots results, and returns general results!!!
analysisRandom=function(pop,speciesL=NULL,meteoVars,
                        organism,measure,LTERsite){
  
  meteoVars2=paste0(meteoVars,"_2")
  if(is.null(speciesL)) speciesL=unique(pop$species)
  
  newPop=NULL
  for(i in 1:length(speciesL)) { 
    tmp=subset(pop,species==speciesL[i])
    if(length(unique(tmp$site))>2) newPop=rbind(newPop,tmp) 
  }

  #Summary files
  nSpp=length(unique(newPop$species))
  results=data.frame(R2=rep(NA,nSpp),quadratic=rep(NA,nSpp),curvature=rep(NA,nSpp),
                     bestAICweight=rep(NA,nSpp),bestPredictor=rep(NA,nSpp),repT=rep(NA,nSpp),
                     repS=rep(1,nSpp),species=rep(NA,nSpp),reps=rep(NA,nSpp),LTER=rep(LTERsite,nSpp),
                     taxa=rep(NA,nSpp),skewPred=rep(NA,nSpp),measure=rep(measure,nSpp),
                     sensM=rep(NA,nSpp),sensV=rep(0,nSpp))
  if(length(organism) == 1) results$taxa=rep(organism,nSpp)
  if(length(organism) > 1) results$taxa=organism[which(speciesL %in% unique(newPop$species))]
  
  #Loop across species
  for(i in 1:length(unique(newPop$species))){
    
    grD=subset(pop,species==unique(newPop$species)[i])

    meteo=grD[,c("year",meteoVars)]
    meteo2=grD[,c("year",meteoVars2)]
    
    mod=list()
    nMetVars=(ncol(meteo)-1)
    for(colID in 1:nMetVars){ #linear models
      eval(parse(n=1,text=paste0("mod[[colID]]=lme(logGr ~ ",meteoVars[colID],",random = ~ 1 | year,data=grD)")))
    }
    for(colID in 1:nMetVars){ #quadratic models
      eval(parse(n=1,text=paste0("mod[[colID+nMetVars]]=lme(logGr ~ ",meteoVars[colID]," + ",
                                 meteoVars2[colID],", random = ~ 1 | year,data=grD)")))
    }
    mod[[length(mod)+1]]=lme(logGr ~ 1, random = ~ 1 | year,data=grD)
    
    #Model results
    names(mod)=c(names(meteo)[-1],names(meteo2)[-1],"NULL")
    bestMod=AICtable(mod)$model_n[1]
    bestVariab=gsub("_2","",names(mod)[bestMod])
    xI=grep(paste0("\\b",bestVariab,"\\b"),names(meteo)) #"\\b" gets the exact match
    
    if(sum(xI)==0) { xI=2}
    colorz=as.numeric(grD$site)
    plot(grD$logGr ~ meteo[,xI],pch=16,main=unique(newPop$species)[i],
         col=colorz,xlab=names(meteo)[xI],ylab=expression("log("*lambda*")"))
    xSeq=seq(min(meteo[,xI]),max(meteo[,xI]),by=0.1)
    betas=as.numeric(fixef(mod[[bestMod]]))
    
    #Graphs (and sensitivity analysis)
    cM=mean(meteo[,xI],na.rm=T)
    cV=var(meteo[,xI],na.rm=T)
    sensMean=sensVariance=0 #for null models
    #best model is linear
    if(bestMod < (nMetVars+1)){ 
      yMean=betas[1] + betas[2]*xSeq ; lines(xSeq,yMean,lwd=2) 
      sensMean      <-  (cM*0.01) / cM
      sensVariance  <-  0 
    } 
    #best model is quadratic
    if(bestMod > nMetVars & bestMod < length(mod)){ 
      yMean=betas[1] + betas[2]*xSeq + betas[3]*xSeq^2 ; lines(xSeq,yMean,lwd=2) 
      sensMean        <- (((2*betas[3]*cM + betas[2])*cM) / (((betas[2]*cM)+(betas[3]*cM^2)) + (betas[3]*cV))) * ((cM*0.01) / cM)
      sensVariance    <- ( (betas[3]*cV) / (((betas[2]*cM)+(betas[3]*cM^2)) + (betas[3]*cV)) ) * ((cV*0.01) / cV)
    }
    
    #Write up results
    #results$R2[i]=summary(mod[[bestMod]])$adj.r.squared
    if(bestMod==length(mod)) results$quadratic[i]="NULL"
    if(bestMod < (nMetVars+1)) results$quadratic[i]="linear"
    if(bestMod > nMetVars & bestMod < length(mod)) {
      results$quadratic[i]="quadratic" ; results$curvature[i]=fixef(mod[[bestMod]])[3]
    }
    results$species[i]=speciesL[i]
    results$bestAICweight[i]=AICtable(mod)$weights[1]
    results$bestPredictor[i]=bestVariab
    results$repT[i]=length(unique(grD$year))
    results$repS[i]=mean(table(grD$year))
    results$reps[i]=nrow(grD)
    results$sensM[i]=sensMean
    results$sensV[i]=sensVariance
    if(results$bestPredictor[i]!="NULL") results$skewPred[i]=skewness(meteo[,results$bestPredictor[i]])
  }
  
  return(results)
  
}

