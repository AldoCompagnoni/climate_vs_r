setwd("C:/Users/Aldo/MEGA/Projects/RICE/LTER/")

#Identify results files
files=list.files("Results1")[grep("result",list.files("Results1"))]
#Put it together
metaResults=NULL #
for(i in 1:length(files)){
  metaResults=rbind(metaResults,read.csv(paste0("Results1/",files[i])))
}

#ANALYSIS#######################################################################

metaResults=metaResults[-which(metaResults$taxa=="-99999"),]

#Overall counts
metaResults$quadraticSign=metaResults$quadratic
metaResults$quadraticSign[metaResults$curvature>0]="quadPos"
metaResults$quadraticSign[metaResults$curvature<0]="quadNeg"
categories=summarise(group_by(metaResults,quadraticSign),count=n())
categories$proportions=categories$count/501  


#Ectotherm/Endotherms
metaResults$ectotherms="ectotherm"
metaResults$ectotherms[metaResults$taxa=="mammal" | metaResults$taxa=="bird"]="endotherm"
ectoEndo=summarise(group_by(metaResults,ectotherms,quadraticSign),count=n())
ectoEndo$proportion=ectoEndo$count/c(rep(sum(ectoEndo$count[1:4]),4),rep(sum(ectoEndo$count[5:8]),4))

head(ectoEndo)
tmp=as.table(rbind(ectoEndo$count[1:4],
               ectoEndo$count[5:8]))
dimnames(tmp)=list(organism=c("ecto","endo"),
                   relationship=c("Linear","NULL","quadNeg","quadPos"))
(Xsq <- chisq.test(tmp))


#Plants
plantProp=summarise(group_by(subset(metaResults,taxa=="plant"),quadraticSign),count=n())
plantProp$proportions=plantProp$count/sum(plantProp$count)

#Plants
metaResults$plankton="plankton"
metaResults$plankton[which(metaResults$taxa!="zooplankton" & metaResults$taxa!="plankton")]="nonPlankton"

#Plankton
planktonCount=summarise(group_by(metaResults,plankton,quadraticSign),count=n())
planktonCount$proportion=planktonCount$count/c(rep(sum(planktonCount$count[1:4]),4),rep(sum(planktonCount$count[5:8]),4))

round(planktonCount$proportion[5:8],2)*100

planktonCount$proportion=

  
tmp=as.table(rbind(planktonCount$count[1:4],
                   planktonCount$count[5:8]))
dimnames(tmp)=list(organism=c("nonPla","Pla"),
                   relationship=c("Linear","NULL","quadNeg","quadPos"))
(Xsq <- chisq.test(tmp))


par(mfrow=c(1,2))
boxplot(sensM ~ plankton,data=metaResults)
boxplot(curvature ~ plankton,data=metaResults)



#Proportions
metaResults$environment="terrestrial"
metaResults$environment[metaResults$LTER=="NorthTemperateLakes"]="fresh water"
metaResults$environment[metaResults$LTER=="SantaBarbara"]="marine"
metaResults$environment[metaResults$LTER=="Palmer" & metaResults$taxa=="zooplankton"]="marine"
metaResults$environment[metaResults$LTER=="Konza" & metaResults$taxa=="fish"]="fresh water"

remove=which(abs(metaResults$sens)>sd(abs(metaResults$sens),na.rm=T))

boxplot(abs(sens) ~ environment,data=metaResults[-remove,])
summary(lm(abs(sens) ~ environment,data=metaResults[-remove,]))

par(mfrow=c(1,1))
plot(metaResults$reps,metaResults$R2,pch=16)
mod=lm(metaResults$R2 ~ metaResults$repT)
abline(mod)
boxplot(R2 ~ quadratic,data=metaResults)
plot(metaResults$bestAICweight,metaResults$R2)

boxplot(abs(sens) ~ ectotherms,data=metaResults[-remove,])
summary(lm(abs(sens) ~ ectotherms,data=metaResults[-remove,]))



hist(metaResults$sens)



metaResults

summary(lm(metaResults$repT ~ quadratic,data=metaResults))

boxplot(quadratic ~ taxa,data=summarise(
  group_by(subset(metaResults,taxa=="bird" | taxa=="mammal"),taxa,quadratic),count=n()))


unique(metaResults$taxa)


metaResults$taxa[metaResults$taxa=="birds"]="bird"
metaResults$taxa[metaResults$taxa=="smallMammals"]="mammal"
metaResults$taxa[metaResults$taxa=="small mammal"]="mammal"
metaResults$taxa[metaResults$taxa=="small mammal"]="mammal"


unique(metaResults$taxa)
