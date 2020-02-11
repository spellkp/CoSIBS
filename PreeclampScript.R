#Libraries needed
#setwd("D:/CoSIBS/Project")
library(reshape2)
library(data.table)
library(randomForest)
library(rfPermute)
library(spatstat)
library("lattice")
library(ggplot2)
#Read in the data
data<-read.csv("CoSIBS preterm proteomics.csv", header = T, sep = ",")

#Subset for if used in the analysis (there are 102 obs) 2 have missing PH at day 7 so n = 100
data1=subset(data,data$cross_sec_day7==1)

#Subset Data into two useful pieces of information, clinical and protein
clinical.data=data1[,10:60]
protein.data=cbind(data1[,1:9],data1$SampID)
colnames(protein.data)=c(colnames(data1)[1:9],"SampID")

#Pick One Protein label ( I chose Target) and one Measurement variable (I chose lmeas, could also use meas)
Measurement.data=as.data.frame(cbind(as.character(protein.data$Target),protein.data$lmeas,protein.data$SampID))
colnames(Measurement.data)=c("Target","lmeas","SampID")
Measurement.data$lmeas=as.numeric(as.character(Measurement.data$lmeas))

#Create a wide data frame
wide.data=reshape(Measurement.data, idvar ="SampID", timevar = "Target", direction = "wide")

#Subset the clinical data to only contain unique IDs
clinic_unique=clinical.data[!duplicated(clinical.data$SampID),]
clinic_unique$SampID=factor(clinic_unique$SampID)

#Merge the two together; I used the data table package you could also use merge in base R
clinic.dt=as.data.table(clinic_unique)
wide.dt=as.data.table(wide.data)
setkey(clinic.dt,SampID)
setkey(wide.dt,SampID)
protein.clinical=clinic.dt[wide.dt]
protein.clinical=as.data.frame(protein.clinical)

#subset data to exclude missing values in ph7crit3
ph7=protein.clinical[is.na(protein.clinical$ph7crit3)==FALSE,]
predicts=cbind(ph7$ph7crit3,ph7[,52:1172])
predicts [ ,1:10]

#Replace random characters in column names with periods so that random forest will run
colnames(predicts)=c("ph7crit3",colnames(ph7)[52:1172])
colnames(predicts)=gsub("-",".",colnames(predicts))
colnames(predicts)=gsub(" ","",colnames(predicts))
colnames(predicts)=gsub(":",".",colnames(predicts))
colnames(predicts)=gsub("/",".",colnames(predicts))
colnames(predicts)=gsub(",",".",colnames(predicts))
colnames(predicts)[252]="lmeas.Ubiquitin.1"

#Setting seeds for random forest and running the forest

set.seed(52)
random=randomForest(formula=factor(predicts$ph7crit3)~.,data=predicts,importance=T)

#Heat map of all proteins
HM=impHeatmap(random)

mostimp=HM$data[HM$data[3]<10,]

# Heat maps of proteins in groups of 10
for(i in 0:10){
  
  mostimp=HM$data[HM$data[,3]>i*10 & HM$data[,3]<(i+1)*10,]
  
  ranks=T
  g <- ggplot(mostimp, aes_string("class", "predictor")) + geom_raster(aes_string(fill = "value"))
  g <- g + if (ranks) {
    scale_fill_gradient2("Rank", low = "#a50026", mid = "#ffffbf", 
                         high = "#313695", midpoint = mean(range(mostimp$value)), 
                         guide = guide_colorbar(reverse = TRUE))
  }
  
  g <- g + if (is.null(xlab)) 
    theme(axis.title.x = element_blank())
  
  g <- g + if (is.null(ylab)) 
    theme(axis.title.y = element_blank())
  
  plot(g)
}

set.seed(52)
library(AUCRF)
y= factor(predicts[,1]-1)
predicts2=cbind(y,predicts[,2:1122])  #0 is yes, 1 is no


fit=AUCRF(y~.,data=predicts2, trees=5000)
summary(fit)


random2=randomForest(formula=factor(predicts$ph7crit3)~ lmeas.Survivin + lmeas.Endocan +
                       lmeas.STK16 + lmeas.SCGF.alpha + lmeas.SCGF.beta + lmeas.BoneproteoglycanII + lmeas.GP114 +
                       lmeas.MMP.7 + lmeas.MSP + lmeas.MSPR + lmeas.NCK1 + lmeas.ON + lmeas.PDE5A + lmeas.PDGF.AA +
                       lmeas.PDGF.BB + lmeas.PF.4 + lmeas.TIMP.3 + lmeas.amyloidprecursorprotein,data=predicts,importance=T)
varImpPlot(random2, main='Variable Importance Plot' )
partialPlot(random2,predicts, x.var = "lmeas.Survivin")
partialPlot(random2,predicts, x.var = "lmeas.STK16")
partialPlot(random2,predicts, "lmeas.Endocan")
partialPlot(random2,predicts, "lmeas.GP114")
partialPlot(random2,predicts, "lmeas.MSP")
partialPlot(random2,predicts, "lmeas.BoneproteoglycanII")
partialPlot(random2,predicts, "lmeas.NCK1")
partialPlot(random2,predicts, "lmeas.PF.4")

#figure for letter
par(mfrow=c(3,2))
partialPlot(random2,predicts, "lmeas.GP114",1, main=" ", xlab= "GP114", ylab="Logit Pr (PVD)")
partialPlot(random2,predicts, "lmeas.MSP",1, main=" ", xlab= "MSP", ylab="Logit Pr (PVD)")
partialPlot(random2,predicts, "lmeas.Endocan",1, main=" ", xlab= "Endocan", ylab="Logit Pr (PVD)")
partialPlot(random2,predicts, x.var = "lmeas.Survivin",1, main=" ", xlab= "Survivin", ylab="Logit Pr (PVD)")
partialPlot(random2,predicts, "lmeas.BoneproteoglycanII",1, main=" ", xlab= "Bone proteoglycan II", ylab="Logit Pr (PVD)")
partialPlot(random2,predicts, "lmeas.NCK1",1, main=" ", xlab= "NCK1", ylab="Logit Pr (PVD)")



ivh <- clinic.dt[clinic.dt$ivh==1]
count(ivh$bpd1)

######## Using clinic.dt to make table ########

## Gender counts ##
library(plyr)
none$gender <- as.factor(none$gender)
count(none, 'gender')

## Mean, Sd, and Range Gestational Age ##
mean(none$gestational_age)
sd(none$gestational_age)
range(none$gestational_age)

## Mean, Sd, and Range Birth Weight ##
mean(none$weight)
sd(none$weight)
range(none$weight)

## Mean, Sd, and Range Birth Weight Z-score ##
mean(none$zscore)
sd(none$zscore)
range(none$zscore)

## Ventilation ##
count(none$ventilation)

## Surfactant ##
count(none$treatmentsurfact)

## Gestational Diabetes ##
count(none$gesdiabetes)

## Mode of Delivery ##
count(clinic.dt$csection)

## Mean, Sd, and Range Maternal Age ##
mean(none$mom_age)
sd(none$mom_age)
range(none$mom_age)

## Smoker ##
count(none$smoking)

## CSection ##
count(none$csection)

## Corticosteroids ##
count(none$corticosteroids)

## PPROM ##
count(none$prematrupmembran)

## Multiple Gestation ##
count(none$multgest)

## Multiple Infections ##
allinf <- vector(mode = "numeric", length = 64)
multinf <- none[,c(34:37)]
n = 63

for(i in 1:n){
  rowsum <- sum(multinf[i,])
  ifelse(rowsum == 4, allinf[i] <- 4,
         ifelse(rowsum == 3, allinf[i] <- 3,
                ifelse(rowsum == 2, allinf[i] <- 2,
                       ifelse(rowsum == 1, allinf[i] <- 1,
                              allinf[i] <- 0))))
}
multinf <- cbind(multinf, allinf)
oneinf <- multinf[multinf$allinf == 1]
  # Number of infections #
  count(allinf)

# Multiple Count
count(oneinf$pneumonia)
count(oneinf$necro)
count(oneinf$sepsis)
count(oneinf$disoxygen)



## My data frame ##
## Preeclampsia ##
preec <- cbind(clinic.dt$preeclampsia,protein.clinical[,52:1172])
colnames(preec)[1] <- "preeclampsia"
preec <- preec[!(preec$preeclampsia=="99"),]

colnames(preec)=gsub("-",".",colnames(preec))
colnames(preec)=gsub(" ","",colnames(preec))
colnames(preec)=gsub(":",".",colnames(preec))
colnames(preec)=gsub("/",".",colnames(preec))
colnames(preec)=gsub(",",".",colnames(preec))
colnames(preec)[252]="lmeas.Ubiquitin.1"

set.seed(49)
random=randomForest(formula=factor(preec$preeclampsia)~.,data=preec,importance=T)

HM=impHeatmap(random)

mostimp=HM$data[HM$data[3]<10,]

# Heat maps of proteins in groups of 10
for(i in 0:10){
  
  mostimp=HM$data[HM$data[,3]>i*10 & HM$data[,3]<(i+1)*10,]
  
  ranks=T
  g <- ggplot(mostimp, aes_string("class", "predictor")) + geom_raster(aes_string(fill = "value"))
  g <- g + if (ranks) {
    scale_fill_gradient2("Rank", low = "#a50026", mid = "#ffffbf", 
                         high = "#313695", midpoint = mean(range(mostimp$value)), 
                         guide = guide_colorbar(reverse = TRUE))
  }
  
  g <- g + if (is.null(xlab)) 
    theme(axis.title.x = element_blank())
  
  g <- g + if (is.null(ylab)) 
    theme(axis.title.y = element_blank())
  
  plot(g)
}

set.seed(52)
library(AUCRF)
y = factor(preec[,1]-1)
predicts2=cbind(y,preec[,2:1122])  #0 is yes, 1 is no


fit=AUCRF(y~.,data=predicts2, trees=5000)
summary(fit)

set.seed(52)
random2=randomForest(formula=factor(preec$preeclampsia)~ lmeas.sTie.2 + lmeas.GRB2.relatedadapterprotein2 +
                     lmeas.LAG.3 + lmeas.Dynactinsubunit2 + lmeas.Ficolin.3 + lmeas.Epithelialcellkinase + lmeas.MASP3 +
                     lmeas.FER + lmeas.FABP + lmeas.Cystatin.S + lmeas.DMP1 + lmeas.IL.11RA + lmeas.GRN + lmeas.Noggin,
                     data=preec,importance=T)
varImpPlot(random2, main='Variable Importance Plot' )


par(mfrow=c(3,2))
partialPlot(random2,preec, "lmeas.sTie.2",1, main=" ", xlab= "sTie.2", ylab="Logit Pr (Preeclamp)")
partialPlot(random2,preec, "lmeas.GRB2.relatedadapterprotein2",1, main=" ", xlab= "GRB2.relatedadapterprotein2", ylab="Logit Pr (Preeclamp)")
partialPlot(random2,preec, "lmeas.Noggin",1, main=" ", xlab= "Noggin", ylab="Logit Pr (Preeclamp)")
partialPlot(random2,preec, x.var = "lmeas.Cystatin.S",1, main=" ", xlab= "Cystatin.S", ylab="Logit Pr (Preeclamp)")
partialPlot(random2,preec, "lmeas.FER",1, main=" ", xlab= "FER", ylab="Logit Pr (Preeclamp)")
partialPlot(random2,preec, "lmeas.Ficolin.3",1, main=" ", xlab= "Ficolin.3", ylab="Logit Pr (Preeclamp)")



## Boxplots for Preeclamp

#sTie.2
par(mfrow=c(1,1))
cols <- rainbow(2, s = 0.5)
boxplot(lmeas.sTie.2 ~ preeclampsia, data = preec,
        at = c(1:2), col = cols, xaxs = FALSE, names = c("Yes", "No"), 
        xlab = "Preeclampsia", ylab = "lmeas")
#legend("topleft", fill = cols, legend = c("Yes","No"), horiz = T)
title(main= "sTie.2")

#Cystatin.S
par(mfrow=c(1,1))
cols <- rainbow(2, s = 0.5)
boxplot(lmeas.Cystatin.S ~ preeclampsia, data = preec,
        at = c(1:2), col = cols, xaxs = FALSE, names = c("Yes", "No"), 
        xlab = "Preeclampsia", ylab = "lmeas")
#legend("topleft", fill = cols, legend = c("Yes","No"), horiz = T)
title(main= "Cystatin.S")




