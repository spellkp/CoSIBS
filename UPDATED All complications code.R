#Libraries needed
setwd("~/Desktop/CoSIBS/Project")
library(reshape2)
library(data.table)
library(randomForest)
library("spatstat")
library("lattice")
library(ggplot2)
library(rfPermute)
library(descr)
library(Hmisc)

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

############

##Random Forest

clinic.forest<-rbind(clinic.dt)
protein.forest<-rbind(protein.clinical)

## New complications variable
#Complication
#1=preeclampsia
#2=chorioamionitis 
#0=none

clinic.forest$complication[clinic.forest$preeclampsia==1]<-1
clinic.forest$complication[clinic.forest$chorioamnionitis==1]<-2
clinic.forest$complication[clinic.forest$chorioamnionitis != 1 & clinic.forest$preeclampsia != 1]<-0

complications<- cbind(clinic.forest$complication, protein.forest[,52:1172])
complications <- na.omit(complications)

#Replace random characters in column names with periods so that random forest will run
colnames(complications)=gsub("-",".",colnames(complications))
colnames(complications)=gsub(" ","",colnames(complications))
colnames(complications)=gsub(":",".",colnames(complications))
colnames(complications)=gsub("/",".",colnames(complications))
colnames(complications)=gsub(",",".",colnames(complications))
colnames(complications)[252]="lmeas.Ubiquitin.1"


set.seed(52)
random=randomForest(formula=factor(complications$`clinic.forest$complication`)~.,data=complications,importance=T)

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
y= factor(complications[,1]-1) 
complications2=cbind(y,complications[,2:1122])


fit=AUCRF(y~.,data=complications2, trees=5000)
summary(fit)

# boxplot: differences across two group 
#Unipro# for each - google if they have been seen before
random2complications=randomForest(formula=factor(complications$`clinic.forest$complication`)~ lmeas.BAFF+lmeas.FABP+lmeas.BMP.7+lmeas.GP114+lmeas.SLAF6+lmeas.CYTD+ lmeas.Caspase.3+lmeas.GIB+lmeas.IL.1sRI+lmeas.BoneproteoglycanII+  
                                    lmeas.Endoglin+lmeas.CSF.1+lmeas.Secretin+lmeas.YES+lmeas.EphA1+lmeas.CRK+lmeas.Epithelialcellkinase+lmeas.ECM1+lmeas.CD109+lmeas.OX2G+                
                                    lmeas.Carbonicanhydrase9+lmeas.LIN7B+lmeas.SLIK5+lmeas.CTLA.4+lmeas.Cystatin.S+lmeas.Dtk+lmeas.VEGFsR3+lmeas.MAPKAPK3+lmeas.PARK7 ,data=complications,importance=T)

varImpPlot(random2complications, main='Variable Importance Plot' )

#var accuracy partial plot (first 6 proteins)
par(mfrow=c(3,2))
partialPlot(random2complications, complications, x.var = "lmeas.GIB")
partialPlot(random2complications, complications, x.var = "lmeas.ECM1")
partialPlot(random2complications, complications, x.var = "lmeas.Secretin")
partialPlot(random2complications, complications, x.var = "lmeas.SLAF6")
partialPlot(random2complications, complications, x.var = "lmeas.IL.1sRI")
partialPlot(random2complications, complications, x.var = "lmeas.MAPKAPK3")


#figure for letter (var accuracy)
par(mfrow=c(3,2))
partialPlot(random2complications,complications, "lmeas.GIB",1, main=" ", xlab= "GIB", ylab="Logit Pr (All)")
partialPlot(random2complications,complications, "lmeas.ECM1",1, main=" ", xlab= "ECM1", ylab="Logit Pr (All)")
partialPlot(random2complications,complications, "lmeas.Secretin",1, main=" ", xlab= "Secretin", ylab="Logit Pr (All)")
partialPlot(random2complications,complications, x.var = "lmeas.SLAF6",1, main=" ", xlab= "SLAF6", ylab="Logit Pr (All)")
partialPlot(random2complications,complications, "lmeas.IL.1sRI",1, main=" ", xlab= "IL.1sRI", ylab="Logit Pr (All)")
partialPlot(random2complications,complications, "lmeas.MAPKAPK3",1, main=" ", xlab= "MAPKAPK3", ylab="Logit Pr (All)")





