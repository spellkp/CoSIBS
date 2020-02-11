#Libraries needed
#setwd("~/Desktop/CoSIBS/Project")
library(reshape2)
library(data.table)
library(randomForest)
library("spatstat")
library("lattice")
library(ggplot2)
library(rfPermute)
#library(descr)
#library(Hmisc)

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
#ph7=protein.clinical[is.na(protein.clinical$ph7crit3)==FALSE,]
#predicts=cbind(ph7$ph7crit3,ph7[,52:1172])
#predicts [ ,1:10]

#Replace random characters in column names with periods so that random forest will run
#colnames(predicts)=c("ph7crit3",colnames(ph7)[52:1172])
#colnames(predicts)=gsub("-",".",colnames(predicts))
#colnames(predicts)=gsub(" ","",colnames(predicts))
#colnames(predicts)=gsub(":",".",colnames(predicts))
#colnames(predicts)=gsub("/",".",colnames(predicts))
#colnames(predicts)=gsub(",",".",colnames(predicts))
#colnames(predicts)[252]="lmeas.Ubiquitin.1"

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
set.seed(52)
random2complications=randomForest(formula=factor(complications$`clinic.forest$complication`)~ lmeas.Ficolin.3+lmeas.SET+lmeas.Cystatin.S
                                    +lmeas.sTie.2+lmeas.HeparincofactorII+lmeas.GRB2.relatedadapterprotein2+ lmeas.NPS.PLA2
                                    +lmeas.DLL4+lmeas.Cytochromec+lmeas.Albumin
                                    +lmeas.GDF2+lmeas.HGFA+lmeas.CK.MB+lmeas.GranzymeB+lmeas.IL.1RAcP+lmeas.BNP.32
                                    +lmeas.CD27+lmeas.NCK1+lmeas.Glypican3+lmeas.GFRa.3               
                                    +lmeas.MP2K2+lmeas.TSP4+lmeas.Dkk.4+lmeas.vWF+lmeas.HRG+lmeas.HCG
                                    +lmeas.HSP70protein8+lmeas.a2.Antiplasmin+lmeas.C1r +lmeas.MASP3 + lmeas.FER
                                    + lmeas.CadherinE + lmeas.ENTP3 + lmeas.DKK1 + lmeas.PDK1 + lmeas.BoneproteoglycanII
                                    + lmeas.ADAMTS.5,data=complications,importance=T)

varImpPlot(random2complications, main='Variable Importance Plot' )

#var accuracy partial plot (first 6 proteins)
par(mfrow=c(3,2))
partialPlot(random2complications, complications, x.var = "lmeas.HRG")
partialPlot(random2complications, complications, x.var = "lmeas.IL.1RAcP")
partialPlot(random2complications, complications, "lmeas.sTie.2")
partialPlot(random2complications, complications, "lmeas.CK.MB")
partialPlot(random2complications, complications, "lmeas.GRB2.relatedadapterprotein2")
partialPlot(random2complications, complications, "lmeas.vWF")


#figure for letter (var accuracy)
par(mfrow=c(3,2))
partialPlot(random2complications,complications, "lmeas.HRG",1, main=" ", xlab= "HRG", ylab="Logit Pr (All)")
partialPlot(random2complications,complications, "lmeas.IL.1RAcP",1, main=" ", xlab= "IL.1RAcP", ylab="Logit Pr (All)")
partialPlot(random2complications,complications, "lmeas.sTie.2",1, main=" ", xlab= "sTie.2", ylab="Logit Pr (All)")
partialPlot(random2complications,complications, x.var = "lmeas.CK.MB",1, main=" ", xlab= "CK.MB", ylab="Logit Pr (All)")
partialPlot(random2complications,complications, "lmeas.GRB2.relatedadapterprotein2",1, main=" ", xlab= "GRB2.relatedadapterprotein2", ylab="Logit Pr (All)")
partialPlot(random2complications,complications, "lmeas.vWF",1, main=" ", xlab= "vWF", ylab="Logit Pr (All)")


## Boxplot with All Maternal Complications

#sTie.2
par(mfrow=c(1,1))
cols <- rainbow(3, s = 0.5)
boxplot(lmeas.sTie.2 ~ `clinic.forest$complication`, data = complications,
        at = c(1:3), col = cols, xaxs = FALSE, names = c("Other", "Preeclampsia", "Chorioamnionitis"), 
        xlab = "Maternal Complications", ylab = "lmeas")
#legend("topleft", fill = cols, legend = c("Yes","No"), horiz = T)
title(main= "sTie.2")

#GRB2.relatedadapterprotein2
par(mfrow=c(1,1))
cols <- rainbow(3, s = 0.5)
boxplot(lmeas.GRB2.relatedadapterprotein2 ~ `clinic.forest$complication`, data = complications,
        at = c(1:3), col = cols, xaxs = FALSE, names = c("Other", "Preeclampsia", "Chorioamnionitis"), 
        xlab = "Maternal Complications", ylab = "lmeas")
#legend("topleft", fill = cols, legend = c("Yes","No"), horiz = T)
title(main= "GRB2.relatedadapterprotein2")


#HRG
par(mfrow=c(1,1))
cols <- rainbow(3, s = 0.5)
boxplot(lmeas.HRG ~ `clinic.forest$complication`, data = complications,
        at = c(1:3), col = cols, xaxs = FALSE, names = c("Other", "Preeclampsia", "Chorioamnionitis"), 
        xlab = "Maternal Complications", ylab = "lmeas")
#legend("topleft", fill = cols, legend = c("Yes","No"), horiz = T)
title(main= "HRG")
############################################################################################################
#HRG THIS IS THE MONEY$$$$$$$$$$$$$$$$$$
par(mfrow=c(1,1))
cols <- rainbow(3, s = 0.5)
boxplot(lmeas.HRG ~ `clinic.forest$complication`, data = complications,
        at = c(1:3), col = cols, xaxs = FALSE, names = c("Other", "Preeclampsia", "Chorioamnionitis"), 
        xlab = "Maternal Complications", ylab = "lmeas")
#legend("topleft", fill = cols, legend = c("Yes","No"), horiz = T)
title(main= "HRG")
############################################################################################################
#### Sample Tree from Random Forest  #########
library("party")
x <- ctree(lmeas.HRG ~ ., data=complications)
plot(x, type="simple")

#########################################
options(repos='http://cran.rstudio.org')
have.packages <- installed.packages()
cran.packages <- c('devtools','plotrix','randomForest','tree','fifer')
to.install <- setdiff(cran.packages, have.packages[,1])
if(length(to.install)>0) install.packages(to.install)
library(devtools)
if(!('reprtree' %in% installed.packages())){
  install_github('araastat/reprtree')
}
for(p in c(cran.packages, 'reprtree')) eval(substitute(library(pkg), list(pkg=p)))


library(reprtree)
library(randomForest)
require(fifer)

library("randomForest", lib.loc="~/R/win-library/3.4")
rf_train <- random2complications

tr_struct = getTree(rf_train, 5, labelVar=TRUE)
reprtree:::plot.getTree(rf_train,k = 5)



#for(i in 496:500){
#  tr_struct = getTree(rf_train, i, labelVar=TRUE)
#  print(tr_struct)
#  print(complications[1,])
#  reprtree:::plot.getTree(rf_train,k = i)
#}

