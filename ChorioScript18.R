
#Libraries needed
#setwd("C:/Users/brandie/Documents/grants/NHLBI R25")
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
##ph7=protein.clinical[is.na(protein.clinical$ph7crit3)==FALSE,]
##predicts=cbind(ph7$ph7crit3,ph7[,52:1172])
##predicts [ ,1:10]

## Data Frame for Maternal Complication (Chorio)
chorio <- cbind(clinic.dt$chorioamnionitis, protein.clinical[,52:1172])
colnames(chorio)[1] <- "chorioamnionitis"
chorio <- chorio[!(chorio$chorioamnionitis == "99"),] #; chorio

#Replace random characters in column names with periods so that random forest will run
#colnames(chorio)=c("ph7crit3",colnames(ph7)[52:1172])
colnames(chorio)=gsub("-",".",colnames(chorio))
colnames(chorio)=gsub(" ","",colnames(chorio))
colnames(chorio)=gsub(":",".",colnames(chorio))
colnames(chorio)=gsub("/",".",colnames(chorio))
colnames(chorio)=gsub(",",".",colnames(chorio))
colnames(chorio)[252]="lmeas.Ubiquitin.1"

#Setting seeds for random forest and running the forest

set.seed(52)
random=randomForest(formula=factor(chorio$chorioamnionitis)~.,data=chorio,importance=T)

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
y= factor(chorio[,1]-1) # 0-yes 1-no
chorio2=cbind(y,chorio[,2:1122])


fit=AUCRF(y~.,data=chorio2, trees=5000)
summary(fit)

## Most Important Proteins
#colMax <- function(fit) sapply(fit, max, na.rm = TRUE)
#colSort <- function(fit, ...) sapply(fit, sort, ...)


# boxplot: differences across two group 
#Unipro# for each - google if they have been seen before
set.seed(52)
random2chorio=randomForest(formula=factor(chorio$chorioamnionitis)~ lmeas.IL.1sRI + lmeas.HRG +
                       lmeas.GranzymeB + lmeas.M.CSFR + lmeas.BNP.32 + lmeas.TAFI + lmeas.YES +
                       lmeas.SMAC + lmeas.NET4 + lmeas.Thymidinekinase + lmeas.BGN + lmeas.IgE + lmeas.SAA + lmeas.ASAH2 +
                       lmeas.ATS13 + lmeas.MP2K2 + lmeas.SLAF6 + lmeas.RBP,data=chorio,importance=T)
varImpPlot(random2chorio, main='Variable Importance Plot' )
partialPlot(random2chorio, chorio, x.var = "lmeas.HRG")
partialPlot(random2chorio, chorio, x.var = "lmeas.IL.1sRI")
partialPlot(random2chorio, chorio, "lmeas.BGN")
partialPlot(random2chorio, chorio, "lmeas.RBP")
partialPlot(random2chorio, chorio, "lmeas.SLAF6")
partialPlot(random2chorio, chorio, "lmeas.SMAC")
#partialPlot(random2chorio, chorio, "lmeas.YES")
#partialPlot(random2chorio, chorio, "lmeas.SMAC")

#figure for letter
par(mfrow=c(3,2))
partialPlot(random2chorio, chorio, "lmeas.HRG",1, main=" ", xlab= "HRG", ylab="Logit Pr (Chorio)")
partialPlot(random2chorio, chorio, "lmeas.IL.1sRI",1, main=" ", xlab= "IL.1sRI", ylab="Logit Pr (Chorio)")
partialPlot(random2chorio, chorio, "lmeas.BGN",1, main=" ", xlab= "BGN", ylab="Logit Pr (Chorio)")
partialPlot(random2chorio, chorio, x.var = "lmeas.RBP",1, main=" ", xlab= "RBP", ylab="Logit Pr (Chorio)")
partialPlot(random2chorio, chorio, "lmeas.SLAF6",1, main=" ", xlab= "SLAF6", ylab="Logit Pr (Chorio)")
partialPlot(random2chorio, chorio, "lmeas.SMAC",1, main=" ", xlab= "SMAC", ylab="Logit Pr (Chorio)")


#treetime <- getTree(random2chorio, k = 1, labelVar = TRUE)
#plot.getTree(treetime)
#library(reprtree)
#reprtree:::plot.getTree(random2chorio)

## Boxplot with Chorio Maternal Complication

#SLAF6
par(mfrow=c(1,1))
cols <- rainbow(2, s = 0.5)
boxplot(lmeas.SLAF6 ~ chorioamnionitis, data = chorio,
        at = c(1:2), col = cols, xaxs = FALSE, names = c("Yes", "No"), 
        xlab = "Chorioamnionitis", ylab = "lmeas")
#legend("topleft", fill = cols, legend = c("Yes","No"), horiz = T)
title(main= "SLAF6")

#HRG
par(mfrow=c(1,1))
cols <- rainbow(2, s = 0.5)
boxplot(lmeas.HRG ~ chorioamnionitis, data = chorio,
        at = c(1:2), col = cols, xaxs = FALSE, names = c("Yes", "No"), 
        xlab = "Chorioamnionitis", ylab = "lmeas")
#legend("topleft", fill = cols, legend = c("Yes","No"), horiz = T)
title(main= "HRG")



## Data Frame for Maternal Complication 
## Chorio 
chorio <- cbind(clinic.dt$chorioamnionitis, protein.clinical[,52:1172])
colnames(chorio)[1] <- "chorio"
chorio <- chorio[!(chorio$chorio == "99"),]; chorio

########## New Dataframe for Chorio #########
chorio.newdt <- clinic.dt[clinic.dt$chorioamnionitis == 1]
View(chorio.newdt)

## Gender counts ##
library(plyr)
chorio.newdt$gender <- as.factor(chorio.newdt$gender)
count(chorio.newdt, 'gender')

## Maternal Complications ##
#count(chorio.newdt, 'preeclampsia')
#count(chorio.newdt, 'chorioamnionitis')
#count(clinic.dt, 'ivh')

## Mean, Sd, and Range Gestational Age ##
mean(chorio.newdt$gestational_age)
sd(chorio.newdt$gestational_age)
range(chorio.newdt$gestational_age)

## Mean, Sd, and Range Birth Weight ##
mean(chorio.newdt$weight)
sd(chorio.newdt$weight)
range(chorio.newdt$weight)

## Mean, Sd, and Range Birth Weight Z-score ##
mean(chorio.newdt$zscore)
sd(chorio.newdt$zscore)
range(chorio.newdt$zscore)

## Ventilation ##
count(chorio.newdt$ventilation)

## Surfactant ##
count(chorio.newdt$treatmentsurfact)

## Gestational Diabetes ##
count(chorio.newdt$gesdiabetes)

## Mode of Delivery ##
count(chorio.newdt$csection)

## Multiple Infections ##
allinf <- vector(mode = "numeric", length = 14)
multinf <- chorio.newdt[,c(33:36)]
n = length(multinf$pneumonia)

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

## Smoker ##
count(chorio.newdt$smoking)

## Maternal Age ##
mean(chorio.newdt$mom_age)
sd(chorio.newdt$mom_age)
range(chorio.newdt$mom_age)


count(chorio.newdt$corticosteroids)
count(chorio.newdt$prolongedrupmemb)
count(chorio.newdt$multgest)

########### New Dataframe for bpd ####################
bpd1.newdt <- clinic.dt[clinic.dt$bpd1 != 1]
#View(bpd1.newdt)

## Gender counts ##
library(plyr)
bpd1.newdt$gender <- as.factor(bpd1.newdt$gender)
count(bpd1.newdt, 'gender')

## Maternal Complications ##
#count(chorio.newdt, 'preeclampsia')
#count(chorio.newdt, 'chorioamnionitis')
#count(clinic.dt, 'ivh')

## Mean, Sd, and Range Gestational Age ##
mean(bpd1.newdt$gestational_age)
sd(bpd1.newdt$gestational_age)
range(bpd1.newdt$gestational_age)

## Mean, Sd, and Range Birth Weight ##
mean(bpd1.newdt$weight)
sd(bpd1.newdt$weight)
range(bpd1.newdt$weight)

## Mean, Sd, and Range Birth Weight Z-score ##
mean(bpd1.newdt$zscore)
sd(bpd1.newdt$zscore)
range(bpd1.newdt$zscore)

## Ventilation ##
count(bpd1.newdt$ventilation)

## Surfactant ##
count(bpd1.newdt$treatmentsurfact)

## Gestational Diabetes ##
count(bpd1.newdt$gesdiabetes)

## Mode of Delivery ##
count(bpd1.newdt$csection)

## Multiple Infections ##
allinf <- vector(mode = "numeric", length = 14)
multinf <- bpd1.newdt[,c(33:36)]
n = length(multinf$pneumonia)

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

## Smoker ##
count(bpd1.newdt$smoking)

## Maternal Age ##
mean(bpd1.newdt$mom_age)
sd(bpd1.newdt$mom_age)
range(bpd1.newdt$mom_age)

count(bpd1.newdt$corticosteroids)
count(bpd1.newdt$prolongedrupmemb)
count(bpd1.newdt$multgest)


