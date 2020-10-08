#### ULKAR IMAMVERDIYEVA AND PATRICK SHEA ###
### Female Leaders and Foreign Policy     ###
### Part 2: MULTIPLE TREATMENT ANALYSIS   ###
#############################################

rm(list=ls())
library(Synth)
library(foreign)
library(tictoc)


#########  MILITARY SPENDING 1960-2014 ###################
## Multiple Treatment Estimation ##


###1.  DATA clean up ###
rm(list=ls())
setwd("~/Dropbox/Women Leaders and Conflict/Empirical Analysis/Rdata and scripts")
library(readstata13)
mt_ms1 <- read.dta13("mt_ms 2695.dta")
mt_ms1 <- mt_ms1[mt_ms1$year<2014,]
mt_ms1<-mt_ms1[mt_ms1$cown!=732,]
colnames(mt_ms1)[4]=c("female")

# 1977 recoded as male=0, as Gandhi was in office till March 23#
mt_ms1$female[(mt_ms1$cown==750 & mt_ms1$year==1977)] <- 0
#Meir only was in office till June 03, so I have recoded 1974 as having a male leader
mt_ms1$female[(mt_ms1$cown==666 & mt_ms1$year==1974)] <- 0
# canada only one year treatment coded as non-treatment
mt_ms1$female[(mt_ms1$cown==20 & mt_ms1$year==1993)] <- 0
# Ecuador only one year treatment coded as non-treatment
mt_ms1$female[(mt_ms1$cown==130 & mt_ms1$year==1997)] <- 0
# delete switzerland
mt_ms1<-mt_ms1[mt_ms1$cown!=225,]
# norway one year treatment coded as zero , 1981
mt_ms1$female[(mt_ms1$cown==385 & mt_ms1$year==1981)] <- 0
mt_ms1$female[(mt_ms1$cown==385 & mt_ms1$year==2013)] <- 0
# India second treatment coded as zero #
mt_ms1$female[(mt_ms1$cown==750 & mt_ms1$year>1979 & mt_ms1$year<1985)] <- 0
# Sri lanka coded female, if the president is female 
mt_ms1$female[(mt_ms1$cown==780)] <- 0
mt_ms1$female[(mt_ms1$cown==780 & mt_ms1$year>1994 & mt_ms1$year<2005)] <- 1
# Pakistan second treatment coded as zero
mt_ms1$female[(mt_ms1$cown==770 & mt_ms1$year>1992 & mt_ms1$year<1997)] <- 0

### Duplicate countries with two treatments
# Original Argentina is cown 160, and a replica Argentina cown 161
df160 <- data.frame(matrix(ncol = 496, nrow = 54))
df160 <- mt_ms1[mt_ms1$cown==160,]
df160$cown=161
df160$countrycode_pw <- "Arg-new"
df160$female[(df160$year>1973 & df160$year<1977)] <- 0
mt_ms1$female[(mt_ms1$cown==160 & mt_ms1$year>2006)] <- 0

# Append df160
mt_ms1 <- rbind(mt_ms1, df160)

# Original Philippines is cown 840, and a replica Philippines is cown 841
df840 <- data.frame(matrix(ncol = dim(mt_ms1)[2], nrow = 54))
df840 <- mt_ms1[mt_ms1$cown==840,]
df840$cown=841
df840$countrycode_pw <- "Phil-new"
df840$female[(df840$year>1985 & df840$year<1993)] <- 0
mt_ms1$female[(mt_ms1$cown==840 & mt_ms1$year>2000)] <- 0

# Append df840
mt_ms1 <- rbind(mt_ms1, df840)

### 2. Data set-up for synth ###
#### Identify treated countries ###
macro <- aggregate(female ~ cown, data = mt_ms1, sum)
macro$donor <- 0
macro$treated <- 0
macro$donor[macro$female==0]=1
macro$treated[macro$female>0]=1
macro <- subset(macro, select=-c(female))

treated <- macro$cown[macro$treated==1]
donor <- macro$cown[macro$donor==1]

tr.start <- numeric(length(treated))
tr.end <- numeric(length(treated))

for (i in 1:length(treated))
{
    tr.years <- mt_ms1$year[(mt_ms1$female==1 & mt_ms1$cown==treated[i])]
    tr.start[i] <- min(tr.years)
    tr.end[i] <- max(tr.years)
}

mt_ms1 <- merge(mt_ms1, macro, by="cown")

# Covaraites
predictors=c("polity2_20", "logrgdp_f2", "rgdppc", "cinc_c", "logtrade", "logthreat_f2", "allies", "issues", "rival")

#### 3. Running Multiple Treatment SCM ###
tic()

for (i in 1:length(treated))
{
    dataprep.out<- 
        dataprep( 
            foo = mt_ms1, 
            predictors = predictors,
            predictors.op = "mean", 
            special.predictors = list(list("logsmilex_i",tr.start[i]-1, c("mean"))),
            dependent = "logsmilex_i",  
            unit.variable = "cown", 
            time.variable = "year",  
            treatment.identifier = treated[i], 
            controls.identifier = donor, 
            time.predictors.prior = c((max(1960,tr.start[i]-10)):tr.start[i]), 
            time.optimize.ssr = c((max(1960,tr.start[i]-10)):tr.start[i]), 
            unit.names.variable = "countrycode_pw", 
            time.plot = c((max(1960,tr.start[i]-10)):tr.end[i]) 
        )
    
    synth.out <- synth(dataprep.out)
    synth.tables <- synth.tab( dataprep.res = dataprep.out, synth.res = synth.out) 
    assign(paste("W",treated[i], sep=""), (synth.tables$tab.w))
    assign(paste("V",treated[i], sep=""), (synth.out$solution.v))
    assign(paste("Pred",treated[i], sep=""), (synth.tables$tab.pred))
    
    assign(paste("Y1",treated[i], sep="."), (dataprep.out$Y1))
    assign(paste("Y1s",treated[i], sep="."), (dataprep.out$Y0plot%*%synth.out$solution.w))
    assign(paste("Alpha",treated[i], sep="."), (dataprep.out$Y1-dataprep.out$Y0plot%*%synth.out$solution.w))
    
    ### Run placebo for each treated unit
    for (j in 1:length(donor))
    {
        
        dataprep.out<- 
            dataprep( 
                foo = mt_ms1, 
                predictors = predictors,
                special.predictors = list(list("logsmilex_i",tr.start[i]-1, c("mean"))),
                predictors.op = "mean", 
                dependent = "logsmilex_i",  
                unit.variable = "cown", 
                time.variable = "year",  
                treatment.identifier = donor[j], 
                controls.identifier = donor[-j], 
                time.predictors.prior = c((max(1960,tr.start[i]-10)):tr.start[i]), 
                time.optimize.ssr = c((max(1960,tr.start[i]-10)):tr.start[i]), 
                unit.names.variable = "countrycode_pw", 
                time.plot = c((max(1960,tr.start[i]-10)):tr.end[i]) 
            )
        
        synth.out <- synth(dataprep.out)
        synth.tables <- synth.tab( dataprep.res = dataprep.out, synth.res = synth.out) 
        assign(paste("Y1",treated[i],donor[j], sep="."), (dataprep.out$Y1))
        assign(paste("Y1s",treated[i],donor[j], sep="."), (dataprep.out$Y0plot%*%synth.out$solution.w))
        assign(paste("Alpha",treated[i],donor[j], sep="."), (dataprep.out$Y1-dataprep.out$Y0plot%*%synth.out$solution.w))
        
    } #j loop
    
    #}
} # i loop
toc()

### Merging back outcomes with centered at 0(treatment year)
for (i in 1: length(treated))
{
    temp <- get(paste("Alpha",treated[i], sep="."))
    year <- as.numeric(rownames(temp))
    year <- year-tr.start[i]+1
    
    assign(paste("alpha",treated[i], sep=""), (data.frame(date=year, temp)))
    temp3 <- get(paste("alpha",treated[i], sep=""))
    
    for (j in 1: length(donor))
    {
        temp2 <- get(paste("Alpha",treated[i],donor[j], sep="."))
        temp3 <- cbind(temp3, temp2)
        colnames(temp3)[j+2] <- paste("Alpha", treated[i],donor[j], sep=".")
    }
    assign(paste("alpha",treated[i], sep=""), (temp3))
}


for (i in 1: length(treated))
{
    temp <- get(paste("Y1",treated[i], sep="."))
    year <- as.numeric(rownames(temp))
    year <- year-tr.start[i]+1
    
    assign(paste("y1",treated[i], sep="."), (data.frame(date=year, temp)))
    temp3 <- get(paste("y1",treated[i], sep="."))
    
    for (j in 1: length(donor))
    {
        temp2 <- get(paste("Y1",treated[i],donor[j], sep="."))
        temp3 <- cbind(temp3, temp2)
        colnames(temp3)[j+2] <- paste("Y1", treated[i],donor[j], sep=".")
    }
    assign(paste("y1",treated[i], sep="."), (temp3))
}


for (i in 1: length(treated))
{
    temp <- get(paste("Y1s",treated[i], sep="."))
    year <- as.numeric(rownames(temp))
    year <- year-tr.start[i]+1
    
    assign(paste("y1s",treated[i], sep="."), (data.frame(date=year, temp)))
    temp3 <- get(paste("y1s",treated[i], sep="."))
    
    for (j in 1: length(donor))
    {
        temp2 <- get(paste("Y1s",treated[i],donor[j], sep="."))
        temp3 <- cbind(temp3, temp2)
        colnames(temp3)[j+2] <- paste("Y1s", treated[i],donor[j], sep=".")
    }
    colnames(temp3)[2] <- paste("w", treated[i], sep="")
    assign(paste("y1s",treated[i], sep="."), (temp3))
}

alpha <- merge(alpha94, alpha140, by="date", all=T)
alpha <- merge(alpha, alpha155, by="date", all=T)
alpha <- merge(alpha, alpha160, by="date", all=T)
alpha <- merge(alpha, alpha161, by="date", all=T)
alpha <- merge(alpha, alpha200, by="date", all=T)
alpha <- merge(alpha, alpha375, by="date", all=T)
alpha <- merge(alpha, alpha385, by="date", all=T)
alpha <- merge(alpha, alpha390, by="date", all=T)
alpha <- merge(alpha, alpha640, by="date", all=T)
alpha <- merge(alpha, alpha666, by="date", all=T)
alpha <- merge(alpha, alpha750, by="date", all=T)
alpha <- merge(alpha, alpha770, by="date", all=T)
alpha <- merge(alpha, alpha780, by="date", all=T)
alpha <- merge(alpha, alpha800, by="date", all=T)
alpha <- merge(alpha, alpha840, by="date", all=T)
alpha <- merge(alpha, alpha841, by="date", all=T)
alpha <- merge(alpha, alpha900, by="date", all=T)
alpha <- merge(alpha, alpha920, by="date", all=T)

coln.treated <- character(length(treated))
coln.donor <- matrix(0,length(donor),length(treated))

for (i in 1:length(treated))
{
    coln.treated[i] <-paste("X", treated[i], sep="") 
    
    for (j in 1:length(donor))
    {
        coln.donor[j,i] <-paste("Alpha", treated[i], donor[j], sep=".")
    }
}
coln.donor <- as.vector(coln.donor)

alpha.treated <- alpha[,colnames(alpha) %in% coln.treated]
alpha.treated$date <- alpha$date

alpha.donor <- alpha[,colnames(alpha) %in% coln.donor]
alpha.donor$date <- alpha$date

alpha.treated$ave.tr <- rowMeans(alpha.treated[,1:length(treated)], na.rm=TRUE)
alpha.donor$ave.d <- rowMeans(alpha.donor[,1:length(coln.donor)], na.rm=TRUE)

## 5. Pre RMSPE comparison ##
## RMSPE comparison ##
rmspe.tr <- numeric(length(treated))
rmspe.d <- matrix(0, length(donor),length(treated))


for (i in 1:length(treated))
{
    rmspe.tr[i] <- sqrt(mean(alpha.treated[1:10, i]^2, na.rm=TRUE))
    
    for (j in 1:length(donor))
    {
        temp <- get(paste("alpha", treated[i], sep=""))
        rmspe.d[j,i] <- sqrt(mean(temp[temp$date<0, 2+j]^2,na.rm=TRUE ))
    }	
}

## Post RMSPE
postrmspe.tr <- numeric(length(treated))
postrmspe.d <- matrix(0, length(donor),length(treated))
for (i in 1:length(treated))
{
    postrmspe.tr[i] <- sqrt(mean(alpha.treated[11:23, i]^2, na.rm=TRUE))
    
    for (j in 1:length(donor))
    {
        temp <- get(paste("alpha", treated[i], sep=""))
        postrmspe.d[j,i] <- sqrt(mean(temp[temp$date>0, 2+j]^2,na.rm=TRUE ))
    }	
}

# 5.1 e(avg_pre_rmspe_p) A measure of fit. Concerning if significant.
## The proportion of placebos that have a pre- treatment RMSPE at least as large as the average of the treated units. 
rmspe.fit=numeric(length(treated))
for (i in 1: length(treated)) {
    
    rmspe.fit[i] <- length(which(rmspe.d[,i]>rmspe.tr[i]))/length(donor)
    
}
postestimation<-data.frame(cown=treated,rmspe.fit=rmspe.fit)

## 5.2  e(pval_joint_post_t)
## The proportion of placebos that have a ratio of post-treatment RMSPE over pre-treatment RMSPE at least as large as the average ratio for the treated units.
postestimation$pval.joint.post_t<-numeric(length(treated))
for (i in 1: length(treated)) {
    
    postestimation$pval.joint.post_t[i] <- length(which(postrmspe.d[,i]/rmspe.d[,i]>postrmspe.tr[i]/rmspe.tr[i]))/length(donor)
}

### Calculating P values, with all sample ###
# Random sample of 10000 #
ssize <- 10000
lead <- 10
alpha.pl <- matrix(0,ssize,lead)

for (t in 1:lead)
{
    for (i in 1:ssize)
    {
        temp <- sample(1:30, length(treated), replace=FALSE)
        alpha.pl[i,t] <- mean( c(
                                 alpha94[t+10,2+temp[1]],
                                 alpha140[t+10,2+temp[2]],
                                 alpha155[t+10,2+temp[3]],
                                 alpha160[t+10,2+temp[4]],
                                 alpha161[t+10,2+temp[5]],
                                 alpha200[t+10,2+temp[6]],
                                 alpha375[t+10,2+temp[7]],
                                 alpha385[t+10,2+temp[8]],
                                 alpha390[t+10,2+temp[9]],
                                 alpha640[t+10,2+temp[10]],
                                 alpha666[t+10,2+temp[11]],
                                 alpha750[t+10,2+temp[12]],
                                 alpha770[t+10,2+temp[13]],
                                 alpha780[t+10,2+temp[14]],
                                 alpha800[t+10,2+temp[15]],
                                 alpha840[t+10,2+temp[16]],
                                 alpha841[t+10,2+temp[17]],
                                 alpha900[t+10,2+temp[18]],
                                 alpha920[t+10,2+temp[19]]), na.rm=TRUE)
        
    }
}

alpha.l <- alpha.treated[11:(10+lead), "ave.tr"]

for (t in 1:lead)
{
    alpha.l[t] <- mean( c(
                          alpha94[t+10,2],
                          
                          alpha140[t+10,2],
                          alpha155[t+10,2],
                          alpha160[t+10,2],
                          alpha161[t+10,2],
                          alpha200[t+10,2],
                          alpha375[t+10,2],
                          alpha385[t+10,2],
                          alpha390[t+10,2],
                          alpha640[t+10,2],
                          alpha666[t+10,2],
                          alpha750[t+10,2],
                          alpha770[t+10,2],
                          alpha780[t+10,2],
                          alpha800[t+10,2],
                          alpha840[t+10,2],
                          alpha841[t+10,2],
                          alpha900[t+10,2],
                          alpha920[t+10,2]), na.rm=TRUE)
}

pvalue <- numeric(lead)
for (t in 1:lead)
{
    pvalue[t] <- length(which(alpha.pl[,t]>alpha.l[t]))/ssize
}
pvalue

### 6. Plotting p-values
plot(1:10, pvalue, ylim=c(0,1), axes = FALSE,xlab = "", ylab = "",  pch=20, bg="grey", lwd=5, cex=2.5)
box(); grid()
axis(1, at = 1:10, labels = seq(1, 10, 1), cex.axis = 1.25)
axis(2, at = seq(0, 1, 0.1), cex.axis = 1, las = 2)
title(xlab = "Number of years after treatment", cex.lab = 1.5, line = 2.4)
title(ylab = "Average lead specific P-Values", cex.lab = 1.7, line = 2.6)
title("Lead Specific Significance Level (P-Values) for Military Spending")




######### FOREIGN POLICY  1960-2012  ###################
## Multiple Treatment Estimation ##

rm(list=ls())

### 1.  DATA clean up ###
setwd("~/Dropbox/Women Leaders and Conflict/Empirical Analysis/Rdata and scripts")
library(readstata13)
mt_fp1 <- read.dta13("mt_fp_2(60-14)_2915.dta")
colnames(mt_fp1)[4]=c("female")

# 1977 recoded as male=0, as Gandhi was in office till March 23 #
mt_fp1$female[(mt_fp1$cown==750 & mt_fp1$year==1977)] <- 0
#Meir only was in office till June 03, so I have recoded 1974 as having a male leader
mt_fp1$female[(mt_fp1$cown==666 & mt_fp1$year==1974)] <- 0
# canada only one year treatment coded as non-treatment
mt_fp1$female[(mt_fp1$cown==20 & mt_fp1$year==1993)] <- 0
# ecuador coded as zero 1997 
mt_fp1$female[(mt_fp1$cown==130 & mt_fp1$year==1997)] <- 0
# norway one year treatment coded as zero , 1981
mt_fp1$female[(mt_fp1$cown==385 & mt_fp1$year==1981)] <- 0
# India second treatment coded as zero #
mt_fp1$female[(mt_fp1$cown==750 & mt_fp1$year>1979 & mt_fp1$year<1985)] <- 0
# Pakistan second treatment coded as zero
mt_fp1$female[(mt_fp1$cown==770 & mt_fp1$year>1992 & mt_fp1$year<1997)] <- 0
# Sri lanka coded female, if the president is female 
mt_fp1$female[(mt_fp1$cown==780)] <- 0
mt_fp1$female[(mt_fp1$cown==780 & mt_fp1$year>1994 & mt_fp1$year<2005)] <- 1


### Duplicate countries with two treatments

# Original Argentina is cown 160, and a replica Argentina cown 161
df160 <- data.frame(matrix(ncol = 485, nrow = 53))
df160 <- mt_fp1[mt_fp1$cown==160,]
df160$cown=161
df160$countrycode_pw <- "Arg-new"
df160$female[(df160$year>1973 & df160$year<1977)] <- 0
mt_fp1$female[(mt_fp1$cown==160 & mt_fp1$year>2006)] <- 0

# Append df160
mt_fp1 <- rbind(mt_fp1, df160)


# Original Philippines is cown 840, and a replica Philippines is cown 841
df840 <- data.frame(matrix(ncol = 485, nrow = 53))
df840 <- mt_fp1[mt_fp1$cown==840,]
df840$cown=841
df840$countrycode_pw <- "Phil-new"
df840$female[(df840$year>1985 & df840$year<1993)] <- 0
mt_fp1$female[(mt_fp1$cown==840 & mt_fp1$year>2000)] <- 0

# Append df840
mt_fp1 <- rbind(mt_fp1, df840)

### 2. Data set-up for synth ###
#### Identify treated countries ###
macro <- aggregate(female ~ cown, data = mt_fp1, sum)
macro$donor <- 0
macro$treated <- 0
macro$donor[macro$female==0]=1
macro$treated[macro$female>0]=1
macro <- subset(macro, select=-c(female))

treated <- macro$cown[macro$treated==1]
donor <- macro$cown[macro$donor==1]

tr.start <- numeric(length(treated))
tr.end <- numeric(length(treated))

for (i in 1:length(treated))
{
    tr.years <- mt_fp1$year[(mt_fp1$female==1 & mt_fp1$cown==treated[i])]
    tr.start[i] <- min(tr.years)
    tr.end[i] <- max(tr.years)
}

mt_fp1 <- merge(mt_fp1, macro, by="cown")

predictors=c("polity2_20", "loggdppc_pw", "gdpgrowth_pw", "cinc_c", "lntenure_a")
mt_fp1$fpdiff<-numeric(dim(mt_fp1)[1])
mt_fp1$fpdiff=mt_fp1$idealpoint2_v-mt_fp1$idealpoint2_v[mt_fp1$cown==2]

#### 3. Running MULTIPLE TREATMENT SYNTH
tic()
for (i in 1:length(treated))
{
    dataprep.out<- 
        dataprep( 
            foo = mt_fp1, 
            predictors = predictors,
            predictors.op = "mean", 
            dependent = "fpdiff", 
            unit.variable = "cown", 
            time.variable = "year",  
            treatment.identifier = treated[i], 
            controls.identifier = donor, 
            time.predictors.prior = c((max(1960,tr.start[i]-5)):tr.start[i]), 
            time.optimize.ssr = c((max(1960,tr.start[i]-5)):tr.start[i]), 
            unit.names.variable = "countrycode_pw", 
            time.plot = c((max(1960,tr.start[i]-5)):tr.end[i]) 
        )
    
    synth.out <- synth(dataprep.out)
    synth.tables <- synth.tab( dataprep.res = dataprep.out, synth.res = synth.out) 
    assign(paste("W",treated[i], sep=""), (synth.tables$tab.w))
    assign(paste("V",treated[i], sep=""), (synth.out$solution.v))
    assign(paste("Pred",treated[i], sep=""), (synth.tables$tab.pred))
    
    assign(paste("Y1",treated[i], sep="."), (dataprep.out$Y1))
    assign(paste("Y1s",treated[i], sep="."), (dataprep.out$Y0plot%*%synth.out$solution.w))
    assign(paste("Alpha",treated[i], sep="."), (dataprep.out$Y1-dataprep.out$Y0plot%*%synth.out$solution.w))
    
    ### Run placebo for each treated unit
    for (j in 1:length(donor))
    {
       
        dataprep.out<- 
            dataprep( 
                foo = mt_fp1, 
                predictors = predictors,
                predictors.op = "mean", 
                dependent = "fpdiff", 
                unit.variable = "cown", 
                time.variable = "year",  
                treatment.identifier = donor[j], 
                controls.identifier = donor[-j], 
                time.predictors.prior = c((max(1960,tr.start[i]-5)):tr.start[i]), 
                time.optimize.ssr = c((max(1960,tr.start[i]-5)):tr.start[i]), 
                unit.names.variable = "countrycode_pw", 
                time.plot = c((max(1960,tr.start[i]-5)):tr.end[i]) 
            )
        
        synth.out <- synth(dataprep.out)
        synth.tables <- synth.tab( dataprep.res = dataprep.out, synth.res = synth.out) 
        assign(paste("Y1",treated[i],donor[j], sep="."), (dataprep.out$Y1))
        assign(paste("Y1s",treated[i],donor[j], sep="."), (dataprep.out$Y0plot%*%synth.out$solution.w))
        assign(paste("Alpha",treated[i],donor[j], sep="."), (dataprep.out$Y1-dataprep.out$Y0plot%*%synth.out$solution.w))
        
    } ## j loop
    #}
    
} ## i loop
toc()


### 4. Merging back outcomes with centered at 0(treatment year)
for (i in 1: length(treated))
{
    temp <- get(paste("Alpha",treated[i], sep="."))
    year <- as.numeric(rownames(temp))
    year <- year-tr.start[i]+1
    
    assign(paste("alpha",treated[i], sep=""), (data.frame(date=year, temp)))
    temp3 <- get(paste("alpha",treated[i], sep=""))
    
    for (j in 1: length(donor))
    {
        temp2 <- get(paste("Alpha",treated[i],donor[j], sep="."))
        temp3 <- cbind(temp3, temp2)
        colnames(temp3)[j+2] <- paste("Alpha", treated[i],donor[j], sep=".")
    }
    
    assign(paste("alpha",treated[i], sep=""), (temp3))
}

## Y
for (i in 1: length(treated))
{
    temp <- get(paste("Y1",treated[i], sep="."))
    year <- as.numeric(rownames(temp))
    year <- year-tr.start[i]+1
    
    assign(paste("y1",treated[i], sep="."), (data.frame(date=year, temp)))
    temp3 <- get(paste("y1",treated[i], sep="."))
    
    for (j in 1: length(donor))
    {
        temp2 <- get(paste("Y1",treated[i],donor[j], sep="."))
        temp3 <- cbind(temp3, temp2)
        colnames(temp3)[j+2] <- paste("Y1", treated[i],donor[j], sep=".")
    }
    
    assign(paste("y1",treated[i], sep="."), (temp3))
}

## Ys
for (i in 1: length(treated))
{
    temp <- get(paste("Y1s",treated[i], sep="."))
    year <- as.numeric(rownames(temp))
    year <- year-tr.start[i]+1
    
    assign(paste("y1s",treated[i], sep="."), (data.frame(date=year, temp)))
    temp3 <- get(paste("y1s",treated[i], sep="."))
    
    for (j in 1: length(donor))
    {
        temp2 <- get(paste("Y1s",treated[i],donor[j], sep="."))
        temp3 <- cbind(temp3, temp2)
        colnames(temp3)[j+2] <- paste("Y1s", treated[i],donor[j], sep=".")
    }
    
    colnames(temp3)[2] <- paste("w", treated[i], sep="")
    assign(paste("y1s",treated[i], sep="."), (temp3))
}

## All alphas
alpha <- merge(alpha93, alpha94, by="date", all=T)
alpha <- merge(alpha, alpha95, by="date", all=T)
alpha <- merge(alpha, alpha140, by="date", all=T)
alpha <- merge(alpha, alpha155, by="date", all=T)
alpha <- merge(alpha, alpha160, by="date", all=T)
alpha <- merge(alpha, alpha161, by="date", all=T)
alpha <- merge(alpha, alpha200, by="date", all=T)
alpha <- merge(alpha, alpha375, by="date", all=T)
alpha <- merge(alpha, alpha385, by="date", all=T)
alpha <- merge(alpha, alpha390, by="date", all=T)
alpha <- merge(alpha, alpha640, by="date", all=T)
alpha <- merge(alpha, alpha666, by="date", all=T)
alpha <- merge(alpha, alpha750, by="date", all=T)
alpha <- merge(alpha, alpha770, by="date", all=T)
alpha <- merge(alpha, alpha780, by="date", all=T)
alpha <- merge(alpha, alpha800, by="date", all=T)
alpha <- merge(alpha, alpha840, by="date", all=T)
alpha <- merge(alpha, alpha841, by="date", all=T)
alpha <- merge(alpha, alpha900, by="date", all=T)
alpha <- merge(alpha, alpha920, by="date", all=T)

coln.treated <- character(length(treated))
coln.donor <- matrix(0,length(donor),length(treated))

for (i in 1:length(treated))
{
    coln.treated[i] <-paste("X", treated[i], sep="") 
    
    for (j in 1:length(donor))
    {
        coln.donor[j,i] <-paste("Alpha", treated[i], donor[j], sep=".")
    }
}
coln.donor <- as.vector(coln.donor)

alpha.treated <- alpha[,colnames(alpha) %in% coln.treated]
alpha.treated$date <- alpha$date

alpha.donor <- alpha[,colnames(alpha) %in% coln.donor]
alpha.donor$date <- alpha$date

alpha.treated$ave.tr <- rowMeans(alpha.treated[,1:length(treated)], na.rm=TRUE)
alpha.donor$ave.d <- rowMeans(alpha.donor[,1:length(coln.donor)], na.rm=TRUE)

## 5. Pre RMSPE comparison ##
## RMSPE comparison ##
rmspe.tr <- numeric(length(treated))
rmspe.d <- matrix(0, length(donor),length(treated))


for (i in 1:length(treated))
{
    rmspe.tr[i] <- sqrt(mean(alpha.treated[alpha.treated$date<0, i]^2, na.rm=TRUE))
    
    for (j in 1:length(donor))
    {
        temp <- get(paste("alpha", treated[i], sep=""))
        rmspe.d[j,i] <- sqrt(mean(temp[temp$date<=0, 2+j]^2,na.rm=TRUE ))
    }	
}

## Post RMSPE
postrmspe.tr <- numeric(length(treated))
postrmspe.d <- matrix(0, length(donor),length(treated))
for (i in 1:length(treated))
{
    postrmspe.tr[i] <- sqrt(mean(alpha.treated[1:5, i]^2, na.rm=TRUE))
    
    for (j in 1:length(donor))
    {
        temp <- get(paste("alpha", treated[i], sep=""))
        postrmspe.d[j,i] <- sqrt(mean(temp[temp$date>0, 2+j]^2,na.rm=TRUE ))
    }	
}

# 5.1 e(avg_pre_rmspe_p) A measure of fit. Concerning if significant.
## The proportion of placebos that have a pre- treatment RMSPE at least as large as the average of the treated units. 
rmspe.fit=numeric(length(treated))
for (i in 1: length(treated)) {
    
    rmspe.fit[i] <- length(which(rmspe.d[,i]>rmspe.tr[i]))/length(donor)
    
}
postestimation<-data.frame(cown=treated,rmspe.fit=rmspe.fit)

## 5.2  e(pval_joint_post_t)
## The proportion of placebos that have a ratio of post-treatment RMSPE over pre-treatment RMSPE at least as large as the average ratio for the treated units.
postestimation$pval.joint.post_t<-numeric(length(treated))
for (i in 1: length(treated)) {
    
    postestimation$pval.joint.post_t[i] <- length(which(postrmspe.d[,i]/rmspe.d[,i]>postrmspe.tr[i]/rmspe.tr[i]))/length(donor)
}


### Calculating P values, with all sample(16) ###
# Random sample of 10000 #
ssize <- 10000
lead <- 10
alpha.pl <- matrix(0,ssize,lead)

for (t in 1:10)
{
    for (i in 1:ssize)
    {
        temp <- sample(1:35, length(treated), replace=FALSE)
        alpha.pl[i,t] <- mean( c(alpha93[t+5,2+temp[1]], ## change to t+10 if 10 year pretreatment
                                 alpha94[t+5,2+temp[2]],
                                 alpha95[t+5,2+temp[3]],
                                 alpha140[t+5,2+temp[4]],
                                 alpha155[t+5,2+temp[5]],
                                 alpha160[t+5,2+temp[6]],
                                 alpha161[t+5,2+temp[7]],
                                 alpha200[t+5,2+temp[8]],
                                 alpha375[t+5,2+temp[9]],
                                 alpha385[t+5,2+temp[10]],
                                 alpha390[t+5,2+temp[11]],
                                 alpha640[t+5,2+temp[12]],
                                 alpha666[t+5,2+temp[13]],
                                 alpha750[t+5,2+temp[14]],
                                 alpha770[t+5,2+temp[15]],
                                 alpha780[t+5,2+temp[16]],
                                 alpha800[t+5,2+temp[17]],
                                 alpha840[t+5,2+temp[18]],
                                 alpha841[t+5,2+temp[19]],
                                 alpha900[t+5,2+temp[20]],
                                 alpha920[t+5,2+temp[21]]), na.rm=TRUE)
        
    }
}

alpha.l <- alpha.treated[6:(5+lead), "ave.tr"]

for (t in 1:lead)
{
    alpha.l[t] <- mean( c(alpha93[t+5,2],
                          alpha94[t+5,2],
                          alpha95[t+5,2],
                          alpha140[t+5,2],
                          alpha155[t+5,2],
                          alpha160[t+5,2],
                          alpha161[t+5,2],
                          alpha200[t+5,2],
                          alpha375[t+5,2],
                          alpha385[t+5,2],
                          alpha390[t+5,2],
                          alpha640[t+5,2],
                          alpha666[t+5,2],
                          alpha750[t+5,2],
                          alpha780[t+5,2],
                          alpha800[t+5,2],
                          alpha840[t+5,2],
                          alpha841[t+5,2],
                          alpha900[t+5,2],
                          alpha920[t+5,2]), na.rm=TRUE)
}

pvalue <- numeric(lead)
for (t in 1:lead)
{
    pvalue[t] <- length(which(alpha.pl[,t]>alpha.l[t]))/ssize
}
pvalue

### 6. Plotting p-values
plot(1:10, pvalue, ylim=c(0,1), axes = FALSE,xlab = "", ylab = "",  pch=20, bg="grey", lwd=5, cex=2.5)
box(); grid()
axis(1, at = 1:10, labels = seq(1, 10, 1), cex.axis = 1.25)
axis(2, at = seq(0, 1, 0.1), cex.axis = 1, las = 2)
title(xlab = "Number of years after treatment", cex.lab = 1.5, line = 2.4)
title(ylab = "Average lead specific P-Values", cex.lab = 1.7, line = 2.6)
title("Lead Specific Significance Level (P-Values) for Ideal-Point estimaes")

