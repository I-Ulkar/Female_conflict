####       ULKAR IMAMVERDI AND PATRICK SHEA    ###
###       Female Leaders and Foreign Policy    ###
###         Part 3: Appendix                   ###
##################################################

rm(list=ls())
library(Synth)
library(foreign)
library(tictoc)
library(readstata13)
library(ggplot2)
########Synthetic Control Analysis of Gender and Conflict          (Appendix C) ####### 
rm(list=ls())
setwd("~/Dropbox/Women Leaders and Conflict/Empirical Analysis/Rdata and scripts")

library(haven)
ukcb1 <- read_dta("UK_cb1.dta")

library(Synth)
library(foreign)

### 1. Prepping data
macro <- aggregate(female ~ cown, data = ukcb1, sum)
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
    tr.years <- ukcb1$year[(ukcb1$female==1 & ukcb1$cown==treated[i])]
    tr.start[i] <- min(tr.years)
    tr.end[i] <- max(tr.years)
}
ukcb1 <- merge(ukcb1, macro, by="cown")
ukcb1$partycontrol<-as.numeric(ukcb1$gtm_pr_k)
predictors=c("mean_weis_confcc2_cc1_k", "m_woman2_k", "gtm_pr_k", "newleftright_k",   "parliament_k", "allies_k", "majpow_k", "cap_k")
### 2. Estimating Synthetic Counterfactual for UK 
ukcb.dp <- dataprep(
    foo=ukcb1, 
    predictors = predictors,
    special.predictors = list(list("mean_weis_confcc1_cc2_k", 1970:1979, c("mean"))),
    predictors.op = "mean", 
    dependent = "mean_weis_confcc1_cc2_k", 
    unit.variable = "cown", 
    time.variable = "year",  
    treatment.identifier = 200, 
    controls.identifier = donor, 
    time.predictors.prior = c(1970:1978), 
    time.optimize.ssr = c(1970:1978), 
    unit.names.variable = "country", 
    time.plot = c(1970:tr.end[i]) 
)

ukcb.synth <- synth(ukcb.dp, "ALL")

synthfp.tables <- synth.tab( dataprep.res = ukcb.dp, synth.res = ukcb.synth) 
ukcb.pred <- ukcb.dp$Y0plot %*% ukcb.synth$solution.w ####Y1s in Previous code
ukcb.alpha<- ukcb.dp$Y1plot - ukcb.pred

###. 3 Plotting the results UK CB
plot(1970:1990, ukcb.dp$Y1plot, ylim = c(-0.5, 2.5),  type = "n", lwd = 2, xlab = "", ylab = "")
box(); grid()
abline(v = 1979, lwd = 1)  
lines(1970:1990, ukcb.dp$Y1plot, lwd = 2)
lines(1970:1990, ukcb.pred, lwd = 2, lty = 2)
title(xlab = "Year", cex.lab = 1.5, line = 2)
title(ylab = "Conflict Behavior", cex.lab = 2.5, line = 2)
rug(ukcb1$mean_weis_confcc1_cc2_k, side = 2, lwd=1, ticksize = .025)
legend(list(x = 1970,y = 2.9), bty = "n", c("UK", "Synthetic UK"), lwd = 2, col = "black", lty = c(1, 2, -1), cex = 1.3, y.intersp=0.3)
text(1974, -0.2, "First year of Thatcher's tenure", cex = 1.2)
arrows(1978, -0.2, 1979, -0.2, length = .2)
text(1987, 2.2, bquote(paste("MSPE" == .(round(ukcb.synth$loss.v, 5)))), cex = 1.5)


### 4. Placebo cases ##
for (j in 1:length(donor))
{
    ukcb_pl.dp<- 
        dataprep( 
            foo = ukcb1, 
            predictors = predictors,
            special.predictors = list(list("mean_weis_confcc1_cc2_k", 1978, c("mean"))),
            predictors.op = "mean", 
            dependent = "mean_weis_confcc1_cc2_k", 
            unit.variable = "cown", 
            time.variable = "year",  
            treatment.identifier = donor[j], 
            controls.identifier = donor[-j], 
            time.predictors.prior = c(1970:1978), 
            time.optimize.ssr = c(1970:1978), 
            unit.names.variable = "country", 
            time.plot = c(1970:tr.end[i]) 
        ) 
    ukcb_synth_pl <- synth(ukcb_pl.dp, "ALL")
    synthfp_pl.tables <- synth.tab( dataprep.res = ukcb_pl.dp, synth.res = ukcb_synth_pl, round.digit = 5) 
    assign(paste("W",treated[i],donor[j], sep="."), (synthfp_pl.tables$tab.w))
    assign(paste("V",treated[i],donor[j], sep="."), (ukcb_synth_pl$solution.v))
    assign(paste("Pred",treated[i],donor[j], sep="."), (synthfp_pl.tables$tab.pred))
    
    assign(paste("Y1",treated[i],donor[j], sep="."), (ukcb_pl.dp$Y1))
    assign(paste("Y1s",treated[i],donor[j], sep="."), (ukcb_pl.dp$Y0plot%*%ukcb_synth_pl$solution.w))
    assign(paste("Alpha",treated[i],donor[j], sep="."), (ukcb_pl.dp$Y1-ukcb_pl.dp$Y0plot%*%ukcb_synth_pl$solution.w))
}

### 5. Post Estimation results ###
## pre-RMSPE comparison ##
uk.prermspe <- sqrt(mean(ukcb.alpha[1:9]^2, na.rm=TRUE))
uk.postrmspe <- sqrt(mean(ukcb.alpha[10:21]^2, na.rm=TRUE))
uk.prermspe.d <- numeric(length(donor))
uk.postrmspe.d <- numeric(length(donor))

for (j in 1:22)
{
    temp <- get(paste("Alpha", treated[i], donor[j], sep="."))
    uk.prermspe.d[j] <- sqrt(mean(temp[1:9]^2,na.rm=TRUE ))
    uk.postrmspe.d[j] <- sqrt(mean(temp[10:21]^2,na.rm=TRUE ))
}

# 5.1 e(avg_pre_rmspe_p) A measure of fit. Concerning if significant.
## The proportion of placebos that have a pre- treatment RMSPE at least as large as the average of the treated units. 
rmspe.fit <- length(which(uk.prermspe.d>uk.prermspe))/length(donor)


## 5.2  e(pval_joint_post_t)
## The proportion of placebos that have a ratio of post-treatment RMSPE over pre-treatment RMSPE at least as large as the average ratio for the treated units.
pval.joint.post_t <- length(which(uk.postrmspe.d/uk.prermspe.d>uk.postrmspe/uk.prermspe))/length(donor)

alphas <-ukcb.alpha
for (j in 1: length(donor))
{
    temp1 <- get(paste("Alpha",treated[i],donor[j], sep="."))
    alphas <- cbind(alphas, temp1)
    colnames(alphas)[j+1] <- paste(treated[i],donor[j], sep=".")
}

## 5.3  e(pvals_t)
# A vector of the proportions of placebo pseudo t-statistics (unit’s effect divided by its pre-treatment RMSPE) that are at least as largeas the main pseudo t-statistic for each post-treatment period.
pvals_t <- numeric(12) 
for (t in 1:12)
{
    pvals_t[t] <-  length(which(abs(alphas[t+8,-1]/uk.prermspe.d)>abs(alphas[t+8,1]/uk.prermspe)))/length(donor)
}

plot(1:12, pvals_t, ylim=c(0,1), axes = FALSE,xlab = "", ylab = "",  pch=20, bg="grey", lwd=5, cex=2.5)
box(); grid()
axis(1, at = 1:12, labels = seq(1, 12, 1), cex.axis = 1.25)
axis(2, at = seq(0, 1, 0.1), cex.axis = 1, las = 2)
title(xlab = "Number of years after treatment", cex.lab = 1.5, line = 2.4)
title(ylab = "Lead specific P-Values", cex.lab = 1.9, line = 2.6)
text(pvals_t, labels = round(pvals_t, 2), pos = 3)

#with(LifeCycleSavings[1:9,], 
### 6. Plotting placebo cases and UK CB
plot(1970:1990, ukcb.alpha, type = "n", ylim = c(-2, 2), lwd = 2, xlab = "", ylab = "", axes = FALSE) 
box(); grid()
for(j in 1: 22){
    lines(1970:1990, alphas[, j], type = "l", col = "gray80", lwd = 2)}
lines(1970:1990, alphas[,1 ], lwd = 2) #INDIA
abline(v=1979, lwd = 1)
abline(h=0, lwd = 1)
axis(1, at = 1970:1990, labels = seq(1970, 1990, 1), cex.axis = 1.25)
axis(2, at = seq(-2, 2, 1), cex.axis = 1.25, las = 2)
title(xlab = "Year", cex.lab = 1.5, line = 2)
title(ylab = "Difference in Conflict Behavior", cex.lab = 1.5, line = 2)
text(1975, -1.5, "Treatment time", cex = 1.4)
arrows(1977.2, -1.5, 1979, -1.5, length = .1)
legend(list(x = 1982,y = -0.5), ncol = 1, bty = "n", c("UK", "Placebo Countries"), lty=c(1,1), lwd = 2, col = c("black", "gray80"), cex = 1.25, y.intersp=0.3)


#country weights Table C10 
synthfp.tables$tab.w
#Balance Table C11
synthfp.tables$tab.pred





######### 1. Single treatment analysis  #################

### ### ISRAEL Military Spending with Time Covariates             (Appendix D.1 )  ############
rm(list=ls())
library(haven)
isrms1 <- read_dta("isr53-73fem_ms504r.dta")
library(Synth)

colnames(isrms1)[4]=c("female")
#isrfp1[isrfp1$cown == 666, c("logsmilex_f2", "lagdv")]

### 1. Prepping data
macro <- aggregate(female ~ cown, data = isrms1, sum)
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
    tr.years <- isrms1$year[(isrms1$female==1 & isrms1$cown==treated[i])]
    tr.start[i] <- min(tr.years)
    tr.end[i] <- max(tr.years)
}
isrms1 <- merge(isrms1, macro, by="cown")
predictors = c( "polity2_20","wyrs", "cyrs" ,"gdpgrowth_pw", "gdppc_pw", "women_legst", "cinc_c", "logtrade", "logthreat", "allies", "issues", "rival")

### 2. Estimating Synthetic Counterfactual for Israel 
isrms.dp <- dataprep(
    foo=isrms1, 
    predictors = predictors,
    predictors.op = "mean", 
    dependent = "logsmilex_f2", 
    special.predictors = list(list("logsmilex_f2",1968, c("mean"))),
    unit.variable = "cown", 
    time.variable = "year",  
    treatment.identifier = 666, 
    controls.identifier = donor, 
    time.predictors.prior = c(1956:1968), 
    time.optimize.ssr = c(1956:1968), 
    unit.names.variable = "country", 
    time.plot = c(1956:tr.end[i]) 
)
isrms.synth <- synth(isrms.dp, "ALL")

synthms.tables <- synth.tab( dataprep.res = isrms.dp, synth.res = isrms.synth) 
isrms.pred <- isrms.dp$Y0plot %*% isrms.synth$solution.w ####Y1s in Previous code
isrms.alpha<- isrms.dp$Y1plot - isrms.pred

## 3. Plot results
plot(1956:1973, isrms.dp$Y1plot, ylim=c(3, 13), type = "n", lwd = 2, xlab = "", ylab = "")
box(); grid()
abline(v = 1969, lwd = 1)

lines(1956:1973, isrms.dp$Y1plot, lwd = 2)
lines(1956:1973, isrms.pred, lwd = 2, lty = 2)
title(xlab = "Year", cex.lab = 1.5, line = 2)
title(ylab = "Military Spending", cex.lab = 1, line = 2)
legend(list(x = 1955,y = 13), bty = "n", c("Israel", "Synthetic Israel"), lwd = 2, col = "black", lty = c(1, 2, -1), cex = 1.5, y.intersp=0.5)
text(1963, 5, "First year of Meir's tenure", cex = 1.7)
text(1971.8, 11.8, bquote(paste("MSPE" == .(round(isrms.synth$loss.v, 4)))), cex = 1.5)
arrows(1967, 5, 1969, 5, length = .1)
rug(isrms1$logsmilex_f2, side = 2, lwd=1, ticksize = .025)

# #pdf("isr_ms.time.pdf")  

### 4. Placebo cases ##
### Run placebo for each donor unit
for (j in 1:length(donor))
{
    isrms_pl.dp<- 
        dataprep( 
            foo = isrms1, 
            predictors = predictors,
            predictors.op = "mean", 
            dependent = "logsmilex_f2", 
            special.predictors = list(list("logsmilex_f2",1968, c("mean"))),
            unit.variable = "cown", 
            time.variable = "year",  
            treatment.identifier = donor[j], 
            controls.identifier = donor[-j], 
            time.predictors.prior = c(1956:1968), 
            time.optimize.ssr = c(1956:1968), 
            unit.names.variable = "country", 
            time.plot = c(1956:tr.end[i]) 
        ) 
    
    isrms_synth_pl <- synth(isrms_pl.dp)
    synthms_pl.tables <- synth.tab( dataprep.res = isrms_pl.dp, synth.res = isrms_synth_pl, round.digit = 5) 
    assign(paste("W",treated[i],donor[j], sep="."), (synthms_pl.tables$tab.w))
    assign(paste("V",treated[i],donor[j], sep="."), (isrms_synth_pl$solution.v))
    assign(paste("Pred",treated[i],donor[j], sep="."), (synthms_pl.tables$tab.pred))
    
    assign(paste("Y1",treated[i],donor[j], sep="."), (isrms_pl.dp$Y1))
    assign(paste("Y1s",treated[i],donor[j], sep="."), (isrms_pl.dp$Y0plot%*%isrms_synth_pl$solution.w))
    assign(paste("Alpha",treated[i],donor[j], sep="."), (isrms_pl.dp$Y1-isrms_pl.dp$Y0plot%*%isrms_synth_pl$solution.w))
}


### 5. Post Estimation results ###

## pre-RMSPE comparison ##
isr.prermspe <- sqrt(mean(isrms.alpha[1:13]^2, na.rm=TRUE))
isr.postrmspe <- sqrt(mean(isrms.alpha[14:18]^2, na.rm=TRUE))
isr.prermspe.d <- numeric(length(donor))
isr.postrmspe.d <- numeric(length(donor))

for (j in 1:length(donor))
{
    temp <- get(paste("Alpha", treated[i], donor[j], sep="."))
    isr.prermspe.d[j] <- sqrt(mean(temp[1:13]^2,na.rm=TRUE ))
    isr.postrmspe.d[j] <- sqrt(mean(temp[14:18]^2,na.rm=TRUE ))
}

# 5.1 e(avg_pre_rmspe_p) A measure of fit. Concerning if significant.
## The proportion of placebos that have a pre- treatment RMSPE at least as large as the average of the treated units. 
rmspe.fit <- length(which(isr.prermspe.d>isr.prermspe))/length(donor)


## 5.2  e(pval_joint_post_t)
## The proportion of placebos that have a ratio of post-treatment RMSPE over pre-treatment RMSPE at least as large as the average ratio for the treated units.
pval.joint.post_t <- length(which(isr.postrmspe.d/isr.prermspe.d>isr.postrmspe/isr.prermspe))/length(donor)

alphas <-isrms.alpha
for (j in 1: length(donor))
{
    temp1 <- get(paste("Alpha",treated[i],donor[j], sep="."))
    alphas <- cbind(alphas, temp1)
    colnames(alphas)[j+1] <- paste(treated[i],donor[j], sep=".")
}

## 5.3  e(pvals_t)
# A vector of the proportions of placebo pseudo t-statistics (unit’s effect divided by its pre-treatment RMSPE) that are at least as largeas the main pseudo t-statistic for each post-treatment period.
pvals_t <- numeric(5) 
for (t in 1:5)
{
    pvals_t[t] <-  length(which(abs(alphas[t+13,-1]/isr.prermspe.d)>abs(alphas[t+13,1]/isr.prermspe)))/length(donor)
}

plot(1:5, pvals_t, ylim=c(0,1), axes = FALSE,xlab = "", ylab = "",  pch=20, bg="grey", lwd=5, cex=2.5)
box(); grid()
axis(1, at = 1:5, labels = seq(1, 5, 1), cex.axis = 1.25)
axis(2, at = seq(0, 1, 0.1), cex.axis = 1, las = 2)
title(xlab = "Number of years after treatment", cex.lab = 1.5, line = 2.4)
title(ylab = "Lead specific P-Values", cex.lab = 2.2, line = 2.6)
text(pvals_t, labels = round(pvals_t, 2), pos = 3, cex = 1.8)



### 6. Plotting placebo cases and Israel
plot(1956:1973, isrms.alpha, type = "n", ylim = c(-3, 3), lwd = 2, xlab = "", ylab = "", axes = FALSE) 
box(); grid()
for(j in 3: 23){
    lines(1956:1973, alphas[, j], type = "l", col = "gray80", lwd = 2)
} 
lines(1956:1973, alphas[,1 ], lwd = 2) #Israel 
abline(v=1969, lwd = 1)
abline(h=0, lwd = 1)
axis(1, at = 1956:1973, labels = seq(1956, 1973, 1), cex.axis = 1.25)
axis(2, at = seq(-3, 3, 1), cex.axis = 1.25, las = 2)
title(xlab = "Year", cex.lab = 1.5, line = 2)
title(ylab = "Difference in Military Spending", cex.lab = 1.5, line = 2)
text(1964, -2.5, "Treatment time", cex = 1.2)
arrows(1966, -2.5, 1969, -2.5, length = .1)
legend(list(x = 1956,y = 3.3), ncol = 1, bty = "n", c("Israel", "Placebo Countries"), lty=c(1,1), lwd = 2, col = c("black", "gray80"), cex = 1.25, y.intersp=0.5)

#country weights Table
synthms.tables$tab.w
#Balance Table
synthms.tables$tab.pred





### ### Denmark Military spending (Female 2011, 2012, 2013, 2014) (Appendix D.2) ### ### ####

#setwd("~/Dropbox/Women Leaders and Conflict/Empirical Analysis/Rdata and scripts")
mt_ms1 <- read.dta13("mt_ms 2695.dta")

#print(mt_ms1[is.na(mt_ms1$women_legst),c("cown", "year")])

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


#Delete countries with female leaders at any time in sample
mt_ms1[(mt_ms1$female==1), c("cown", "year")]

mtms <- mt_ms1[(mt_ms1$female==0 | mt_ms1$cown==390), ]
m <- aggregate(female ~ cown, data = mt_ms1 , sum)
m <- m[m$female==0,]


mtms <- mt_ms1[mt_ms1$cown %in% c(m$cown, 390),]

#drop years before 2000
mtms <- mtms[mtms$year>2000,]

#Check which countries missing data for women_leg
print(mtms[is.na(mtms$women_legst),c("cown", "year")])
mtms <- mtms[mtms$cown!=775,]
mtms <- mtms[mtms$cown!=101,]

### 1. Prepping data
macro <- aggregate(female ~ cown, data = mtms, sum)
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
    tr.years <- mtms$year[(mtms$female==1 & mtms$cown==treated[i])]
    tr.start[i] <- min(tr.years)
    tr.end[i] <- max(tr.years)
}
mtms <- merge(mtms, macro, by="cown")

### 2. Estimating Synthetic Counterfactual for Denmark (390) 

predictors2 = c( "polity2_20", "gdpgrowth_pw", "gdppc_pw", "women_legst", "cinc_c", "logtrade", "logthreat_f2", "allies", "issues", "rival")
mtms.dp <- dataprep(
    foo=mtms, 
    predictors = predictors2,
    special.predictors = list(list("logsmilex_f2", 2010, c("mean"))),
    predictors.op = "mean", 
    dependent = "logsmilex_f2", 
    unit.variable = "cown", 
    time.variable = "year",  
    treatment.identifier = treated, 
    controls.identifier = donor, 
    time.predictors.prior = c(2001:2010), 
    time.optimize.ssr = c(2001:2010), 
    unit.names.variable = "country", 
    time.plot = c(2001:tr.end[i]) 
)
mtms.synth <- synth(mtms.dp)

synthms2.tables <- synth.tab( dataprep.res = mtms.dp, synth.res = mtms.synth) 
mtms.pred <- mtms.dp$Y0plot %*% mtms.synth$solution.w ####Y1s in Previous code
mtms.alpha<- mtms.dp$Y1plot - mtms.pred

## 3. Plot results
plot(2001:2014, mtms.dp$Y1plot, ylim=c(3, 13), type = "n", lwd = 2, xlab = "", ylab = "")
box(); grid()
abline(v = 2011, lwd = 1)
lines(2001:2014, mtms.dp$Y1plot, lwd = 2)
lines(2001:2014, mtms.pred, lwd = 2, lty = 2)
title(xlab = "Year", cex.lab = 1.5, line = 2)
title(ylab = "Military Spending", cex.lab = 1, line = 2)
legend(list(x = 2000,y = 13), bty = "n", c("Denmark", "Synthetic Denmark"), lwd = 2, col = "black", lty = c(1, 2, -1), cex = 1.5, y.intersp=0.5)
text(2006, 5, "First year of Helle's tenure", cex = 1.4)
text(2012.7, 11.8, bquote(paste("MSPE" == .(round(mtms.synth$loss.v, 3)))), cex = 1.3)
arrows(1967, 5, 1969, 5, length = .1)
rug(mtms$logsmilex_f2, side = 2, lwd=1, ticksize = .025)
#pdf("Denmark_ms1_(406).pdf") 

### 4. Placebo cases ##
### Run placebo for each donor unit
for (j in 1:length(donor))
{
    dnmms_pl.dp<- 
        dataprep( 
            foo = mtms, 
            predictors = predictors2,
            predictors.op = "mean", 
            dependent = "logsmilex_f2", 
            special.predictors = list(list("logsmilex_f2",2010, c("mean"))),
            unit.variable = "cown", 
            time.variable = "year",  
            treatment.identifier = donor[j], 
            controls.identifier = donor[-j], 
            time.predictors.prior = c(2001:2010), 
            time.optimize.ssr = c(2001:2010), 
            unit.names.variable = "country", 
            time.plot = c(2001:tr.end[i]) 
        ) 
    
    dnmms_synth_pl <- synth(dnmms_pl.dp)
    synthms_pl.tables <- synth.tab( dataprep.res = dnmms_pl.dp, synth.res = dnmms_synth_pl, round.digit = 5) 
    assign(paste("W",treated[i],donor[j], sep="."), (synthms_pl.tables$tab.w))
    assign(paste("V",treated[i],donor[j], sep="."), (dnmms_synth_pl$solution.v))
    assign(paste("Pred",treated[i],donor[j], sep="."), (synthms_pl.tables$tab.pred))
    
    assign(paste("Y1",treated[i],donor[j], sep="."), (dnmms_pl.dp$Y1))
    assign(paste("Y1s",treated[i],donor[j], sep="."), (dnmms_pl.dp$Y0plot%*%dnmms_synth_pl$solution.w))
    assign(paste("Alpha",treated[i],donor[j], sep="."), (dnmms_pl.dp$Y1-dnmms_pl.dp$Y0plot%*%dnmms_synth_pl$solution.w))
}

### 5. Post Estimation results ###

## pre-RMSPE comparison ##
dnm.prermspe <- sqrt(mean(mtms.alpha[1:10]^2, na.rm=TRUE))
dnm.postrmspe <- sqrt(mean(mtms.alpha[11:14]^2, na.rm=TRUE))
dnm.prermspe.d <- numeric(length(donor))
dnm.postrmspe.d <- numeric(length(donor))

for (j in 1:length(donor))
{
    temp <- get(paste("Alpha", treated[i], donor[j], sep="."))
    dnm.prermspe.d[j] <- sqrt(mean(temp[1:10]^2,na.rm=TRUE ))
    dnm.postrmspe.d[j] <- sqrt(mean(temp[11:14]^2,na.rm=TRUE ))
}


# 5.1 e(avg_pre_rmspe_p) A measure of fit. Concerning if significant.
## The proportion of placebos that have a pre- treatment RMSPE at least as large as the average of the treated units. 
rmspe.fit <- length(which(dnm.prermspe.d>dnm.prermspe))/length(donor)

## 5.2  e(pval_joint_post_t)
## The proportion of placebos that have a ratio of post-treatment RMSPE over pre-treatment RMSPE at least as large as the average ratio for the treated units.
pval.joint.post_t <- length(which(dnm.postrmspe.d/dnm.prermspe.d>dnm.postrmspe/dnm.prermspe))/length(donor)

alphas <-mtms.alpha
for (j in 1: length(donor))
{
    temp1 <- get(paste("Alpha",treated[i],donor[j], sep="."))
    alphas <- cbind(alphas, temp1)
    colnames(alphas)[j+1] <- paste(treated[i],donor[j], sep=".")
}

## 5.3  e(pvals_t)
# A vector of the proportions of placebo pseudo t-statistics (unit’s effect divided by its pre-treatment RMSPE) that are at least as largeas the main pseudo t-statistic for each post-treatment period.
pvals_t <- numeric(4) 
for (t in 1:4)
{
    pvals_t[t] <-  length(which(abs(alphas[t+10,-1]/dnm.prermspe.d)>abs(alphas[t+10,1]/dnm.prermspe)))/length(donor)
}

plot(1:4, pvals_t, ylim=c(0,1), axes = FALSE,xlab = "", ylab = "",  pch=20, bg="grey", lwd=5, cex=2.5)
box(); grid()
axis(1, at = 1:4, labels = seq(1, 4, 1), cex.axis = 1.25)
axis(2, at = seq(0, 1, 0.1), cex.axis = 1, las = 2)
title(xlab = "Number of years after treatment", cex.lab = 1.5, line = 2.4)
title(ylab = "Lead specific P-Values", cex.lab = 2.2, line = 2.6)
text(pvals_t, labels = round(pvals_t, 2), pos = 1, cex = 1.8)
#pdf("Denmark_ms1_pvalue(406).pdf")

### 6. Plotting placebo cases and Denmark
plot(2001:2014, mtms.alpha, type = "n", ylim = c(-3, 3), lwd = 2, xlab = "", ylab = "", axes = FALSE) 
box(); grid()
for(j in 3: 29){
    lines(2001:2014, alphas[, j], type = "l", col = "gray80", lwd = 2)
} 
lines(2001:2014, alphas[,1 ], lwd = 2) #DENMARK
abline(v=2011, lwd = 1)
abline(h=0, lwd = 1)
axis(1, at = 2001:2014, labels = seq(2001, 2014, 1), cex.axis = 1.25)
axis(2, at = seq(-3, 3, 1), cex.axis = 1.25, las = 2)
title(xlab = "Year", cex.lab = 1.5, line = 2)
title(ylab = "Difference in Military Spending", cex.lab = 1.5, line = 2)
text(2005, -2.5, "Treatment time", cex = 1.4)
arrows(2007, -2.5, 2010, -2.5, length = .1)
legend(list(x = 2001,y = 3.3), ncol = 1, bty = "n", c("Denmark", "Placebo Countries"), lty=c(1,1), lwd = 2, col = c("black", "gray80"), cex = 1.25, y.intersp=0.5)
#pdf("Denmark_ms1pl_(406).pdf")

### 7. TABLES

synthms_pl.tables$tab.w #(country weights)
synthms2.tables$tab.pred #(Balance table)

### ### Norway (385) Military spending (Female 1986-1996)         (Appendix D.2) ### ### ####
#setwd("~/Dropbox/Women Leaders and Conflict/Empirical Analysis/Rdata and scripts")
mt_ms1 <- read.dta13("mt_ms 2695.dta")


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


#Delete countries with female leaders at any time in sample
mt_ms1[(mt_ms1$female==1), c("cown", "year")]

tmp <- aggregate(female ~ cown, data = mt_ms1 , sum)
m <- tmp[tmp$female==0,]
norw.ms<- mt_ms1[mt_ms1$cown %in% c(m$cown, 385),]

#drop years before 2000
norw.ms <- norw.ms[norw.ms$year>1975 & norw.ms$year<1997,]

#Check which countries missing data for women_leg in all pretreatment years
print(norw.ms[is.na(norw.ms$women_legst),c("cown", "year")])
norw.ms <- norw.ms[norw.ms$cown!=475,] 
norw.ms <- norw.ms[norw.ms$cown!=101,] #gdp is missing

### 1. Prepping data
macro <- aggregate(female ~ cown, data = norw.ms, sum)
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
    tr.years <- norw.ms$year[(norw.ms$female==1 & norw.ms$cown==treated[i])]
    tr.start[i] <- min(tr.years)
    tr.end[i] <- max(tr.years)
}
norw.ms <- merge(norw.ms, macro, by="cown")

### 2. Estimating Synthetic Counterfactual for Denmark (390) 

predictors2 = c( "polity2_20", "gdpgrowth_pw", "gdppc_pw", "women_legst", "cinc_c", "logtrade", "logthreat_f2", "allies", "issues", "rival")
norw.ms.dp <- dataprep(
    foo=norw.ms, 
    predictors = predictors2,
    special.predictors = list(list("logsmilex_f2", 1985, c("mean"))),
    predictors.op = "mean", 
    dependent = "logsmilex_f2", 
    unit.variable = "cown", 
    time.variable = "year",  
    treatment.identifier = treated, 
    controls.identifier = donor, 
    time.predictors.prior = c(1976:1985), 
    time.optimize.ssr = c(1976:1985), 
    unit.names.variable = "country", 
    time.plot = c(1976:tr.end[i]) 
)
norw.ms.synth <- synth(norw.ms.dp)

synthms.tables <- synth.tab( dataprep.res = norw.ms.dp, synth.res = norw.ms.synth) 
norw.ms.pred <- norw.ms.dp$Y0plot %*% norw.ms.synth$solution.w ####Y1s in Previous code
norw.ms.alpha<- norw.ms.dp$Y1plot - norw.ms.pred

## 3. Plot results
plot(1976:1996, norw.ms.dp$Y1plot, ylim=c(3, 13), type = "n", lwd = 2, xlab = "", ylab = "")
box(); grid()
abline(v = 1986, lwd = 1)

lines(1976:1996, norw.ms.dp$Y1plot, lwd = 2)
lines(1976:1996, norw.ms.pred, lwd = 2, lty = 2)
title(xlab = "Year", cex.lab = 1.5, line = 2)
title(ylab = "Military Spending", cex.lab = 1, line = 2)
legend(list(x = 1976,y = 14.5), bty = "n", c("Norway", "Synthetic Norway"), lwd = 2, col = "black", lty = c(1, 2, -1), cex = 1.5, y.intersp=0.5)
text(1980, 5, "First year of Gro Harlem's tenure", cex = 1.3)
text(1994, 12.8, bquote(paste("MSPE" == .(round(norw.ms.synth$loss.v, 5)))), cex = 1.5)
arrows(1985, 5, 1986, 5, length = .1)
rug(norw.ms$logsmilex_f2, side = 2, lwd=1, ticksize = .025)
#pdf("Norway_ms1_(609).pdf") 


### 4. Placebo cases ##
### Run placebo for each donor unit
for (j in 1:length(donor))
{
    norw.ms_pl.dp<- 
        dataprep( 
            foo = norw.ms, 
            predictors = predictors2,
            predictors.op = "mean", 
            dependent = "logsmilex_f2", 
            special.predictors = list(list("logsmilex_f2",1985, c("mean"))),
            unit.variable = "cown", 
            time.variable = "year",  
            treatment.identifier = donor[j], 
            controls.identifier = donor[-j], 
            time.predictors.prior = c(1976:1985), 
            time.optimize.ssr = c(1976:1985), 
            unit.names.variable = "country", 
            time.plot = c(1976:tr.end[i]) 
        ) 
    
    norw.ms_synth_pl <- synth(norw.ms_pl.dp)
    synthms_pl.tables <- synth.tab( dataprep.res = norw.ms_pl.dp, synth.res = norw.ms_synth_pl, round.digit = 5) 
    assign(paste("W",treated[i],donor[j], sep="."), (synthms_pl.tables$tab.w))
    assign(paste("V",treated[i],donor[j], sep="."), (norw.ms_synth_pl$solution.v))
    assign(paste("Pred",treated[i],donor[j], sep="."), (synthms_pl.tables$tab.pred))
    
    assign(paste("Y1",treated[i],donor[j], sep="."), (norw.ms_pl.dp$Y1))
    assign(paste("Y1s",treated[i],donor[j], sep="."), (norw.ms_pl.dp$Y0plot%*%norw.ms_synth_pl$solution.w))
    assign(paste("Alpha",treated[i],donor[j], sep="."), (norw.ms_pl.dp$Y1-norw.ms_pl.dp$Y0plot%*%norw.ms_synth_pl$solution.w))
}

### 5. Post Estimation results ###

## pre-RMSPE comparison ##
norw.prermspe <- sqrt(mean(norw.ms.alpha[1:10]^2, na.rm=TRUE))
norw.postrmspe <- sqrt(mean(norw.ms.alpha[11:21]^2, na.rm=TRUE))
norw.prermspe.d <- numeric(length(donor))
norw.postrmspe.d <- numeric(length(donor))

for (j in 1:length(donor))
{
    temp <- get(paste("Alpha", treated[i], donor[j], sep="."))
    norw.prermspe.d[j] <- sqrt(mean(temp[1:10]^2,na.rm=TRUE ))
    norw.postrmspe.d[j] <- sqrt(mean(temp[11:21]^2,na.rm=TRUE ))
}


# 5.1 e(avg_pre_rmspe_p) A measure of fit. Concerning if significant.
## The proportion of placebos that have a pre- treatment RMSPE at least as large as the average of the treated units. 
rmspe.fit <- length(which(norw.prermspe.d>norw.prermspe))/length(donor)

## 5.2  e(pval_joint_post_t)
## The proportion of placebos that have a ratio of post-treatment RMSPE over pre-treatment RMSPE at least as large as the average ratio for the treated units.
pval.joint.post_t <- length(which(norw.postrmspe.d/norw.prermspe.d>norw.postrmspe/norw.prermspe))/length(donor)

alphas <-norw.ms.alpha
for (j in 1: length(donor))
{
    temp1 <- get(paste("Alpha",treated[i],donor[j], sep="."))
    alphas <- cbind(alphas, temp1)
    colnames(alphas)[j+1] <- paste(treated[i],donor[j], sep=".")
}

## 5.3  e(pvals_t)
# A vector of the proportions of placebo pseudo t-statistics (unit’s effect divided by its pre-treatment RMSPE) that are at least as largeas the main pseudo t-statistic for each post-treatment period.
pvals_t <- numeric(11) 
for (t in 1:11)
{
    pvals_t[t] <-  length(which(abs(alphas[t+10,-1]/norw.prermspe.d)>abs(alphas[t+10,1]/norw.prermspe)))/length(donor)
}

plot(1:11, pvals_t, ylim=c(0,1), axes = FALSE,xlab = "", ylab = "",  pch=20, bg="grey", lwd=5, cex=2.5)
box(); grid()
axis(1, at = 1:11, labels = seq(1, 11, 1), cex.axis = 1.25)
axis(2, at = seq(0, 1, 0.1), cex.axis = 1, las = 2)
title(xlab = "Number of years after treatment", cex.lab = 1.5, line = 2.4)
title(ylab = "Lead specific P-Values", cex.lab = 2.2, line = 2.6)
text(pvals_t, labels = round(pvals_t, 2), pos = 3, cex = 1.8)
#pdf("Norway_ms1_pvalue(609).pdf")

### 6. Plotting placebo cases and Norway
plot(1976:1996, norw.ms.alpha, type = "n", ylim = c(-3, 3), lwd = 2, xlab = "", ylab = "", axes = FALSE) 
box(); grid()
for(j in 3: 29){
    lines(1976:1996, alphas[, j], type = "l", col = "gray80", lwd = 2)
} 
lines(1976:1996, alphas[,1 ], lwd = 2) #NORWAY
abline(v=1986, lwd = 1)
abline(h=0, lwd = 1)
axis(1, at = 1976:1996, labels = seq(1976, 1996, 1), cex.axis = 1.25)
axis(2, at = seq(-3, 3, 1), cex.axis = 1.25, las = 2)
title(xlab = "Year", cex.lab = 1.5, line = 2)
title(ylab = "Difference in Military Spending", cex.lab = 1.5, line = 2)
text(1980, -2.5, "Treatment time", cex = 1.2)
arrows(1984, -2.5, 1986, -2.5, length = .1)
legend(list(x = 1975,y = 3.3), ncol = 1, bty = "n", c("Norway", "Placebo Countries"), lty=c(1,1), lwd = 2, col = c("black", "gray80"), cex = 1.25, y.intersp=0.5)
#pdf("Norway_ms1pl_(609).pdf")

### 7. TABLES

synthms.tables$tab.w #(country weights)
synthms.tables$tab.pred #(Balance table)



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 




######### 2. MT with  Multiple Imputation DATA 1                  (Appendix D.4 )   #####
rm(list=ls())
library(mice)
library(VIM)
datak61 <- read_dta("~/Dropbox/Women Leaders and Conflict/Empirical Analysis/Rdata and scripts/Data_091419_k61.dta")
vars = c("cown", "year", "female_a", "country_pw", "smilex_f2", "polity2", "gdpgrowth_pw", "gdppc_pw", "women_legst", "cinc_c", "trade_cowt", "threat_f2", "allies", "issues", "rival")
data <- datak61[vars]


### 1. Checking the pattern of missing values
dataforplot <- data[,5:15]
colnames(dataforplot)[1:8] <- c("Military Spending", "Democracy", "GDP growth", "GDP pc", "Women in Parliment", "Capabilities", "Trade", "Threat")

mice_plot <- aggr(dataforplot, col=c('blue','yellow'),
                    numbers=TRUE, sortVars=F ,
                    labels=names(dataforplot), cex.axis=.5,horiz=TRUE,
                    gap=5, ylab=c("Missing data","Pattern"))

# Getting imputed data
imputed_Data <- mice(data, m=5, maxit = 50, method = c('pmm'), seed = 500)

# Extracting data
impdata1 <- complete(imputed_Data, 1)
impdata2 <- complete(imputed_Data, 2)
impdata3 <- complete(imputed_Data, 3)
impdata4 <- complete(imputed_Data, 4)
impdata5 <- complete(imputed_Data, 5)

### # MULTIPLE TREATMENT WITH FIRST IMPUTED DATA
colnames(impdata1)[3]=c("female")

### 2.  Cleaning up the data

# Sri lanka coded female, if the president is female 
impdata1$female[(impdata1$cown==780)] <- 0
impdata1$female[(impdata1$cown==780 & impdata1$year>1994 & impdata1$year<2005)] <- 1
# 1977 recoded as male=0, as Gandhi was in office till March 23#
impdata1$female[(impdata1$cown==750 & impdata1$year==1977)] <- 0
# India second treatment coded as zero #
impdata1$female[(impdata1$cown==750 & impdata1$year>1979 & impdata1$year<1985)] <- 0
#Meir only was in office till June 03, so I have recoded 1974 as having a male leader
impdata1$female[(impdata1$cown==666 & impdata1$year==1974)] <- 0
# canada only one year treatment coded as non-treatment
impdata1$female[(impdata1$cown==20 & impdata1$year==1993)] <- 0
# Ecuador only one year treatment coded as non-treatment
impdata1$female[(impdata1$cown==130 & impdata1$year==1997)] <- 0
# delete switzerland
impdata1<-impdata1[impdata1$cown!=225,]
# norway one year treatment coded as zero , 1981
impdata1$female[(impdata1$cown==385 & impdata1$year==1981)] <- 0
# Chile, the second term after 4 year break (2014 and 2015) is coded as zero
impdata1$female[(impdata1$cown==155 & impdata1$year>=2014)] <- 0
# Pakistan second treatment coded as zero
impdata1$female[(impdata1$cown==770 & impdata1$year>1992 & impdata1$year<1997)] <- 0
#Liberia First treatment coded as zero
impdata1$female[(impdata1$cown==450 & impdata1$year>1995 & impdata1$year<1998)] <- 0

### Duplicate countries with two treatments
# Original Argentina is cown 160, and a replica Argentina cown 161
df160 <- data.frame(matrix(ncol = dim(impdata1)[2], nrow = length(unique(impdata1$year))))
df160 <- impdata1[impdata1$cown==160,]
df160$cown=161
df160$country_pw <- "Arg-new"
df160$female[(df160$year>1973 & df160$year<1977)] <- 0
impdata1$female[(impdata1$cown==160 & impdata1$year>2006)] <- 0

# Append df160
impdata1 <- rbind(impdata1, df160)

# Original Norway is cown 385, and a replica Norway cown 386
df385 <- data.frame(matrix(ncol = dim(impdata1)[2], nrow = length(unique(impdata1$year))))
df385 <- impdata1[impdata1$cown==385,]
df385$cown=386
df385$country_pw <- "Norw-new"
df385$female[(df385$year>2012)] <- 0
impdata1$female[(impdata1$cown==385 & impdata1$year<2012)] <- 0

# Append df385
impdata1 <- rbind(impdata1, df385)

# Original Philippines is cown 840, and a replica Philippines is cown 841
df840 <- data.frame(matrix(ncol = dim(impdata1)[2], nrow = length(unique(impdata1$year))))
df840 <- impdata1[impdata1$cown==840,]
df840$cown=841
df840$country_pw <- "Philp-new"
df840$female[(df840$year>1985 & df840$year<1993)] <- 0
impdata1$female[(impdata1$cown==840 & impdata1$year>2000)] <- 0

# Append df840
impdata1 <- rbind(impdata1, df840)


### 3. Data set-up for synth ###
#### Identify treated countries ###
macro <- aggregate(female ~ cown, data = impdata1, sum)
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
    tr.years <- impdata1$year[(impdata1$female==1 & impdata1$cown==treated[i])]
    tr.start[i] <- min(tr.years)
    tr.end[i] <- max(tr.years)
}

impdata1 <- merge(impdata1, macro, by="cown")

predictors= c("polity2", "gdpgrowth_pw", "gdppc_pw", "women_legst", "cinc_c", "trade_cowt", "threat_f2", "allies", "issues", "rival")
impdata1$logsmilex <- log(impdata1$smilex_f2+0.01)


#### 4. Running Multiple Treatment SCM ###
tic()

for (i in 1:length(treated))
{
    dataprep.out<- 
        dataprep( 
            foo = impdata1, 
            predictors = predictors,
            predictors.op = "mean", 
            special.predictors = list(list("logsmilex",tr.start[i]-1, c("mean"))),
            dependent = "logsmilex",  
            unit.variable = "cown", 
            time.variable = "year",  
            treatment.identifier = treated[i], 
            controls.identifier = donor, 
            time.predictors.prior = c((max(1955,tr.start[i]-10)):tr.start[i]), 
            time.optimize.ssr = c((max(1955,tr.start[i]-10)):tr.start[i]), 
            unit.names.variable = "country_pw", 
            time.plot = c((max(1955,tr.start[i]-10)):tr.end[i]) 
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
    for (j in 28:length(donor))
    {    
        if (i==2 & j==29 | i==4 & j== 10 | i==4 & j== 28 | i==6 & j== 24| i==7 & j== 43 | i==10 & j== 2
            | i==10 & j== 16 | i==11 & j== 16 | i==14 & j== 10 | i==14 & j== 28 | i==21 & j== 20 | i==22 & j== 10 
            | i==26 & j== 29)  #not converging cases
        {
            assign(paste("Y1",treated[i],donor[j], sep="."), (NA))
            assign(paste("Y1s",treated[i],donor[j], sep="."), (NA))
            assign(paste("Alpha",treated[i],donor[j], sep="."), (NA))
        }
        else
        {
        dataprep.out<- 
            dataprep( 
                foo = impdata1, 
                predictors = predictors,
                special.predictors = list(list("logsmilex",tr.start[i]-1, c("mean"))),
                predictors.op = "mean", 
                dependent = "logsmilex",  
                unit.variable = "cown", 
                time.variable = "year",  
                treatment.identifier = donor[j], 
                controls.identifier = donor[-j], 
                time.predictors.prior = c((max(1955,tr.start[i]-10)):tr.start[i]), 
                time.optimize.ssr = c((max(1955,tr.start[i]-10)):tr.start[i]), 
                unit.names.variable = "country_pw", 
                time.plot = c((max(1955,tr.start[i]-10)):tr.end[i]) 
            )
        
        synth.out <- synth(dataprep.out)
        synth.tables <- synth.tab( dataprep.res = dataprep.out, synth.res = synth.out) 
        assign(paste("Y1",treated[i],donor[j], sep="."), (dataprep.out$Y1))
        assign(paste("Y1s",treated[i],donor[j], sep="."), (dataprep.out$Y0plot%*%synth.out$solution.w))
        assign(paste("Alpha",treated[i],donor[j], sep="."), (dataprep.out$Y1-dataprep.out$Y0plot%*%synth.out$solution.w))
        
    } #j loop
    
    } #NA loop
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

alpha <- merge(alpha93, alpha94, by="date", all=T)
alpha <- merge(alpha, alpha95, by="date", all=T)
alpha <- merge(alpha, alpha140, by="date", all=T)
alpha <- merge(alpha, alpha155, by="date", all=T)
alpha <- merge(alpha, alpha160, by="date", all=T)
alpha <- merge(alpha, alpha161, by="date", all=T)
alpha <- merge(alpha, alpha200, by="date", all=T)
alpha <- merge(alpha, alpha260, by="date", all=T)
alpha <- merge(alpha, alpha355, by="date", all=T)
alpha <- merge(alpha, alpha375, by="date", all=T)
alpha <- merge(alpha, alpha385, by="date", all=T)
alpha <- merge(alpha, alpha386, by="date", all=T)
alpha <- merge(alpha, alpha390, by="date", all=T)
alpha <- merge(alpha, alpha450, by="date", all=T)
alpha <- merge(alpha, alpha640, by="date", all=T)
alpha <- merge(alpha, alpha666, by="date", all=T)
alpha <- merge(alpha, alpha732, by="date", all=T)
alpha <- merge(alpha, alpha750, by="date", all=T)
alpha <- merge(alpha, alpha770, by="date", all=T)
alpha <- merge(alpha, alpha780, by="date", all=T)
alpha <- merge(alpha, alpha800, by="date", all=T)
alpha <- merge(alpha, alpha840, by="date", all=T)
alpha <- merge(alpha, alpha850, by="date", all=T)
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

## 6. Pre RMSPE comparison ##
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
            alpha93[t+10,2+temp[1]],
            alpha94[t+10,2+temp[2]],
            alpha95[t+10,2+temp[3]],
            alpha140[t+10,2+temp[4]],
            alpha155[t+10,2+temp[5]],
            alpha160[t+10,2+temp[6]],
            alpha161[t+10,2+temp[7]],
            alpha200[t+10,2+temp[8]],
            alpha260[t+10,2+temp[9]],
            alpha355[t+10,2+temp[10]],
            alpha375[t+10,2+temp[11]],
            alpha385[t+10,2+temp[12]],
            alpha386[t+10,2+temp[13]],
            alpha390[t+10,2+temp[14]],
            alpha450[t+10,2+temp[15]],
            alpha640[t+10,2+temp[16]],
            alpha666[t+10,2+temp[17]],
            alpha732[t+10,2+temp[18]],
            alpha750[t+10,2+temp[19]],
            alpha770[t+10,2+temp[20]],
            alpha780[t+10,2+temp[21]],
            alpha800[t+10,2+temp[22]],
            alpha840[t+10,2+temp[23]],
            alpha841[t+10,2+temp[24]],
            alpha850[t+10,2+temp[25]],
            alpha900[t+10,2+temp[26]],
            alpha920[t+10,2+temp[27]]), na.rm=TRUE)
        
    }
}

alpha.l <- alpha.treated[11:(10+lead), "ave.tr"]

for (t in 1:lead)
{
    alpha.l[t] <- mean( c(
        alpha93[t+10,2],
        alpha94[t+10,2],
        alpha95[t+10,2],
        alpha140[t+10,2],
        alpha155[t+10,2],
        alpha160[t+10,2],
        alpha161[t+10,2],
        alpha200[t+10,2],
        alpha260[t+10,2],
        alpha355[t+10,2],
        alpha375[t+10,2],
        alpha385[t+10,2],
        alpha386[t+10,2],
        alpha390[t+10,2],
        alpha450[t+10,2],
        alpha640[t+10,2],
        alpha666[t+10,2],
        alpha732[t+10,2],
        alpha750[t+10,2],
        alpha770[t+10,2],
        alpha780[t+10,2],
        alpha800[t+10,2],
        alpha840[t+10,2],
        alpha841[t+10,2],
        alpha850[t+10,2],
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
text(pvalue, labels = round(pvalue, 2), pos = 1, cex = 1.5, offset = 1)
title("Lead Specific Significance Level (P-Values) for Military Spending")
#mtms_imp1_pvalue.pdf


library(gdata)
keep(alpha, alpha.donor, alpha.treated, impdata2, impdata3, impdata4, impdata5, pvalue, postestimation, sure=TRUE)



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 


#### 3. Robustness test UK (leave one out)                        (Appendix D.5) ####
rm(list=ls())
setwd("~/Dropbox/Women Leaders and Conflict/Empirical Analysis/Rdata and scripts")
#### July 24th Robustness test UK (leave one out)


library(readstata13)
ukms1 <- read.dta13("UK_ms.dta")
library(Synth)
library(foreign)
library(tictoc)
colnames(ukms1)[4]=c("female")

# 1. Run the main model #

macro <- aggregate(female ~ cown, data = ukms1, sum)
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
    tr.years <- ukms1$year[(ukms1$female==1 & ukms1$cown==treated[i])]
    tr.start[i] <- min(tr.years)
    tr.end[i] <- max(tr.years)
}
ukms <- merge(ukms1, macro, by="cown")

predictors = c("polity2_20", "gdppc_pw", "gdpgrowth_pw","women_legst", "cinc_c", "logtrade", "logthreat", "allies", "issues", "rival")

ukms_all.dp <- dataprep(
    foo=ukms, 
    predictors = predictors,
    special.predictors = list(list("logsmilex_f2",1978, c("mean"))),
    predictors.op = "mean", 
    dependent = "logsmilex_f2", 
    unit.variable = "cown", 
    time.variable = "year",  
    treatment.identifier = 200, 
    controls.identifier = donor, 
    time.predictors.prior = c(1960:1978), 
    time.optimize.ssr = c(1960:1978), 
    unit.names.variable = "country", 
    time.plot = c(1960:1990) 
)
ukms_all.synth <- synth(ukms_all.dp, "ALL")

synthms.tables <- synth.tab(dataprep.res = ukms_all.dp, synth.res = ukms_all.synth) 
ukms_all.pred <- ukms_all.dp$Y0plot %*% ukms_all.synth$solution.w ####Y1s in Previous code
ukms_all.alpha<- ukms_all.dp$Y1plot - ukms_all.pred

# Check country weights 
a <- synthms.tables$tab.w

storegaps <- matrix(NA, length(1960:1990), 6)
mspe <- matrix(NA, length(1), 6 )

colnames(storegaps) <- c(2,20,211,220,380,900)
colnames(mspe) <- c(2,20,211,220,380,900)

weights <- data.frame(cown=macro[macro$donor==1,1])

for(k in 1:6){
    omit <- c(2,20,211,220,380,900)[k]
    ukms.dp <- dataprep(
        foo=ukms, 
        predictors = predictors,
        special.predictors = list(list("logsmilex_f2",1978, c("mean"))),
        predictors.op = "mean", 
        dependent = "logsmilex_f2", 
        unit.variable = "cown", 
        time.variable = "year",  
        treatment.identifier = 200, 
        controls.identifier = donor[-which(donor==omit)], 
        time.predictors.prior = c(1960:1978), 
        time.optimize.ssr = c(1960:1978), 
        unit.names.variable = "country", 
        time.plot = c(1960:1990) 
    )
    ukms.synth <- synth(ukms.dp, "ALL")
    storegaps[,k] <- (ukms.dp$Y0%*%ukms.synth$solution.w)
    mspe[,k] <- (ukms.synth$loss.v)
    synth_k.tables <- synth.tab(dataprep.res = ukms.dp, synth.res = ukms.synth) 
    n=which(weights$cown==omit)
    weights1 <- rep(NA, 33)
    weights1[1:(n-1)] <- synth_k.tables$tab.w[1:(n-1),1]
    weights1[(n+1):33] <- synth_k.tables$tab.w[n:32,1]
    weights <- cbind(weights, weights1) 
} # close loop over leave one outs
colnames(weights) <- c("cown", "2", "20", "211", "220", "380", "900")

# Plotting 
plot(1960:1990, ukms.dp$Y1plot, ylim = c(8, 13),  type = "n", lwd = 2, xlab = "", ylab = "", cex.axis=1.3)
box(); grid()
abline(v = 1979, lwd = 1)  
lines(1960:1990, ukms_all.dp$Y1plot, col="black", lty="solid")
lines(1960:1990, ukms_all.pred, col="black", lty="dashed",lwd=2)
for(i in 1:6){
    lines(1960:1990,storegaps[,i],col="blue",lty="solid")
}
title(xlab = "Year", cex.lab = 1.5, line = 2)
title(ylab = "Military Spending", cex.lab = 2.5, line = 2)
legend(list(x = 1959,y = 13.4), 
       ncol = 1, bty = "n", c("UK", "Synthetic UK", "Synthetic UK (leave-one-out)"), 
       lty=c("solid","dashed","solid"), lwd = 2, col=c("black","black","blue"), 
       cex = 1.3, y.intersp=0.7)

text(1969.3, 8.8, "First year of Thatcher's tenure", cex = 1.4)
arrows(1976, 8.8, 1979, 8.8, length = .15)
text(1985.2, 12.5, bquote(paste("Average MSPE" == .(round(rowMeans(mspe), 5)))), cex = 1.3)
rug(ukms1$logsmilex_f2, side = 2, lwd=1, ticksize = .025)







####### 4. Philippines Lame Duck 1995-2010                        (Appendix E) ######
#rm(list=ls())
library(Synth)
library(foreign)
library(tictoc)
library(mice)
library(VIM)
library(MSCMT)
library(ggplot2)
library(haven)

phil= read_dta("~/Desktop/Ulker desktop/MSCMT/Phil_2nd_mixed_large.dta")
vars = c("cown", "year", "female_a", "country_pw", "logsmilex_f2", "polity2_20", "gdpgrowth_pw", "gdppc_pw", "women_legst", "cinc_c", "logtrade", "logthreat", "allies", "issues", "rival", "treatment")
data <- phil[vars]

# 1. Getting imputed data
imputed_Data <- mice(data, m=5, maxit = 50, method = c('pmm'), seed = 500)
data1 <- complete(imputed_Data, 1)

data <- as.data.frame(data1)
impdata_long <- listFromLong(data, unit.variable="cown", time.variable="year", unit.names.variable="country_pw")

# 2.1 MSCMT analysis Single Treatment
treatment.identifier <- c("Philippines")
controls.identifier  <- setdiff(colnames(impdata_long[[1]]),
                                c(treatment.identifier, "Philippines"))

times.dep  <- cbind("logsmilex_f2" = c(1995,2003))
times.pred <- cbind("polity2_20"   = c(1995,2003),
                    "logsmilex_f2" = c(2003),
                    "gdppc_pw"     = c(1995,2003),
                    "gdpgrowth_pw" = c(1995,2003),
                    "women_legst"  = c(1995,2003),
                    "cinc_c"       = c(1995,2003),
                    "logtrade"     = c(1995,2003),
                    "logthreat"    = c(1995,2003),
                    "allies"       = c(1995,2003),
                    "issues"       = c(1995,2003),
                    "rival"        = c(1995,2003))

phil.fns <- rep("id", ncol(times.pred))


res_phil <- mscmt(impdata_long, treatment.identifier, controls.identifier, times.dep, times.pred,  phil.fns, seed=1)
res_phil

# 2.2 MSCMT analysis GAPS PLOT
ggplot(res_phil, type="gap")


#  3.1 MSCMT analysis Running placebo studies
res_phil_pl <- mscmt(impdata_long, treatment.identifier, controls.identifier, times.dep, times.pred, placebo=TRUE, phil.fns, seed=1)

#  3.2 MSCMT analysis P-value calculation

pvalue1 <- pvalue(res_phil_pl, exclude.ratio = 1, ratio.type = "mspe",  alternative="less")


# 4.1  MSCMT analysis PLOTS ##
library(ggplot2)
library(grid)
phil <- ggplot(res_phil, type="comparison", ylab="Military Spending (logged)", main="Comparison of Military Spending") +
    theme(legend.position = "none")+
    geom_vline(aes(xintercept = -400), colour="blue", linetype = "longdash") + 
    theme(legend.position = "none") + theme(axis.title.y = element_text (size=16))+
    theme(plot.title = element_text(hjust = 0.5, size=16))+
    scale_colour_manual(values=c("red","black"), labels=c("Philippines","Synthesized Philippines"), labs(fill=""))+
    theme(legend.position = c(0.24, 1.1)) + theme(legend.text = element_text(size=16))
phil

MSPE <- grobTree(textGrob("MSPE=0.0003", x=0.7,  y=0.1, hjust=0,
                          gp=gpar(col="black", fontsize=18)))
phil + annotation_custom(MSPE)

# 4.2  MSCMT analysis PLOTS ##
#Plotting lead-specific p-values with exclude ration 1
p1 <- as.data.frame(pvalue1)
year <- seq(2004, 2010, 1)
p1 <- cbind(year, p1)
colnames(p1)[2] <- c("P-value")

ggplot(p1, aes(year, `P-value`, label = round(`P-value`,3)))+ ggtitle("Lead Specific P-values") +
    theme(plot.title = element_text(hjust = 0.5, size=16)) +  theme(axis.title.y = element_text (size=20))+
    geom_point()+geom_text(nudge_y = 0.05, size = 5, fontface="bold") + xlim(2004, 2010) +ylim(0, 0.5)

#Phil_2nd_mixed_large_pvalue1.pdf

# 4.3  MSCMT analysis PLOTS ##
## Placebo Plots
ggplot(res_phil_pl, type="placebo.gaps", main= "Philippines and placebo cases", alternative = c("less"), ratio.type="rmspe", exclude.ratio = 1,  
       treatment.time=1969, ylab="Gaps for Philippines and donor counrties")+ 
    theme(plot.title = element_text(hjust = 0.8,   size = 24, face = "bold"))+
    theme(legend.text = element_text(size=16))+
    theme(legend.key.size = unit(1.2, "cm"))


## 5  MSCMT analysis Donor Weights and balance tables
balancetable <- res_phil$predictor.table
dweight <- res_phil$w



##### 5. ARGENTINA Lame Duck 2003 to 2015.                        (Appendix E) ####
rm(list=ls())
library(haven)
arg<- read_dta("~/Dropbox/Women Leaders and Conflict/Empirical Analysis/Rdata and scripts/Arg_lameduck.dta")
library(Synth)
library(foreign)
library(tictoc)
library(plm)
library(mice)
library(VIM)

### 1. Imputing missing values.
vars = c("cown", "countrycode_pw", "year", "female_a", "logsmilex_f2", "polity2_20", "gdpgrowth_pw", "gdppc_pw", "women_legst", "cinc_c", "logtrade", "allies", "treatment")
data_arg <- arg[vars]



### 1.1 Checking the pattern of missing values
dataforplot <- data_arg[,5:12]
colnames(dataforplot)[1:8] <- c("Military Spending", "Democracy", "GDP growth", "GDP pc", "Women in Parliment", "Capabilities", "Trade", "allies")

mice_plot <- aggr(dataforplot, col=c('blue','yellow'),
                  numbers=TRUE, sortVars=F ,
                  labels=names(dataforplot), cex.axis=.5,horiz=TRUE,
                  gap=5, ylab=c("Missing data","Pattern"))

# 1.2 Getting imputed data
imputed_Data <- mice(data_arg, m=5, maxit = 50, method = c('pmm'), seed = 500)
data1 <- complete(imputed_Data, 1)



data1 = as.data.frame(data1)
arg_long <- listFromLong(data1, unit.variable="cown", time.variable="year", unit.names.variable="countrycode_pw")

# 2.1 MSCMT analysis Single Treatment
treatment.identifier <- c("ARG")
controls.identifier  <- setdiff(colnames(arg_long[[1]]),
                                c(treatment.identifier, "ARG"))
times.dep  <- cbind("logsmilex_f2"  = c(2003,2011))
times.pred <- cbind("polity2_20"       = c(2003,2011),
                    "logsmilex_f2"  = c(2011),
                    "gdppc_pw"      = c(2003,2011),
                    "gdpgrowth_pw"  = c(2003,2011),
                    "women_legst"   = c(2003,2011),
                    "cinc_c"        = c(2003,2011),
                    "logtrade"      = c(2003,2011),
                    "allies"        = c(2003,2011))
#"issues"        = c(2003,2011),
#"rival"         = c(2003,2011))

agg.fns <- rep("id", ncol(times.pred))


res_arg <- mscmt(arg_long, treatment.identifier, controls.identifier, times.dep, times.pred,  agg.fns, seed=1)
res_arg

# 2.2 MSCMT analysis GAPS PLOT
#ggplot(res_arg, type="gap")


#  3.1 MSCMT analysis Running placebo studies
res_arg_pl <- mscmt(arg_long, treatment.identifier, controls.identifier, times.dep, times.pred, placebo=TRUE, agg.fns, seed=1)

#  3.2 MSCMT analysis P-value calculation
pvalue1 <- pvalue(res_arg_pl, exclude.ratio = 1, ratio.type = "mspe",  alternative="greater")




# 4.1  MSCMT analysis PLOTS ##
library(ggplot2)
library(grid)
arg_p <- ggplot(res_arg, type="comparison", ylab="Military Spending (logged)", main="Comparison of Military Spending") +
    #theme(legend.position = "none")+
    geom_vline(aes(xintercept = -400), colour="blue", linetype = "longdash") + 
    theme(legend.position = "none") + theme(axis.title.y = element_text (size=16))+
    theme(plot.title = element_text(hjust = 0.5, size=16))+
    scale_colour_manual(values=c("red","black"), labels=c("Argentina","Synthesized Argentina"), labs(fill=""))+
    theme(legend.position = c(0.8, -0.01)) + theme(legend.text = element_text(size=16))+
    theme(legend.key.size = unit(1.2, "cm"))

arg_p
MSPE <- grobTree(textGrob("MSPE=0.012", x=0.1,  y=0.9, hjust=0,
                          gp=gpar(col="black", fontsize=18)))
arg_p + annotation_custom(MSPE) 


# 4.2  MSCMT analysis PLOTS ##
#Plotting lead-specific p-values with exclude ration 1
p1 <- as.data.frame(pvalue1)
year <- seq(2012, 2015, 1)
p1 <- cbind(year, p1)
colnames(p1)[2] <- c("P-value")

ggplot(p1, aes(year, `P-value`, label = round(`P-value`,3)))+ ggtitle("Lead Specific P-values") +
    theme(plot.title = element_text(hjust = 0.5, size=16)) +  theme(axis.title.y = element_text (size=20))+
    geom_point()+geom_text(nudge_y = 0.05, size = 5, fontface="bold") + xlim(2012, 2015) +ylim(0, 0.5)


# 4.3  MSCMT analysis PLOTS ##
## Placebo Plots
ggplot(res_arg_pl, type="placebo.gaps", main= "Argentina and placebo cases", alternative = c("greater"), ratio.type="rmspe", exclude.ratio = 1,  
       treatment.time=2011, ylab="Gaps for Argentina and donor counrties")+ 
    #theme(legend.position = "none")+
    theme(plot.title = element_text(hjust = 0.8,   size = 24, face = "bold"))+
    theme(legend.text = element_text(size=12))+
    theme(legend.key.size = unit(0.6, "cm"))


## 5  MSCMT analysis Donor Weights and balance tables
balancetable <- res_arg$predictor.table
dweight <- res_arg$w

#### 6. Argentina Health Expenditure with left countries only     (Appendix J.1 ) ##### 
rm(list=ls())
library(MSCMT)
arg_hl<- read_dta("~/Dropbox/Women Leaders and Conflict/Empirical Analysis/Rdata and scripts/arg_h_left (304).dta")

vars = c("cown", "year", "female_a", "country_pw", "countrycode_pw", "oecd_full", "health_gdp", "gov1rlc_PI", "women_legst", "polity2_20", "gdppc_wbn", "gdp_wbn")
vars2 = c("cown", "year", "country_pw", "health_gdp", "gov1rlc_PI", "women_legst", "polity2_20", "gdppc_wbn", "loggdp_wbn")

arg_hl <- arg_hl[vars2]
arg_hl<- as.data.frame(arg_hl)
arg_hl <- arg_hl[arg_hl$cown!=701, ]
arg_hllong <- listFromLong(arg_hl, unit.variable="cown", time.variable="year", unit.names.variable="country_pw")


treatment.identifier <- c("Argentina")
controls.identifier  <- setdiff(colnames(arg_hllong[[1]]),
                                c(treatment.identifier, "Argentina"))
times.dep  <- cbind("health_gdp"        = c(2000, 2007))
times.pred <- cbind("polity2_20"          = c(2000, 2007),
                    "health_gdp"        = c(2007, 2007),
                    "health_gdp"        = c(2004, 2004),
                    "health_gdp"        = c(2000, 2000),
                    "gdppc_wbn"         = c(2000, 2007),
                    "loggdp_wbn"         = c(2000, 2007),
                    "women_legst"         = c(2000, 2007))
agg.fns <- rep("mean", ncol(times.pred))                      
agg.fns <- rep("id", ncol(times.pred))                      

res <- mscmt(arg_hllong, treatment.identifier, controls.identifier, times.dep, times.pred, agg.fns, seed=1)
res

## Donor Weights and balance tables
balancetable <- res$predictor.table
dweight <- res$w

## Summary Statictics by Group
library(dplyr)

stbycown <- arg_hl %>% 
    group_by(cown) %>% 
    summarize(mean = mean(health_gdp),
              min = min(health_gdp),
              max = max(health_gdp),
              sd = sd(health_gdp))


res_pl <- mscmt(arg_hllong, treatment.identifier, controls.identifier, times.dep, times.pred, placebo=TRUE,  agg.fns, seed=1)
res_pl

ggplot(res_pl, type="p.value", draw.points = TRUE, legend = TRUE)
did(res_pl, range.pre=c(2000, 2007), range.post=c(2008,2015), alternative="greater")
ggplot(res, type="", draw.points = TRUE, legend = TRUE)

ggplot(res, type="comparison", ylab="Health Spending", main="") + geom_vline(aes(xintercept = -380), colour="blue", linetype = "longdash")

pl <- ggplot(res_pl, type="placebo.gaps", main= "Argentina and placebo cases", alternative = c("greater"), ratio.type="rmspe", include.mean = TRUE,  
             treatment.time=2008, ylab="Gaps for Argentina and donor counrties")+ theme(plot.title = element_text(hjust = 0.5)) 





