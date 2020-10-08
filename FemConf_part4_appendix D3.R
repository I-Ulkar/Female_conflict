####         ULKAR IMAMVERDI AND PATRICK SHEA       ###
###          Female Leaders and Foreign Policy      ###
###        Female Empowerment and MID x3 analysis   ###
###        APPENDIX D.3.1, D.3.2  D.3.3 and D.3.4   ###
#######################################################




######### 1.  Interaction GDP growth with GDP per capita (Appendix D.3.1 ) #####
#### ####ISRAEL  Military Spending R1 #### ####
rm(list=ls())
library(readstata13)
setwd("~/Dropbox/Women Leaders and Conflict/Empirical Analysis/Rdata and scripts")
isrms1 <- read.dta13("isr53-73fem_ms504r.dta")
library(Synth)
library(foreign)
library(tictoc)
library(plm)
colnames(isrms1)[4]=c("female")

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

## Creating interaction term
isrms1$gdp_inter <- isrms1$gdpgrowth_pw*isrms1$gdppc_pw

### 2. Estimating Synthetic Counterfactual for Israel 

predictors = c( "polity2_20", "gdpgrowth_pw", "gdppc_pw", "gdp_inter", "women_legst", "cinc_c", "logtrade", "logthreat", "allies", "issues", "rival")
isrms_r1.dp <- dataprep(
    foo=isrms1, 
    predictors = predictors,
    special.predictors = list(list("logsmilex_f2", 1968, c("mean"))),
    predictors.op = "mean", 
    dependent = "logsmilex_f2", 
    unit.variable = "cown", 
    time.variable = "year",  
    treatment.identifier = 666, 
    controls.identifier = donor, 
    time.predictors.prior = c(1956:1968), 
    time.optimize.ssr = c(1956:1968), 
    unit.names.variable = "country", 
    time.plot = c(1956:tr.end[i]) 
)
isrms_r1.synth <- synth(isrms_r1.dp)


synthms.tables <- synth.tab( dataprep.res = isrms_r1.dp, synth.res = isrms_r1.synth) 
isrms_r1.pred <- isrms_r1.dp$Y0plot %*% isrms_r1.synth$solution.w ####Y1s in Previous code
isrms_r1.alpha<- isrms_r1.dp$Y1plot - isrms_r1.pred

## 3. Plot results
plot(1956:1973, isrms_r1.dp$Y1plot, ylim = c(3, 13),  type = "n", lwd = 2, xlab = "", ylab = "")
box(); grid()
abline(v = 1969, lwd = 1)  
lines(1956:1973, isrms_r1.dp$Y1plot, lwd = 2)
lines(1956:1973, isrms_r1.pred, lwd = 2, lty = 2)
title(xlab = "Year", cex.lab = 1.5, line = 2)
title(ylab = "Military Spending", cex.lab = 1, line = 2)
legend(list(x = 1956,y = 14.5), bty = "n", c("Israel", "Synthetic Israel"), lwd = 2, col = "black", lty = c(1, 2, -1), cex = 1.8, y.intersp=0.3)
text(1962.8, 5, "First year of Meir's tenure", cex = 2.3)
text(1971.3, 12.5, bquote(paste("MSPE" == .(round(isrms_r1.synth$loss.v, 5)))), cex = 1.9)
arrows(1967, 5, 1969, 5, length = .1)
rug(isrms1$logsmilex_f2, side = 2, lwd=1, ticksize = .025)
#isr_ms_r1.pdf


### 4. Placebo cases ##
### Run placebo for each donor unit
for (j in 1:length(donor))
{
    isrms_r1_pl.dp<- 
        dataprep( 
            foo = isrms1, 
            predictors = predictors,
            predictors.op = "mean",
            special.predictors = list(list("logsmilex_f2", 1960:1968, c("mean"))),
            dependent = "logsmilex_f2", 
            unit.variable = "cown", 
            time.variable = "year",  
            treatment.identifier = donor[j], 
            controls.identifier = donor[-j], 
            time.predictors.prior = c(1956:1968), 
            time.optimize.ssr = c(1956:1968), 
            unit.names.variable = "country", 
            time.plot = c(1956:tr.end[i]) 
        ) 
    
    isrms_r1_synth_pl <- synth(isrms_r1_pl.dp)
    synthms2_pl.tables <- synth.tab( dataprep.res = isrms_r1_pl.dp, synth.res = isrms_r1_synth_pl, round.digit = 5) 
    assign(paste("W",treated[i],donor[j], sep="."), (synthms2_pl.tables$tab.w))
    assign(paste("V",treated[i],donor[j], sep="."), (isrms_r1_synth_pl$solution.v))
    assign(paste("Pred",treated[i],donor[j], sep="."), (synthms2_pl.tables$tab.pred))
    
    assign(paste("Y1",treated[i],donor[j], sep="."), (isrms_r1_pl.dp$Y1))
    assign(paste("Y1s",treated[i],donor[j], sep="."), (isrms_r1_pl.dp$Y0plot%*%isrms_r1_synth_pl$solution.w))
    assign(paste("Alpha",treated[i],donor[j], sep="."), (isrms_r1_pl.dp$Y1-isrms_r1_pl.dp$Y0plot%*%isrms_r1_synth_pl$solution.w))
}


## pre-RMSPE comparison ##
isr.prermspe2 <- sqrt(mean(isrms_r1.alpha[1:13]^2, na.rm=TRUE))
isr.postrmspe2 <- sqrt(mean(isrms_r1.alpha[14:18]^2, na.rm=TRUE))
isr.prermspe.d2 <- numeric(length(donor))
isr.postrmspe.d2 <- numeric(length(donor))

for (j in 1:length(donor))
{
    temp <- get(paste("Alpha", treated[i], donor[j], sep="."))
    isr.prermspe.d2[j] <- sqrt(mean(temp[1:13]^2,na.rm=TRUE ))
    isr.postrmspe.d2[j] <- sqrt(mean(temp[14:18]^2,na.rm=TRUE ))
}

# 5.1 e(avg_pre_rmspe_p) A measure of fit. Concerning if significant.
## The proportion of placebos that have a pre- treatment RMSPE at least as large as the average of the treated units. 
rmspe.fit2<- length(which(isr.prermspe.d2>isr.prermspe2))/length(donor)


## 5.2  e(pval_joint_post_t)
## The proportion of placebos that have a ratio of post-treatment RMSPE over pre-treatment RMSPE at least as large as the average ratio for the treated units.
pval.joint.post_t2 <- length(which(isr.postrmspe.d2/isr.prermspe.d2>isr.postrmspe2/isr.prermspe2))/length(donor)

alphas2 <-isrms_r1.alpha
for (j in 1: length(donor))
{
    temp12 <- get(paste("Alpha",treated[i],donor[j], sep="."))
    alphas2 <- cbind(alphas2, temp12)
    colnames(alphas2)[j+1] <- paste(treated[i],donor[j], sep=".")
}

## 5.3  e(pvals_t)
# A vector of the proportions of placebo pseudo t-statistics (unit’s effect divided by its pre-treatment RMSPE) that are at least as largeas the main pseudo t-statistic for each post-treatment period.
pvals_t2 <- numeric(5) 
for (t in 1:5)
{
    pvals_t2[t] <-  length(which(abs(alphas2[t+13,-1]/isr.prermspe.d2)>abs(alphas2[t+13,1]/isr.prermspe2)))/length(donor)
}

## 5.3.2 Plotting p values
plot(1:5, pvals_t2, ylim=c(0,1), axes = FALSE,xlab = "", ylab = "",  pch=20, bg="grey", lwd=5, cex=2)
box(); grid()
axis(1, at = 1:5, labels = seq(1, 5, 1), cex.axis = 1.25)
axis(2, at = seq(0, 1, 0.1), cex.axis = 1, las = 2)
title(xlab = "Number of years after treatment", cex.lab = 1.5, line = 2.4)
title(ylab = "Lead specific P-Values", cex.lab = 2, line = 2.6)




plot(1956:1973, isrms_r1.alpha, type = "n", ylim = c(-3, 3), lwd = 2, xlab = "", ylab = "", axes = FALSE) 
box(); grid()
#for(j in 1: 24){
#lines(1956:1973, alphas2[, j], type = "l", col = "gray80", lwd = 2)
#} 

#lines(1953:1973, alphas[,2 ], type = "l", col = "gray80", lwd = 2) #bad fit (2=US)
lines(1956:1973, alphas2[,3 ],type = "l", col = "gray80", lwd = 2)
lines(1956:1973, alphas2[,4 ],type = "l", col = "gray80", lwd = 2)
lines(1956:1973, alphas2[,5 ], type = "l", col = "gray80", lwd = 2)
#lines(1956:1973, alphas2[,6 ],type = "l", col = "gray80", lwd = 2) #bad fit (94=)
lines(1956:1973, alphas2[,7 ],type = "l", col = "gray80", lwd = 2)
lines(1956:1973, alphas2[,8 ],type = "l", col = "gray80", lwd = 2)
lines(1956:1973, alphas2[,9 ],type = "l", col = "gray80", lwd = 2)
lines(1956:1973, alphas2[,10 ],type = "l", col = "gray80", lwd = 2)
lines(1956:1973, alphas2[,11 ],type = "l", col = "gray80", lwd = 2)
lines(1956:1973, alphas2[,12 ],type = "l", col = "gray80", lwd = 2)
lines(1956:1973, alphas2[,13 ],type = "l", col = "gray80", lwd = 2)
lines(1956:1973, alphas2[,14 ],type = "l", col = "gray80", lwd = 2)
lines(1956:1973, alphas2[,15 ],type = "l", col = "gray80", lwd = 2)
lines(1956:1973, alphas2[,16 ],type = "l", col = "gray80", lwd = 2)
lines(1956:1973, alphas2[,17 ],type = "l", col = "gray80", lwd = 2)
lines(1956:1973, alphas2[,18 ],type = "l", col = "gray80", lwd = 2)
lines(1956:1973, alphas2[,19 ],type = "l", col = "gray80", lwd = 2)
lines(1956:1973, alphas2[,20 ],type = "l", col = "gray80", lwd = 2)
lines(1956:1973, alphas2[,21 ],type = "l", col = "gray80", lwd = 2)
lines(1956:1973, alphas2[,22 ],type = "l", col = "gray80", lwd = 2)
lines(1956:1973, alphas2[,23 ],type = "l", col = "gray80", lwd = 2)
lines(1956:1973, alphas2[,24 ],type = "l", col = "gray80", lwd = 2)

lines(1956:1973, alphas2[,1 ], lwd = 2) #Israel 
abline(v=1969, lwd = 1)
abline(h=0, lwd = 1)

axis(1, at = 1956:1973, labels = seq(1956, 1973, 1), cex.axis = 1.25)
axis(2, at = seq(-3, 3, 1), cex.axis = 1.25, las = 2)
title(xlab = "Year", cex.lab = 1.5, line = 2)
title(ylab = "Difference in Military Spending", cex.lab = 1.5, line = 2)
text(1964.5, -2.5, "Treatment time", cex = 1.8)
arrows(1967, -2.5, 1969, -2.5, length = .1)
legend(list(x = 1956,y = 3.5), ncol = 1, bty = "n", c("Israel", "Placebo Countries"), lty=c(1,1), lwd = 2, col = c("black", "gray80"), cex = 1.25, y.intersp=0.5)







#### #### UNITED KINGDOM   Military Spending R1 ####
rm(list=ls())
setwd("~/Dropbox/Women Leaders and Conflict/Empirical Analysis/Rdata and scripts")

library(readstata13)
ukms1 <- read.dta13("UK_ms.dta")

library(Synth)
library(foreign)
library(tictoc)

colnames(ukms1)[4]=c("female")

### 1. Prepping data
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
ukms1 <- merge(ukms1, macro, by="cown")

## Creating interaction term
ukms1$gdp_inter <- ukms1$gdpgrowth_pw*ukms1$gdppc_pw
predictors = c("polity2_20", "gdppc_pw", "gdpgrowth_pw", "gdp_inter", "women_legst", "cinc_c", "logtrade", "logthreat", "allies", "issues", "rival")

### 2. Estimating Synthetic Counterfactual for UK 
ukms_r1.dp <- dataprep(
    foo=ukms1, 
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
    time.plot = c(1960:tr.end[i]) 
)
ukms_r1.synth <- synth(ukms_r1.dp, "ALL")

synthms.tables <- synth.tab( dataprep.res = ukms_r1.dp, synth.res = ukms_r1.synth) 
ukms_r1.pred <- ukms_r1.dp$Y0plot %*% ukms_r1.synth$solution.w ####Y1s in Previous code
ukms_r1.alpha<- ukms_r1.dp$Y1plot - ukms_r1.pred

## PLOTTING SINGLE TREATMENT GRAPH
plot(1960:1990, ukms_r1.dp$Y1plot, ylim = c(3, 13),  type = "n", lwd = 2, xlab = "", ylab = "")
box(); grid()
abline(v = 1979, lwd = 1)  
lines(1960:1990, ukms_r1.dp$Y1plot, lwd = 2)
lines(1960:1990, ukms_r1.pred, lwd = 2, lty = 2)
title(xlab = "Year", cex.lab = 1.5, line = 2)
title(ylab = "Military Spending", cex.lab = 2.5, line = 2)
legend(list(x = 1960,y = 14), bty = "n", c("UK", "Synthetic UK"), lwd = 2, col = "black", lty = c(1, 2, -1), cex = 1.3, y.intersp=0.3)
text(1969, 6, "First year of Thatcher's tenure", cex = 1.7)
arrows(1976, 6, 1979, 6, length = .2)
text(1985, 12.5, bquote(paste("MSPE" == .(round(ukms_r1.synth$loss.v, 5)))), cex = 1.5)
rug(ukms1$logsmilex_f2, side = 2, lwd=1, ticksize = .025)

### 4. Placebo cases ##
### Run placebo for each donor unit
for (j in 1:length(donor))
{
    ukms_r1_pl.dp<- 
        dataprep( 
            foo = ukms1, 
            predictors = predictors,
            predictors.op = "mean",
            special.predictors = list(list("logsmilex_f2",1978, c("mean"))),
            dependent = "logsmilex_f2", 
            unit.variable = "cown", 
            time.variable = "year",  
            treatment.identifier = donor[j], 
            controls.identifier = donor[-j], 
            time.predictors.prior = c(1960:1978), 
            time.optimize.ssr = c(1960:1978), 
            unit.names.variable = "country", 
            time.plot = c(1960:tr.end[i]) 
        ) 
    
    ukms_r1_synth_pl <- synth(ukms_r1_pl.dp)
    synthms_pl.tables <- synth.tab( dataprep.res = ukms_r1_pl.dp, synth.res = ukms_r1_synth_pl, round.digit = 5) 
    assign(paste("W",treated[i],donor[j], sep="."), (synthms_pl.tables$tab.w))
    assign(paste("V",treated[i],donor[j], sep="."), (ukms_r1_synth_pl$solution.v))
    assign(paste("Pred",treated[i],donor[j], sep="."), (synthms_pl.tables$tab.pred))
    
    assign(paste("Y1",treated[i],donor[j], sep="."), (ukms_r1_pl.dp$Y1))
    assign(paste("Y1s",treated[i],donor[j], sep="."), (ukms_r1_pl.dp$Y0plot%*%ukms_r1_synth_pl$solution.w))
    assign(paste("Alpha",treated[i],donor[j], sep="."), (ukms_r1_pl.dp$Y1-ukms_r1_pl.dp$Y0plot%*%ukms_r1_synth_pl$solution.w))
}

### 5. Post Estimation results ###
## pre-RMSPE comparison ##
uk.prermspe <- sqrt(mean(ukms_r1.alpha[1:19]^2, na.rm=TRUE))
uk.postrmspe <- sqrt(mean(ukms_r1.alpha[20:31]^2, na.rm=TRUE))
uk.prermspe.d <- numeric(length(donor))
uk.postrmspe.d <- numeric(length(donor))

for (j in 1:length(donor))
{
    temp <- get(paste("Alpha", treated[i], donor[j], sep="."))
    uk.prermspe.d[j] <- sqrt(mean(temp[1:19]^2,na.rm=TRUE ))
    uk.postrmspe.d[j] <- sqrt(mean(temp[20:31]^2,na.rm=TRUE ))
}

# 5.1 e(avg_pre_rmspe_p) A measure of fit. Concerning if significant.
## The proportion of placebos that have a pre- treatment RMSPE at least as large as the average of the treated units. 
rmspe.fit <- length(which(uk.prermspe.d>uk.prermspe))/length(donor)


## 5.2  e(pval_joint_post_t)
## The proportion of placebos that have a ratio of post-treatment RMSPE over pre-treatment RMSPE at least as large as the average ratio for the treated units.
pval.joint.post_t <- length(which(uk.postrmspe.d/uk.prermspe.d>uk.postrmspe/uk.prermspe))/length(donor)

alphas <-ukms_r1.alpha
for (j in 1: length(donor))
{
    temp1 <- get(paste("Alpha",treated[i],donor[j], sep="."))
    alphas <- cbind(alphas, temp1)
    colnames(alphas)[j+1] <- paste(treated[i],donor[j], sep=".")
}

## 5.3.1  e(pvals_t)
# A vector of the proportions of placebo pseudo t-statistics (unit’s effect divided by its pre-treatment RMSPE) that are at least as largeas the main pseudo t-statistic for each post-treatment period.
pvals_t <- numeric(12) 
for (t in 1:12)
{
    pvals_t[t] <-  length(which(abs(alphas[t+19,-1]/uk.prermspe.d)>abs(alphas[t+19,1]/uk.prermspe)))/length(donor)
}

## 5.3.2 Plotting p values
plot(1:12, pvals_t, ylim=c(0,1), axes = FALSE,xlab = "", ylab = "",  pch=20, bg="grey", lwd=5, cex=2.5)
box(); grid()
axis(1, at = 1:10, labels = seq(1, 10, 1), cex.axis = 1.25)
axis(2, at = seq(0, 1, 0.1), cex.axis = 1, las = 2)
title(xlab = "Number of years after treatment", cex.lab = 1.5, line = 2.4)
title(ylab = "Lead specific P-Values", cex.lab = 1.5, line = 2.6)
#uk_ms_r1_pvalue.pdf

### 6. Plotting placebo cases and UK
plot(1960:1990, ukms_r1.alpha, type = "n", ylim = c(-3, 3), lwd = 2, xlab = "", ylab = "", axes = FALSE) 
box(); grid()
for(j in 3: 33){
    lines(1960:1990, alphas[, j], type = "l", col = "gray80", lwd = 2)
} 
lines(1960:1990, alphas[,1 ], lwd = 2) #UK
abline(v=1979, lwd = 1)
abline(h=0, lwd = 1)
axis(1, at = 1960:1990, labels = seq(1960, 1990, 1), cex.axis = 1.25)
axis(2, at = seq(-3, 3, 1), cex.axis = 1.25, las = 2)
title(xlab = "Year", cex.lab = 1.5, line = 2)
title(ylab = "Difference in Military Spending", cex.lab = 1.5, line = 2)
text(1971, -2.5, "Treatment time", cex = 1.5)
arrows(1976, -2.5, 1979, -2.5, length = .1)
legend(list(x = 1958,y = 3.1), ncol = 1, bty = "n", c("UK", "Placebo Countries"), lty=c(1,1), lwd = 2, col = c("black", "gray80"), cex = 1.5, y.intersp=0.3)

#country weights Table B10
synthms.tables$tab.w
#Balance Table B4
synthms.tables$tab.pred




#### #### INDIA   Military Spending R1 ####
rm(list=ls())
setwd("~/Dropbox/Women Leaders and Conflict/Empirical Analysis/Rdata and scripts")

library(readstata13)
indms1 <- read.dta13("India_ms56-84_725.dta")
indms1 <- indms1[indms1$year<1977,]

library(Synth)
library(foreign)
library(tictoc)

colnames(indms1)[4]=c("female")

### 1. Prepping data
macro <- aggregate(female ~ cown, data = indms1, sum)
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
    tr.years <- indms1$year[(indms1$female==1 & indms1$cown==treated[i])]
    tr.start[i] <- min(tr.years)
    tr.end[i] <- max(tr.years)
}
indms1 <- merge(indms1, macro, by="cown")
indms1$gdp_inter <- indms1$gdpgrowth_pw*indms1$gdppc_pw
predictors = c("polity2_20", "gdppc_pw", "gdpgrowth_pw", "gdp_inter", "women_legst", "cinc_c", "logtrade", "logthreat_f2", "allies", "issues", "rival")

### 2. Estimating Synthetic Counterfactual for INDIA
indms_r1.dp <- dataprep(
    foo=indms1, 
    predictors = predictors,
    special.predictors = list(list("logsmilex_f2",1965, c("mean"))),
    predictors.op = "mean", 
    dependent = "logsmilex_f2", 
    unit.variable = "cown", 
    time.variable = "year",  
    treatment.identifier = 750, 
    controls.identifier = donor, 
    time.predictors.prior = c(1956:1965), 
    time.optimize.ssr = c(1956:1965), 
    unit.names.variable = "country", 
    time.plot = c(1956:tr.end[i]) 
)
indms_r1.synth <- synth(indms_r1.dp, "ALL")

synthms.tables <- synth.tab( dataprep.res = indms_r1.dp, synth.res = indms_r1.synth) 
indms_r1.pred <- indms_r1.dp$Y0plot %*% indms_r1.synth$solution.w ####Y1s in Previous code
indms_r1.alpha<- indms_r1.dp$Y1plot - indms_r1.pred

###. 3 Plotting the results INDIA
plot(1956:1976, indms_r1.dp$Y1plot, ylim = c(3, 13),  type = "n", lwd = 2, xlab = "", ylab = "")
box(); grid()
abline(v = 1966, lwd = 1)  
lines(1956:1976, indms_r1.dp$Y1plot, lwd = 2)
lines(1956:1976, indms_r1.pred, lwd = 2, lty = 2)
title(xlab = "Year", cex.lab = 1.5, line = 2)
title(ylab = "Military Spending", cex.lab = 2.5, line = 2)
legend(list(x = 1956,y = 13), bty = "n", c("India", "Synthetic India"), lwd = 2, col = "black", lty = c(1, 2, -1), cex = 1.3, y.intersp=0.3)
text(1960.5, 5, "First year of Gandhi's tenure", cex = 1.3)
arrows(1964.8, 5, 1966, 5, length = .2)
text(1973, 12.5, bquote(paste("MSPE" == .(round(indms_r1.synth$loss.v, 5)))), cex = 1.5)
rug(indms1$logsmilex_f2, side = 2, lwd=1, ticksize = .025)

### 4. Placebo cases ##
### Run placebo for each donor unit

for (j in 1:length(donor))
{
    indms_r1_pl.dp<- 
        dataprep( 
            foo = indms1, 
            predictors = predictors,
            predictors.op = "mean",
            special.predictors = list(list("logsmilex_f2",1965, c("mean"))),
            dependent = "logsmilex_f2", 
            unit.variable = "cown", 
            time.variable = "year",  
            treatment.identifier = donor[j], 
            controls.identifier = donor[-j], 
            time.predictors.prior = c(1956:1965), 
            time.optimize.ssr = c(1956:1965), 
            unit.names.variable = "country", 
            time.plot = c(1956:tr.end[i]) 
        ) 
    
    indms_r1_synth_pl <- synth(indms_r1_pl.dp)
    synthms_pl.tables <- synth.tab( dataprep.res = indms_r1_pl.dp, synth.res = indms_r1_synth_pl, round.digit = 5) 
    assign(paste("W",treated[i],donor[j], sep="."), (synthms_pl.tables$tab.w))
    assign(paste("V",treated[i],donor[j], sep="."), (indms_r1_synth_pl$solution.v))
    assign(paste("Pred",treated[i],donor[j], sep="."), (synthms_pl.tables$tab.pred))
    
    assign(paste("Y1",treated[i],donor[j], sep="."), (indms_r1_pl.dp$Y1))
    assign(paste("Y1s",treated[i],donor[j], sep="."), (indms_r1_pl.dp$Y0plot%*%indms_r1_synth_pl$solution.w))
    assign(paste("Alpha",treated[i],donor[j], sep="."), (indms_r1_pl.dp$Y1-indms_r1_pl.dp$Y0plot%*%indms_r1_synth_pl$solution.w))
    #bad.fit <- numeric(33)
    #bad.fit[j] <- as.numeric(ifelse(indms_r1_synth_pl$loss.v > indms_r1.synth$loss.v*5, 1, 0))
    
}

### 5. Post Estimation results ###
## pre-RMSPE comparison ##
ind.prermspe <- sqrt(mean(indms_r1.alpha[1:10]^2, na.rm=TRUE))
ind.postrmspe <- sqrt(mean(indms_r1.alpha[11:21]^2, na.rm=TRUE))
ind.prermspe.d <- numeric(length(donor))
ind.postrmspe.d <- numeric(length(donor))

for (j in 1:length(donor))
{
    temp <- get(paste("Alpha", treated[i], donor[j], sep="."))
    ind.prermspe.d[j] <- sqrt(mean(temp[1:10]^2,na.rm=TRUE ))
    ind.postrmspe.d[j] <- sqrt(mean(temp[11:21]^2,na.rm=TRUE ))
}

# 5.1 e(avg_pre_rmspe_p) A measure of fit. Concerning if significant.
## The proportion of placebos that have a pre- treatment RMSPE at least as large as the average of the treated units. 
rmspe.fit <- length(which(ind.prermspe.d>ind.prermspe))/length(donor)


## 5.2  e(pval_joint_post_t)
## The proportion of placebos that have a ratio of post-treatment RMSPE over pre-treatment RMSPE at least as large as the average ratio for the treated units.
pval.joint.post_t <- length(which(ind.postrmspe.d/ind.prermspe.d>ind.postrmspe/ind.prermspe))/length(donor)

alphas <-indms_r1.alpha
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
    pvals_t[t] <-  length(which(abs(alphas[t+10,-1]/ind.prermspe.d)>abs(alphas[t+10,1]/ind.prermspe)))/length(donor)
}

plot(1:11, pvals_t, ylim=c(0,1), axes = FALSE,xlab = "", ylab = "",  pch=20, bg="grey", lwd=5, cex=2.5)
box(); grid()
axis(1, at = 1:11, labels = seq(1, 11, 1), cex.axis = 1.25)
axis(2, at = seq(0, 1, 0.1), cex.axis = 1, las = 2)
title(xlab = "Number of years after treatment", cex.lab = 1.5, line = 2.4)
title(ylab = "Lead specific P-Values", cex.lab = 2.4, line = 2.6)
#ind_ms_r1_pvalue.pdf

### 6. Plotting placebo cases and India
plot(1956:1976, indms_r1.alpha, type = "n", ylim = c(-3, 3), lwd = 2, xlab = "", ylab = "", axes = FALSE) 
box(); grid()

for(j in 3: 24)
{lines(1956:1976, alphas[, j], type = "l", col = "gray80", lwd = 2)}
lines(1956:1976, alphas[,1 ], lwd = 2) #INDIA
abline(v=1966, lwd = 1)
abline(h=0, lwd = 1)
axis(1, at = 1956:1976, labels = seq(1956, 1976, 1), cex.axis = 1.25)
axis(2, at = seq(-3, 3, 1), cex.axis = 1.25, las = 2)
title(xlab = "Year", cex.lab = 1.5, line = 2)
title(ylab = "Difference in Military Spending", cex.lab = 1.5, line = 2)
text(1961, -2.5, "Treatment time", cex = 1.4)
arrows(1964, -2.5, 1966, -2.5, length = .1)
legend(list(x = 1955,y = 3), ncol = 1, bty = "n", c("India", "Placebo Countries"), lty=c(1,1), lwd = 2, col = c("black", "gray80"), cex = 1.25, y.intersp=0.3)

#country weights Table B12 
synthms.tables$tab.w
#Balance Table B6
synthms.tables$tab.pred






######### 2.  Average of military spending as a predictor (Appendix D.3.2 )  #####
#### ####ISRAEL  Military Spending R2 #### ####
rm(list=ls())
library(readstata13)
setwd("~/Dropbox/Women Leaders and Conflict/Empirical Analysis/Rdata and scripts")
isrms1 <- read.dta13("isr53-73fem_ms504r.dta")
library(Synth)
library(foreign)
library(tictoc)
colnames(isrms1)[4]=c("female")

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

### 2. Estimating Synthetic Counterfactual for Israel 

predictors = c( "polity2_20", "gdpgrowth_pw", "gdppc_pw", "women_legst", "cinc_c", "logtrade", "logthreat", "allies", "issues", "rival")
isrms_r2.dp <- dataprep(
    foo=isrms1, 
    predictors = predictors,
    special.predictors = list(list("logsmilex_f2", 1960:1968, c("mean")), list("logsmilex_f2", 1968, c("mean"))),
    predictors.op = "mean", 
    dependent = "logsmilex_f2", 
    unit.variable = "cown", 
    time.variable = "year",  
    treatment.identifier = 666, 
    controls.identifier = donor, 
    time.predictors.prior = c(1956:1968), 
    time.optimize.ssr = c(1956:1968), 
    unit.names.variable = "country", 
    time.plot = c(1956:tr.end[i]) 
)
isrms_r2.synth <- synth(isrms_r2.dp)


synthms2.tables <- synth.tab( dataprep.res = isrms_r2.dp, synth.res = isrms_r2.synth) 
isrms_r2.pred <- isrms_r2.dp$Y0plot %*% isrms_r2.synth$solution.w ####Y1s in Previous code
isrms_r2.alpha<- isrms_r2.dp$Y1plot - isrms_r2.pred

## 3. Plot results
plot(1956:1973, isrms_r2.dp$Y1plot, ylim = c(3, 13),  type = "n", lwd = 2, xlab = "", ylab = "")
box(); grid()
abline(v = 1969, lwd = 1)  
lines(1956:1973, isrms_r2.dp$Y1plot, lwd = 2)
lines(1956:1973, isrms_r2.pred, lwd = 2, lty = 2)
title(xlab = "Year", cex.lab = 1.5, line = 2)
title(ylab = "Military Spending", cex.lab = 2, line = 2)
legend(list(x = 1956,y = 14.5), bty = "n", c("Israel", "Synthetic Israel"), lwd = 2, col = "black", lty = c(1, 2, -1), cex = 1.8, y.intersp=0.3)
text(1962.3, 5, "First year of Meir's tenure", cex = 2)
text(1971.3, 12.5, bquote(paste("MSPE" == .(round(isrms_r2.synth$loss.v, 5)))), cex = 1.5)
arrows(1967, 5, 1969, 5, length = .1)
rug(isrms1$logsmilex_f2, side = 2, lwd=1, ticksize = .025)
#isr_ms_r2.pdf


### 4. Placebo cases ##
### Run placebo for each donor unit
for (j in 1:length(donor))
{
    isrms_r2_pl.dp<- 
        dataprep( 
            foo = isrms1, 
            predictors = predictors,
            predictors.op = "mean",
            special.predictors = list(list("logsmilex_f2", 1960:1968, c("mean")), list("logsmilex_f2", 1968, c("mean"))),
            dependent = "logsmilex_f2", 
            unit.variable = "cown", 
            time.variable = "year",  
            treatment.identifier = donor[j], 
            controls.identifier = donor[-j], 
            time.predictors.prior = c(1956:1968), 
            time.optimize.ssr = c(1956:1968), 
            unit.names.variable = "country", 
            time.plot = c(1956:tr.end[i]) 
        ) 
    
    isrms_r2_synth_pl <- synth(isrms_r2_pl.dp)
    synthms2_pl.tables <- synth.tab( dataprep.res = isrms_r2_pl.dp, synth.res = isrms_r2_synth_pl, round.digit = 5) 
    assign(paste("W",treated[i],donor[j], sep="."), (synthms2_pl.tables$tab.w))
    assign(paste("V",treated[i],donor[j], sep="."), (isrms_r2_synth_pl$solution.v))
    assign(paste("Pred",treated[i],donor[j], sep="."), (synthms2_pl.tables$tab.pred))
    
    assign(paste("Y1",treated[i],donor[j], sep="."), (isrms_r2_pl.dp$Y1))
    assign(paste("Y1s",treated[i],donor[j], sep="."), (isrms_r2_pl.dp$Y0plot%*%isrms_r2_synth_pl$solution.w))
    assign(paste("Alpha",treated[i],donor[j], sep="."), (isrms_r2_pl.dp$Y1-isrms_r2_pl.dp$Y0plot%*%isrms_r2_synth_pl$solution.w))
}


## pre-RMSPE comparison ##
isr.prermspe2 <- sqrt(mean(isrms_r2.alpha[1:13]^2, na.rm=TRUE))
isr.postrmspe2 <- sqrt(mean(isrms_r2.alpha[14:18]^2, na.rm=TRUE))
isr.prermspe.d2 <- numeric(length(donor))
isr.postrmspe.d2 <- numeric(length(donor))

for (j in 1:length(donor))
{
    temp <- get(paste("Alpha", treated[i], donor[j], sep="."))
    isr.prermspe.d2[j] <- sqrt(mean(temp[1:13]^2,na.rm=TRUE ))
    isr.postrmspe.d2[j] <- sqrt(mean(temp[14:18]^2,na.rm=TRUE ))
}

# 5.1 e(avg_pre_rmspe_p) A measure of fit. Concerning if significant.
## The proportion of placebos that have a pre- treatment RMSPE at least as large as the average of the treated units. 
rmspe.fit2<- length(which(isr.prermspe.d2>isr.prermspe2))/length(donor)


## 5.2  e(pval_joint_post_t)
## The proportion of placebos that have a ratio of post-treatment RMSPE over pre-treatment RMSPE at least as large as the average ratio for the treated units.
pval.joint.post_t2 <- length(which(isr.postrmspe.d2/isr.prermspe.d2>isr.postrmspe2/isr.prermspe2))/length(donor)

alphas2 <-isrms_r2.alpha
for (j in 1: length(donor))
{
    temp12 <- get(paste("Alpha",treated[i],donor[j], sep="."))
    alphas2 <- cbind(alphas2, temp12)
    colnames(alphas2)[j+1] <- paste(treated[i],donor[j], sep=".")
}

## 5.3  e(pvals_t)
# A vector of the proportions of placebo pseudo t-statistics (unit’s effect divided by its pre-treatment RMSPE) that are at least as largeas the main pseudo t-statistic for each post-treatment period.
pvals_t2 <- numeric(5) 
for (t in 1:5)
{
    pvals_t2[t] <-  length(which(abs(alphas2[t+13,-1]/isr.prermspe.d2)>abs(alphas2[t+13,1]/isr.prermspe2)))/length(donor)
}

## 5.3.2 Plotting p values
plot(1:5, pvals_t2, ylim=c(0,1), axes = FALSE,xlab = "", ylab = "",  pch=20, bg="grey", lwd=5, cex=2)
box(); grid()
axis(1, at = 1:5, labels = seq(1, 5, 1), cex.axis = 1.25)
axis(2, at = seq(0, 1, 0.1), cex.axis = 1, las = 2)
title(xlab = "Number of years after treatment", cex.lab = 1.5, line = 2.4)
title(ylab = "Lead specific P-Values", cex.lab = 2.1, line = 2.6)
#isr_ms_r2_pvalue.pdf

## 5.3.3 Plotting Place cases and Israel values

plot(1956:1973, isrms_r2.alpha, type = "n", ylim = c(-3, 3), lwd = 2, xlab = "", ylab = "", axes = FALSE) 
box(); grid()
#for(j in 1: 24){
#lines(1956:1973, alphas2[, j], type = "l", col = "gray80", lwd = 2)
#} 

#lines(1953:1973, alphas[,2 ], type = "l", col = "gray80", lwd = 2) #bad fit (2=US)
lines(1956:1973, alphas2[,3 ],type = "l", col = "gray80", lwd = 2)
lines(1956:1973, alphas2[,4 ],type = "l", col = "gray80", lwd = 2)
lines(1956:1973, alphas2[,5 ], type = "l", col = "gray80", lwd = 2)
#lines(1956:1973, alphas2[,6 ],type = "l", col = "gray80", lwd = 2) #bad fit (94=)
lines(1956:1973, alphas2[,7 ],type = "l", col = "gray80", lwd = 2)
lines(1956:1973, alphas2[,8 ],type = "l", col = "gray80", lwd = 2)
lines(1956:1973, alphas2[,9 ],type = "l", col = "gray80", lwd = 2)
lines(1956:1973, alphas2[,10 ],type = "l", col = "gray80", lwd = 2)
lines(1956:1973, alphas2[,11 ],type = "l", col = "gray80", lwd = 2)
lines(1956:1973, alphas2[,12 ],type = "l", col = "gray80", lwd = 2)
lines(1956:1973, alphas2[,13 ],type = "l", col = "gray80", lwd = 2)
lines(1956:1973, alphas2[,14 ],type = "l", col = "gray80", lwd = 2)
lines(1956:1973, alphas2[,15 ],type = "l", col = "gray80", lwd = 2)
lines(1956:1973, alphas2[,16 ],type = "l", col = "gray80", lwd = 2)
lines(1956:1973, alphas2[,17 ],type = "l", col = "gray80", lwd = 2)
lines(1956:1973, alphas2[,18 ],type = "l", col = "gray80", lwd = 2)
lines(1956:1973, alphas2[,19 ],type = "l", col = "gray80", lwd = 2)
lines(1956:1973, alphas2[,20 ],type = "l", col = "gray80", lwd = 2)
lines(1956:1973, alphas2[,21 ],type = "l", col = "gray80", lwd = 2)
lines(1956:1973, alphas2[,22 ],type = "l", col = "gray80", lwd = 2)
lines(1956:1973, alphas2[,23 ],type = "l", col = "gray80", lwd = 2)
lines(1956:1973, alphas2[,24 ],type = "l", col = "gray80", lwd = 2)

lines(1956:1973, alphas2[,1 ], lwd = 2) #Israel 
abline(v=1969, lwd = 1)
abline(h=0, lwd = 1)

axis(1, at = 1956:1973, labels = seq(1956, 1973, 1), cex.axis = 1.25)
axis(2, at = seq(-3, 3, 1), cex.axis = 1.25, las = 2)
title(xlab = "Year", cex.lab = 1.5, line = 2)
title(ylab = "Difference in Military Spending", cex.lab = 1.5, line = 2)
text(1964, -2.5, "Treatment time", cex = 1.8)
arrows(1967, -2.5, 1969, -2.5, length = .1)
text(1971.5, 2.5, "p-value=0.086", cex = 1.8)
legend(list(x = 1956,y = 3.5), ncol = 1, bty = "n", c("Israel", "Placebo Countries"), lty=c(1,1), lwd = 2, col = c("black", "gray80"), cex = 1.25, y.intersp=0.5)







#### #### UNITED KINGDOM   Military Spending R2 ####
rm(list=ls())
setwd("~/Dropbox/Women Leaders and Conflict/Empirical Analysis/Rdata and scripts")

library(readstata13)
ukms1 <- read.dta13("UK_ms.dta")

library(Synth)
library(foreign)
library(tictoc)

colnames(ukms1)[4]=c("female")

### 1. Prepping data
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
ukms1 <- merge(ukms1, macro, by="cown")


predictors = c("polity2_20", "gdppc_pw", "gdpgrowth_pw", "women_legst", "cinc_c", "logtrade", "logthreat", "allies", "issues", "rival")

### 2. Estimating Synthetic Counterfactual for UK 
ukms_r2.dp <- dataprep(
    foo=ukms1, 
    predictors = predictors,
    special.predictors = list(list("logsmilex_f2", 1960:1978, c("mean")), list("logsmilex_f2", 1978, c("mean"))),
    predictors.op = "mean", 
    dependent = "logsmilex_f2", 
    unit.variable = "cown", 
    time.variable = "year",  
    treatment.identifier = 200, 
    controls.identifier = donor, 
    time.predictors.prior = c(1960:1978), 
    time.optimize.ssr = c(1960:1978), 
    unit.names.variable = "country", 
    time.plot = c(1960:tr.end[i]) 
)
ukms_r2.synth <- synth(ukms_r2.dp, "ALL")

synthms.tables <- synth.tab( dataprep.res = ukms_r2.dp, synth.res = ukms_r2.synth) 
ukms_r2.pred <- ukms_r2.dp$Y0plot %*% ukms_r2.synth$solution.w ####Y1s in Previous code
ukms_r2.alpha<- ukms_r2.dp$Y1plot - ukms_r2.pred

## PLOTTING SINGLE TREATMENT GRAPH
plot(1960:1990, ukms_r2.dp$Y1plot, ylim = c(3, 13),  type = "n", lwd = 2, xlab = "", ylab = "")
box(); grid()
abline(v = 1979, lwd = 1)  
lines(1960:1990, ukms_r2.dp$Y1plot, lwd = 2)
lines(1960:1990, ukms_r2.pred, lwd = 2, lty = 2)
title(xlab = "Year", cex.lab = 1.5, line = 2)
title(ylab = "Military Spending", cex.lab = 2.5, line = 2)
legend(list(x = 1960,y = 13.5), bty = "n", c("UK", "Synthetic UK"), lwd = 2, col = "black", lty = c(1, 2, -1), cex = 1, y.intersp=0.8)
text(1968.5, 6, "First year of Thatcher's tenure", cex = 1.3)
arrows(1976, 6, 1979, 6, length = .2)
text(1985, 12.5, bquote(paste("MSPE" == .(round(ukms_r2.synth$loss.v, 5)))), cex = 1.5)
rug(ukms1$logsmilex_f2, side = 2, lwd=1, ticksize = .025)
#uk_ms_r2.pdf

### 4. Placebo cases ##
### Run placebo for each donor unit
for (j in 1:length(donor))
{
    ukms_r2_pl.dp<- 
        dataprep( 
            foo = ukms1, 
            predictors = predictors,
            predictors.op = "mean",
            special.predictors = list(list("logsmilex_f2", 1960:1978, c("mean")), list("logsmilex_f2", 1978, c("mean"))),
            dependent = "logsmilex_f2", 
            unit.variable = "cown", 
            time.variable = "year",  
            treatment.identifier = donor[j], 
            controls.identifier = donor[-j], 
            time.predictors.prior = c(1960:1978), 
            time.optimize.ssr = c(1960:1978), 
            unit.names.variable = "country", 
            time.plot = c(1960:tr.end[i]) 
        ) 
    
    ukms_r2_synth_pl <- synth(ukms_r2_pl.dp)
    synthms_pl.tables <- synth.tab( dataprep.res = ukms_r2_pl.dp, synth.res = ukms_r2_synth_pl, round.digit = 5) 
    assign(paste("W",treated[i],donor[j], sep="."), (synthms_pl.tables$tab.w))
    assign(paste("V",treated[i],donor[j], sep="."), (ukms_r2_synth_pl$solution.v))
    assign(paste("Pred",treated[i],donor[j], sep="."), (synthms_pl.tables$tab.pred))
    
    assign(paste("Y1",treated[i],donor[j], sep="."), (ukms_r2_pl.dp$Y1))
    assign(paste("Y1s",treated[i],donor[j], sep="."), (ukms_r2_pl.dp$Y0plot%*%ukms_r2_synth_pl$solution.w))
    assign(paste("Alpha",treated[i],donor[j], sep="."), (ukms_r2_pl.dp$Y1-ukms_r2_pl.dp$Y0plot%*%ukms_r2_synth_pl$solution.w))
}

### 5. Post Estimation results ###
## pre-RMSPE comparison ##
uk.prermspe <- sqrt(mean(ukms_r2.alpha[1:19]^2, na.rm=TRUE))
uk.postrmspe <- sqrt(mean(ukms_r2.alpha[20:31]^2, na.rm=TRUE))
uk.prermspe.d <- numeric(length(donor))
uk.postrmspe.d <- numeric(length(donor))

for (j in 1:length(donor))
{
    temp <- get(paste("Alpha", treated[i], donor[j], sep="."))
    uk.prermspe.d[j] <- sqrt(mean(temp[1:19]^2,na.rm=TRUE ))
    uk.postrmspe.d[j] <- sqrt(mean(temp[20:31]^2,na.rm=TRUE ))
}

# 5.1 e(avg_pre_rmspe_p) A measure of fit. Concerning if significant.
## The proportion of placebos that have a pre- treatment RMSPE at least as large as the average of the treated units. 
rmspe.fit <- length(which(uk.prermspe.d>uk.prermspe))/length(donor)


## 5.2  e(pval_joint_post_t)
## The proportion of placebos that have a ratio of post-treatment RMSPE over pre-treatment RMSPE at least as large as the average ratio for the treated units.
pval.joint.post_t <- length(which(uk.postrmspe.d/uk.prermspe.d>uk.postrmspe/uk.prermspe))/length(donor)

alphas <-ukms_r2.alpha
for (j in 1: length(donor))
{
    temp1 <- get(paste("Alpha",treated[i],donor[j], sep="."))
    alphas <- cbind(alphas, temp1)
    colnames(alphas)[j+1] <- paste(treated[i],donor[j], sep=".")
}

## 5.3.1  e(pvals_t)
# A vector of the proportions of placebo pseudo t-statistics (unit’s effect divided by its pre-treatment RMSPE) that are at least as largeas the main pseudo t-statistic for each post-treatment period.
pvals_t <- numeric(12) 
for (t in 1:12)
{
    pvals_t[t] <-  length(which(abs(alphas[t+19,-1]/uk.prermspe.d)>abs(alphas[t+19,1]/uk.prermspe)))/length(donor)
}

## 5.3.2 Plotting p values
plot(1:12, pvals_t, ylim=c(0,1), axes = FALSE,xlab = "", ylab = "",  pch=20, bg="grey", lwd=5, cex=2.5)
box(); grid()
axis(1, at = 1:10, labels = seq(1, 10, 1), cex.axis = 1.25)
axis(2, at = seq(0, 1, 0.1), cex.axis = 1, las = 2)
title(xlab = "Number of years after treatment", cex.lab = 1.5, line = 2.4)
title(ylab = "Lead specific P-Values", cex.lab = 1.5, line = 2.6)
#uk_ms_r2_pvalue.pdf

### 6. Plotting placebo cases and UK
plot(1960:1990, ukms_r2.alpha, type = "n", ylim = c(-3, 3), lwd = 2, xlab = "", ylab = "", axes = FALSE) 
box(); grid()
for(j in 3: 33){
    lines(1960:1990, alphas[, j], type = "l", col = "gray80", lwd = 2)
} 
lines(1960:1990, alphas[,1 ], lwd = 2) #UK
abline(v=1979, lwd = 1)
abline(h=0, lwd = 1)
axis(1, at = 1960:1990, labels = seq(1960, 1990, 1), cex.axis = 1.25)
axis(2, at = seq(-3, 3, 1), cex.axis = 1.25, las = 2)
title(xlab = "Year", cex.lab = 1.5, line = 2)
title(ylab = "Difference in Military Spending", cex.lab = 1.5, line = 2)
text(1971, -2.5, "Treatment time", cex = 1.5)
arrows(1976, -2.5, 1979, -2.5, length = .1)
legend(list(x = 1958,y = 3.1), ncol = 1, bty = "n", c("UK", "Placebo Countries"), lty=c(1,1), lwd = 2, col = c("black", "gray80"), cex = 1.5, y.intersp=0.3)

#country weights Table B10
synthms.tables$tab.w
#Balance Table B4
synthms.tables$tab.pred





#### #### INDIA   Military Spending R2 ####
rm(list=ls())
setwd("~/Dropbox/Women Leaders and Conflict/Empirical Analysis/Rdata and scripts")

library(readstata13)
indms1 <- read.dta13("India_ms56-84_725.dta")
indms1 <- indms1[indms1$year<1977,]

library(Synth)
library(foreign)
library(tictoc)

colnames(indms1)[4]=c("female")

### 1. Prepping data
macro <- aggregate(female ~ cown, data = indms1, sum)
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
    tr.years <- indms1$year[(indms1$female==1 & indms1$cown==treated[i])]
    tr.start[i] <- min(tr.years)
    tr.end[i] <- max(tr.years)
}
indms1 <- merge(indms1, macro, by="cown")

predictors = c("polity2_20", "gdppc_pw", "gdpgrowth_pw", "women_legst", "cinc_c", "logtrade", "logthreat_f2", "allies", "issues", "rival")

### 2. Estimating Synthetic Counterfactual for INDIA
indms_r2.dp <- dataprep(
    foo=indms1, 
    predictors = predictors,
    special.predictors = list(list("logsmilex_f2", 1956:1965, c("mean")), list("logsmilex_f2", 1965, c("mean"))),
    predictors.op = "mean", 
    dependent = "logsmilex_f2", 
    unit.variable = "cown", 
    time.variable = "year",  
    treatment.identifier = 750, 
    controls.identifier = donor, 
    time.predictors.prior = c(1956:1965), 
    time.optimize.ssr = c(1956:1965), 
    unit.names.variable = "country", 
    time.plot = c(1956:tr.end[i]) 
)
indms_r2.synth <- synth(indms_r2.dp, "ALL")

synthms.tables <- synth.tab( dataprep.res = indms_r2.dp, synth.res = indms_r2.synth) 
indms_r2.pred <- indms_r2.dp$Y0plot %*% indms_r2.synth$solution.w ####Y1s in Previous code
indms_r2.alpha<- indms_r2.dp$Y1plot - indms_r2.pred

###. 3 Plotting the results INDIA
plot(1956:1976, indms_r2.dp$Y1plot, ylim = c(3, 13),  type = "n", lwd = 2, xlab = "", ylab = "")
box(); grid()
abline(v = 1966, lwd = 1)  
lines(1956:1976, indms_r2.dp$Y1plot, lwd = 2)
lines(1956:1976, indms_r2.pred, lwd = 2, lty = 2)
title(xlab = "Year", cex.lab = 1.5, line = 2)
title(ylab = "Military Spending", cex.lab = 2.5, line = 2)
legend(list(x = 1956,y = 13), bty = "n", c("India", "Synthetic India"), lwd = 2, col = "black", lty = c(1, 2, -1), cex = 1.3, y.intersp=0.3)
text(1960.5, 5, "First year of Gandhi's tenure", cex = 1.3)
arrows(1964.8, 5, 1966, 5, length = .2)
text(1973, 12.5, bquote(paste("MSPE" == .(round(indms_r2.synth$loss.v, 5)))), cex = 1.5)
rug(indms1$logsmilex_f2, side = 2, lwd=1, ticksize = .025)
#ind_ms_r2.pdf

### 4. Placebo cases ##
### Run placebo for each donor unit

for (j in 1:length(donor))
{
    indms_r2_pl.dp<- 
        dataprep( 
            foo = indms1, 
            predictors = predictors,
            predictors.op = "mean",
            special.predictors = list(list("logsmilex_f2", 1956:1965, c("mean")), list("logsmilex_f2", 1965, c("mean"))),
            dependent = "logsmilex_f2", 
            unit.variable = "cown", 
            time.variable = "year",  
            treatment.identifier = donor[j], 
            controls.identifier = donor[-j], 
            time.predictors.prior = c(1956:1965), 
            time.optimize.ssr = c(1956:1965), 
            unit.names.variable = "country", 
            time.plot = c(1956:tr.end[i]) 
        ) 
    
    indms_r2_synth_pl <- synth(indms_r2_pl.dp)
    synthms_pl.tables <- synth.tab( dataprep.res = indms_r2_pl.dp, synth.res = indms_r2_synth_pl, round.digit = 5) 
    assign(paste("W",treated[i],donor[j], sep="."), (synthms_pl.tables$tab.w))
    assign(paste("V",treated[i],donor[j], sep="."), (indms_r2_synth_pl$solution.v))
    assign(paste("Pred",treated[i],donor[j], sep="."), (synthms_pl.tables$tab.pred))
    
    assign(paste("Y1",treated[i],donor[j], sep="."), (indms_r2_pl.dp$Y1))
    assign(paste("Y1s",treated[i],donor[j], sep="."), (indms_r2_pl.dp$Y0plot%*%indms_r2_synth_pl$solution.w))
    assign(paste("Alpha",treated[i],donor[j], sep="."), (indms_r2_pl.dp$Y1-indms_r2_pl.dp$Y0plot%*%indms_r2_synth_pl$solution.w))
    
}

### 5. Post Estimation results ###
## pre-RMSPE comparison ##
ind.prermspe <- sqrt(mean(indms_r2.alpha[1:10]^2, na.rm=TRUE))
ind.postrmspe <- sqrt(mean(indms_r2.alpha[11:21]^2, na.rm=TRUE))
ind.prermspe.d <- numeric(length(donor))
ind.postrmspe.d <- numeric(length(donor))

for (j in 1:length(donor))
{
    temp <- get(paste("Alpha", treated[i], donor[j], sep="."))
    ind.prermspe.d[j] <- sqrt(mean(temp[1:10]^2,na.rm=TRUE ))
    ind.postrmspe.d[j] <- sqrt(mean(temp[11:21]^2,na.rm=TRUE ))
}

# 5.1 e(avg_pre_rmspe_p) A measure of fit. Concerning if significant.
## The proportion of placebos that have a pre- treatment RMSPE at least as large as the average of the treated units. 
rmspe.fit <- length(which(ind.prermspe.d>ind.prermspe))/length(donor)


## 5.2  e(pval_joint_post_t)
## The proportion of placebos that have a ratio of post-treatment RMSPE over pre-treatment RMSPE at least as large as the average ratio for the treated units.
pval.joint.post_t <- length(which(ind.postrmspe.d/ind.prermspe.d>ind.postrmspe/ind.prermspe))/length(donor)

alphas <-indms_r2.alpha
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
    pvals_t[t] <-  length(which(abs(alphas[t+10,-1]/ind.prermspe.d)>abs(alphas[t+10,1]/ind.prermspe)))/length(donor)
}

plot(1:11, pvals_t, ylim=c(0,1), axes = FALSE,xlab = "", ylab = "",  pch=20, bg="grey", lwd=5, cex=2.5)
box(); grid()
axis(1, at = 1:11, labels = seq(1, 11, 1), cex.axis = 1.25)
axis(2, at = seq(0, 1, 0.1), cex.axis = 1, las = 2)
title(xlab = "Number of years after treatment", cex.lab = 1.5, line = 2.4)
title(ylab = "Lead specific P-Values", cex.lab = 2.4, line = 2.6)
#ind_ms_r2_pvalue.pdf

### 6. Plotting placebo cases and India
plot(1956:1976, indms_r2.alpha, type = "n", ylim = c(-3, 3), lwd = 2, xlab = "", ylab = "", axes = FALSE) 
box(); grid()

for(j in 3: 24)
{lines(1956:1976, alphas[, j], type = "l", col = "gray80", lwd = 2)}
lines(1956:1976, alphas[,1 ], lwd = 2) #INDIA
abline(v=1966, lwd = 1)
abline(h=0, lwd = 1)
axis(1, at = 1956:1976, labels = seq(1956, 1976, 1), cex.axis = 1.25)
axis(2, at = seq(-3, 3, 1), cex.axis = 1.25, las = 2)
title(xlab = "Year", cex.lab = 1.5, line = 2)
title(ylab = "Difference in Military Spending", cex.lab = 1.5, line = 2)
text(1961, -2.5, "Treatment time", cex = 1.4)
text(1972, 2.3, bquote(paste("p-value" == .(round(pval.joint.post_t, 3)))), cex = 1.5)
arrows(1964, -2.5, 1966, -2.5, length = .1)
legend(list(x = 1955,y = 3.2), ncol = 1, bty = "n", c("India", "Placebo Countries"), lty=c(1,1), lwd = 2, col = c("black", "gray80"), cex = 1.25, y.intersp=0.3)

#country weights Table B12 
synthms.tables$tab.w
#Balance Table B6
synthms.tables$tab.pred








######### 3.  FEMALE EMPOWERMENT MODELS (Appendix D.3.3 )
#### ISRAEL  Military Spending with Women Political empowerment index #### ####
rm(list=ls())
library(haven)
isrms1 <- read_dta("isr53-73fem_ms504r.dta")

colnames(isrms1)[4]=c("female")

# Merging with MID data
vdem <- read_dta("vdem_gen_isr.dta")
isrms1 <- merge(isrms1, vdem, all.x = TRUE )

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

### 2. Estimating Synthetic Counterfactual for Israel (Women Political Empowerment)

predictors = c( "polity2_20", "gdpgrowth_pw", "gdppc_pw", "v2x_gender", "women_legst", "cinc_c", "logtrade", "logthreat", "allies", "issues", "rival")
isrms.dp <- dataprep(
    foo=isrms1, 
    predictors = predictors,
    special.predictors = list(list("logsmilex_f2", 1968, c("mean"))),
    predictors.op = "mean", 
    dependent = "logsmilex_f2", 
    unit.variable = "cown", 
    time.variable = "year",  
    treatment.identifier = 666, 
    controls.identifier = donor, 
    time.predictors.prior = c(1956:1968), 
    time.optimize.ssr = c(1956:1968), 
    unit.names.variable = "country", 
    time.plot = c(1956:tr.end[i]) 
)

isrms.synth <- synth(isrms.dp)

synthms.tables <- synth.tab( dataprep.res = isrms.dp, synth.res = isrms.synth) 
isrms.pred <- isrms.dp$Y0plot %*% isrms.synth$solution.w ####Y1s in Previous code
isrms.alpha<- isrms.dp$Y1plot - isrms.pred



### 4. Placebo cases ##
### Run placebo for each donor unit
for (j in 1:length(donor))
{
    isrms_pl.dp<- 
        dataprep( 
            foo = isrms1, 
            predictors = predictors,
            predictors.op = "mean",
            special.predictors = list(list("logsmilex_f2", 1968, c("mean"))),
            dependent = "logsmilex_f2", 
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
    synthms2_pl.tables <- synth.tab( dataprep.res = isrms_pl.dp, synth.res = isrms_synth_pl, round.digit = 5) 
    assign(paste("W",treated[i],donor[j], sep="."), (synthms2_pl.tables$tab.w))
    assign(paste("V",treated[i],donor[j], sep="."), (isrms_synth_pl$solution.v))
    assign(paste("Pred",treated[i],donor[j], sep="."), (synthms2_pl.tables$tab.pred))
    
    assign(paste("Y1",treated[i],donor[j], sep="."), (isrms_pl.dp$Y1))
    assign(paste("Y1s",treated[i],donor[j], sep="."), (isrms_pl.dp$Y0plot%*%isrms_synth_pl$solution.w))
    assign(paste("Alpha",treated[i],donor[j], sep="."), (isrms_pl.dp$Y1-isrms_pl.dp$Y0plot%*%isrms_synth_pl$solution.w))
}

# 5.1 e(avg_pre_rmspe_p) A measure of fit. Concerning if significant.
## The proportion of placebos that have a pre- treatment RMSPE at least as large as the average of the treated units. 

## pre-RMSPE comparison ##
isr.prermspe2 <- sqrt(mean(isrms.alpha[1:13]^2, na.rm=TRUE))
isr.postrmspe2 <- sqrt(mean(isrms.alpha[14:18]^2, na.rm=TRUE))
isr.prermspe.d2 <- numeric(length(donor))
isr.postrmspe.d2 <- numeric(length(donor))

for (j in 1:length(donor))
{
    temp <- get(paste("Alpha", treated[i], donor[j], sep="."))
    isr.prermspe.d2[j] <- sqrt(mean(temp[1:13]^2,na.rm=TRUE ))
    isr.postrmspe.d2[j] <- sqrt(mean(temp[14:18]^2,na.rm=TRUE ))
}

rmspe.fit<- length(which(isr.prermspe.d2>isr.prermspe2))/length(donor)
## 5.2  e(pval_joint_post_t)
## The proportion of placebos that have a ratio of post-treatment RMSPE over pre-treatment RMSPE at least as large as the average ratio for the treated units.
pval.joint.post_t2 <- length(which(isr.postrmspe.d2/isr.prermspe.d2>isr.postrmspe2/isr.prermspe2))/length(donor)

alphas2 <-isrms.alpha
for (j in 1: length(donor))
{
    temp12 <- get(paste("Alpha",treated[i],donor[j], sep="."))
    alphas2 <- cbind(alphas2, temp12)
    colnames(alphas2)[j+1] <- paste(treated[i],donor[j], sep=".")
}

## 5.3  e(pvals_t)
# A vector of the proportions of placebo pseudo t-statistics (unit’s effect divided by its pre-treatment RMSPE) that are at least as largeas the main pseudo t-statistic for each post-treatment period.
pvals_t2 <- numeric(5) 
for (t in 1:5)
{
    pvals_t2[t] <-  length(which(abs(alphas2[t+13,-1]/isr.prermspe.d2)>abs(alphas2[t+13,1]/isr.prermspe2)))/length(donor)
}

### 7. Getting country weights and balance tables
#country weights 
donorw <- synthms.tables$tab.w
#Balance Table 
synthms.tables$tab.pred

# PLOT 1 (isr_ms_wp1)
plot(1956:1973, isrms.dp$Y1plot, ylim = c(3, 13),  type = "n", lwd = 2, xlab = "", ylab = "")
box(); grid()
abline(v = 1969, lwd = 1)  
lines(1956:1973, isrms.dp$Y1plot, lwd = 2)
lines(1956:1973, isrms.pred, lwd = 2, lty = 2)
title(xlab = "Year", cex.lab = 1.5, line = 2)
title(ylab = "Military Spending", cex.lab = 1, line = 2)
legend(list(x = 1955,y = 13.5), bty = "n", c("Israel", "Synthetic Israel"), lwd = 2, col = "black", lty = c(1, 2, -1), cex = 1.5, y.intersp=0.5)
text(1962, 5, "First year of Meir's tenure", cex = 2)
text(1971.2, 12.5, bquote(paste("MSPE" == .(round(isrms.synth$loss.v, 3)))), cex = 1.3)
text(1971.2, 3.5, bquote(paste("p-value" == .(round(pval.joint.post_t2, 3)))), cex = 1.3)
arrows(1967, 5, 1969, 5, length = .1)
rug(isrms1$logsmilex_f2, side = 2, lwd=1, ticksize = .025)

#### Israel  Military Spending with Women Civil Liberties Index #### ####

### 2. Estimating Synthetic Counterfactual for Israel

predictors = c( "polity2_20", "gdpgrowth_pw", "gdppc_pw", "v2x_gencs", "women_legst", "cinc_c", "logtrade", "logthreat", "allies", "issues", "rival")
isrms.dp <- dataprep(
    foo=isrms1, 
    predictors = predictors,
    special.predictors = list(list("logsmilex_f2", 1968, c("mean"))),
    predictors.op = "mean", 
    dependent = "logsmilex_f2", 
    unit.variable = "cown", 
    time.variable = "year",  
    treatment.identifier = 666, 
    controls.identifier = donor, 
    time.predictors.prior = c(1956:1968), 
    time.optimize.ssr = c(1956:1968), 
    unit.names.variable = "country", 
    time.plot = c(1956:tr.end[i]) 
)

isrms.synth <- synth(isrms.dp)

synthms.tables <- synth.tab( dataprep.res = isrms.dp, synth.res = isrms.synth) 
isrms.pred <- isrms.dp$Y0plot %*% isrms.synth$solution.w ####Y1s in Previous code
isrms.alpha<- isrms.dp$Y1plot - isrms.pred

#path.plot(dataprep.res = isrms.dp, synth.res = isrms.synth)

### 4. Placebo cases ##
### Run placebo for each donor unit
for (j in 1:length(donor))
{
    isrms_pl.dp<- 
        dataprep( 
            foo = isrms1, 
            predictors = predictors,
            predictors.op = "mean",
            special.predictors = list(list("logsmilex_f2", 1968, c("mean"))),
            dependent = "logsmilex_f2", 
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
    synthms2_pl.tables <- synth.tab( dataprep.res = isrms_pl.dp, synth.res = isrms_synth_pl, round.digit = 5) 
    assign(paste("W",treated[i],donor[j], sep="."), (synthms2_pl.tables$tab.w))
    assign(paste("V",treated[i],donor[j], sep="."), (isrms_synth_pl$solution.v))
    assign(paste("Pred",treated[i],donor[j], sep="."), (synthms2_pl.tables$tab.pred))
    
    assign(paste("Y1",treated[i],donor[j], sep="."), (isrms_pl.dp$Y1))
    assign(paste("Y1s",treated[i],donor[j], sep="."), (isrms_pl.dp$Y0plot%*%isrms_synth_pl$solution.w))
    assign(paste("Alpha",treated[i],donor[j], sep="."), (isrms_pl.dp$Y1-isrms_pl.dp$Y0plot%*%isrms_synth_pl$solution.w))
}

# 5.1 e(avg_pre_rmspe_p) A measure of fit. Concerning if significant.
## The proportion of placebos that have a pre- treatment RMSPE at least as large as the average of the treated units. 

## pre-RMSPE comparison ##
isr.prermspe2 <- sqrt(mean(isrms.alpha[1:13]^2, na.rm=TRUE))
isr.postrmspe2 <- sqrt(mean(isrms.alpha[14:18]^2, na.rm=TRUE))
isr.prermspe.d2 <- numeric(length(donor))
isr.postrmspe.d2 <- numeric(length(donor))

for (j in 1:length(donor))
{
    temp <- get(paste("Alpha", treated[i], donor[j], sep="."))
    isr.prermspe.d2[j] <- sqrt(mean(temp[1:13]^2,na.rm=TRUE ))
    isr.postrmspe.d2[j] <- sqrt(mean(temp[14:18]^2,na.rm=TRUE ))
}

rmspe.fit<- length(which(isr.prermspe.d2>isr.prermspe2))/length(donor)


## 5.2  e(pval_joint_post_t)
## The proportion of placebos that have a ratio of post-treatment RMSPE over pre-treatment RMSPE at least as large as the average ratio for the treated units.
pval.joint.post_t2 <- length(which(isr.postrmspe.d2/isr.prermspe.d2>isr.postrmspe2/isr.prermspe2))/length(donor)

alphas2 <-isrms.alpha
for (j in 1: length(donor))
{
    temp12 <- get(paste("Alpha",treated[i],donor[j], sep="."))
    alphas2 <- cbind(alphas2, temp12)
    colnames(alphas2)[j+1] <- paste(treated[i],donor[j], sep=".")
}

## 5.3  e(pvals_t)
# A vector of the proportions of placebo pseudo t-statistics (unit’s effect divided by its pre-treatment RMSPE) that are at least as largeas the main pseudo t-statistic for each post-treatment period.
pvals_t2 <- numeric(5) 
for (t in 1:5)
{
    pvals_t2[t] <-  length(which(abs(alphas2[t+13,-1]/isr.prermspe.d2)>abs(alphas2[t+13,1]/isr.prermspe2)))/length(donor)
}

### 7. Getting country weights and balance tables
#country weights 
donorw <- synthms.tables$tab.w
#Balance Table 
synthms.tables$tab.pred

# PLOT 2 (isr_ms_wp2)
plot(1956:1973, isrms.dp$Y1plot, ylim = c(3, 13),  type = "n", lwd = 2, xlab = "", ylab = "")
box(); grid()
abline(v = 1969, lwd = 1)  
lines(1956:1973, isrms.dp$Y1plot, lwd = 2)
lines(1956:1973, isrms.pred, lwd = 2, lty = 2)
title(xlab = "Year", cex.lab = 1.5, line = 2)
title(ylab = "Military Spending", cex.lab = 1, line = 2)
legend(list(x = 1955,y = 13.5), bty = "n", c("Israel", "Synthetic Israel"), lwd = 2, col = "black", lty = c(1, 2, -1), cex = 1.5, y.intersp=0.5)
text(1962, 5, "First year of Meir's tenure", cex = 2)
text(1971.2, 12.5, bquote(paste("MSPE" == .(round(isrms.synth$loss.v, 3)))), cex = 1.3)
text(1971.2, 3.5, bquote(paste("p-value" == .(round(pval.joint.post_t2, 3)))), cex = 1.3)
arrows(1967, 5, 1969, 5, length = .1)
rug(isrms1$logsmilex_f2, side = 2, lwd=1, ticksize = .025)

#### UK  Military Spending with Women Political empowerment index ####

rm(list=ls())
library(haven)
ukms1 <- read_dta("UK_ms.dta")
library(Synth)

colnames(ukms1)[4]=c("female")

# Merging with MID data
vdem <- read_dta("vdem_gender_uk.dta")
ukms1 <- merge(ukms1, vdem, all.x = TRUE )

### 1. Prepping data
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
ukms1 <- merge(ukms1, macro, by="cown")

predictors = c("polity2_20", "gdppc_pw", "gdpgrowth_pw", "v2x_gender", "women_legst", "cinc_c", "logtrade", "logthreat", "allies", "issues", "rival")

### 2. Estimating Synthetic Counterfactual for UK 
ukms.dp <- dataprep(
    foo=ukms1, 
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
    time.plot = c(1960:tr.end[i]) 
)
ukms.synth <- synth(ukms.dp)

synthms.tables <- synth.tab( dataprep.res = ukms.dp, synth.res = ukms.synth) 
ukms.pred <- ukms.dp$Y0plot %*% ukms.synth$solution.w ####Y1s in Previous code
ukms.alpha<- ukms.dp$Y1plot - ukms.pred

### 4. Placebo cases ##
### Run placebo for each donor unit
for (j in 1:length(donor))
{
    ukms_pl.dp<- 
        dataprep( 
            foo = ukms1, 
            predictors = predictors,
            predictors.op = "mean",
            special.predictors = list(list("logsmilex_f2",1978, c("mean"))),
            dependent = "logsmilex_f2", 
            unit.variable = "cown", 
            time.variable = "year",  
            treatment.identifier = donor[j], 
            controls.identifier = donor[-j], 
            time.predictors.prior = c(1960:1978), 
            time.optimize.ssr = c(1960:1978), 
            unit.names.variable = "country", 
            time.plot = c(1960:tr.end[i]) 
        ) 
    
    ukms_synth_pl <- synth(ukms_pl.dp)
    synthms_pl.tables <- synth.tab( dataprep.res = ukms_pl.dp, synth.res = ukms_synth_pl, round.digit = 5) 
    assign(paste("W",treated[i],donor[j], sep="."), (synthms_pl.tables$tab.w))
    assign(paste("V",treated[i],donor[j], sep="."), (ukms_synth_pl$solution.v))
    assign(paste("Pred",treated[i],donor[j], sep="."), (synthms_pl.tables$tab.pred))
    
    assign(paste("Y1",treated[i],donor[j], sep="."), (ukms_pl.dp$Y1))
    assign(paste("Y1s",treated[i],donor[j], sep="."), (ukms_pl.dp$Y0plot%*%ukms_synth_pl$solution.w))
    assign(paste("Alpha",treated[i],donor[j], sep="."), (ukms_pl.dp$Y1-ukms_pl.dp$Y0plot%*%ukms_synth_pl$solution.w))
}

### 5. Post Estimation results ###
## pre-RMSPE comparison ##
uk.prermspe <- sqrt(mean(ukms.alpha[1:19]^2, na.rm=TRUE))
uk.postrmspe <- sqrt(mean(ukms.alpha[20:31]^2, na.rm=TRUE))
uk.prermspe.d <- numeric(length(donor))
uk.postrmspe.d <- numeric(length(donor))

for (j in 1:length(donor))
{
    temp <- get(paste("Alpha", treated[i], donor[j], sep="."))
    uk.prermspe.d[j] <- sqrt(mean(temp[1:19]^2,na.rm=TRUE ))
    uk.postrmspe.d[j] <- sqrt(mean(temp[20:31]^2,na.rm=TRUE ))
}

# 5.1 e(avg_pre_rmspe_p) A measure of fit. Concerning if significant.
## The proportion of placebos that have a pre- treatment RMSPE at least as large as the average of the treated units. 
rmspe.fit <- length(which(uk.prermspe.d>uk.prermspe))/length(donor)


## 5.2  e(pval_joint_post_t)
## The proportion of placebos that have a ratio of post-treatment RMSPE over pre-treatment RMSPE at least as large as the average ratio for the treated units.
pval.joint.post_t <- length(which(uk.postrmspe.d/uk.prermspe.d>uk.postrmspe/uk.prermspe))/length(donor)

alphas <-ukms.alpha
for (j in 1: length(donor))
{
    temp1 <- get(paste("Alpha",treated[i],donor[j], sep="."))
    alphas <- cbind(alphas, temp1)
    colnames(alphas)[j+1] <- paste(treated[i],donor[j], sep=".")
}

## 5.3.1  e(pvals_t)
# A vector of the proportions of placebo pseudo t-statistics (unit’s effect divided by its pre-treatment RMSPE) that are at least as largeas the main pseudo t-statistic for each post-treatment period.
pvals_t <- numeric(12) 
for (t in 1:12)
{
    pvals_t[t] <-  length(which(abs(alphas[t+19,-1]/uk.prermspe.d)>abs(alphas[t+19,1]/uk.prermspe)))/length(donor)
}

#country weights =
donorw <- synthms.tables$tab.w
#Balance Table
synthms.tables$tab.pred

###PLOT 1. 
plot(1960:1990, ukms.dp$Y1plot, ylim = c(3, 13),  type = "n", lwd = 2, xlab = "", ylab = "")
box(); grid()
abline(v = 1979, lwd = 1)  
lines(1960:1990, ukms.dp$Y1plot, lwd = 2)
lines(1960:1990, ukms.pred, lwd = 2, lty = 2)
title(xlab = "Year", cex.lab = 1.5, line = 2)
title(ylab = "Military Spending", cex.lab = 2.5, line = 2)
legend(list(x = 1960,y = 13), bty = "n", c("UK", "Synthetic UK"), lwd = 2, col = "black", lty = c(1, 2, -1), cex = 1.3, y.intersp=0.5)
text(1969, 6, "First year of Thatcher's tenure", cex = 1.5)
arrows(1976, 6, 1979, 6, length = .2)
text(1985, 12.5, bquote(paste("MSPE" == .(round(ukms.synth$loss.v, 3)))), cex = 1.5)
text(1985, 4, bquote(paste("p-value" == .(round(pval.joint.post_t, 3)))), cex = 1.5)
rug(ukms1$logsmilex_f2, side = 2, lwd=1, ticksize = .025)

#### UK  Military Spending with Women Civil Liberties Index ####

predictors = c("polity2_20", "gdppc_pw", "gdpgrowth_pw", "v2x_gencl", "women_legst", "cinc_c", "logtrade", "logthreat", "allies", "issues", "rival")

### 2. Estimating Synthetic Counterfactual for UK 
ukms.dp <- dataprep(
    foo=ukms1, 
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
    time.plot = c(1960:tr.end[i]) 
)
ukms.synth <- synth(ukms.dp)

synthms.tables <- synth.tab( dataprep.res = ukms.dp, synth.res = ukms.synth) 
ukms.pred <- ukms.dp$Y0plot %*% ukms.synth$solution.w ####Y1s in Previous code
ukms.alpha<- ukms.dp$Y1plot - ukms.pred



### 4. Placebo cases ##
### Run placebo for each donor unit
for (j in 1:length(donor))
{
    ukms_pl.dp<- 
        dataprep( 
            foo = ukms1, 
            predictors = predictors,
            predictors.op = "mean",
            special.predictors = list(list("logsmilex_f2",1978, c("mean"))),
            dependent = "logsmilex_f2", 
            unit.variable = "cown", 
            time.variable = "year",  
            treatment.identifier = donor[j], 
            controls.identifier = donor[-j], 
            time.predictors.prior = c(1960:1978), 
            time.optimize.ssr = c(1960:1978), 
            unit.names.variable = "country", 
            time.plot = c(1960:tr.end[i]) 
        ) 
    
    ukms_synth_pl <- synth(ukms_pl.dp)
    synthms_pl.tables <- synth.tab( dataprep.res = ukms_pl.dp, synth.res = ukms_synth_pl, round.digit = 5) 
    assign(paste("W",treated[i],donor[j], sep="."), (synthms_pl.tables$tab.w))
    assign(paste("V",treated[i],donor[j], sep="."), (ukms_synth_pl$solution.v))
    assign(paste("Pred",treated[i],donor[j], sep="."), (synthms_pl.tables$tab.pred))
    
    assign(paste("Y1",treated[i],donor[j], sep="."), (ukms_pl.dp$Y1))
    assign(paste("Y1s",treated[i],donor[j], sep="."), (ukms_pl.dp$Y0plot%*%ukms_synth_pl$solution.w))
    assign(paste("Alpha",treated[i],donor[j], sep="."), (ukms_pl.dp$Y1-ukms_pl.dp$Y0plot%*%ukms_synth_pl$solution.w))
}

### 5. Post Estimation results ###
## pre-RMSPE comparison ##
uk.prermspe <- sqrt(mean(ukms.alpha[1:19]^2, na.rm=TRUE))
uk.postrmspe <- sqrt(mean(ukms.alpha[20:31]^2, na.rm=TRUE))
uk.prermspe.d <- numeric(length(donor))
uk.postrmspe.d <- numeric(length(donor))

for (j in 1:length(donor))
{
    temp <- get(paste("Alpha", treated[i], donor[j], sep="."))
    uk.prermspe.d[j] <- sqrt(mean(temp[1:19]^2,na.rm=TRUE ))
    uk.postrmspe.d[j] <- sqrt(mean(temp[20:31]^2,na.rm=TRUE ))
}


# 5.1 e(avg_pre_rmspe_p) A measure of fit. Concerning if significant.
## The proportion of placebos that have a pre- treatment RMSPE at least as large as the average of the treated units. 
rmspe.fit <- length(which(uk.prermspe.d>uk.prermspe))/length(donor)


## 5.2  e(pval_joint_post_t)
## The proportion of placebos that have a ratio of post-treatment RMSPE over pre-treatment RMSPE at least as large as the average ratio for the treated units.
pval.joint.post_t <- length(which(uk.postrmspe.d/uk.prermspe.d>uk.postrmspe/uk.prermspe))/length(donor)

alphas <-ukms.alpha
for (j in 1: length(donor))
{
    temp1 <- get(paste("Alpha",treated[i],donor[j], sep="."))
    alphas <- cbind(alphas, temp1)
    colnames(alphas)[j+1] <- paste(treated[i],donor[j], sep=".")
}

## 5.3.1  e(pvals_t)
# A vector of the proportions of placebo pseudo t-statistics (unit’s effect divided by its pre-treatment RMSPE) that are at least as largeas the main pseudo t-statistic for each post-treatment period.
pvals_t <- numeric(12) 
for (t in 1:12)
{
    pvals_t[t] <-  length(which(abs(alphas[t+19,-1]/uk.prermspe.d)>abs(alphas[t+19,1]/uk.prermspe)))/length(donor)
}

#country weights =
donorw <- synthms.tables$tab.w
#Balance Table
synthms.tables$tab.pred

# PLOT 2 (uk_ms_wp2)
plot(1960:1990, ukms.dp$Y1plot, ylim = c(3, 13),  type = "n", lwd = 2, xlab = "", ylab = "")
box(); grid()
abline(v = 1979, lwd = 1)  
lines(1960:1990, ukms.dp$Y1plot, lwd = 2)
lines(1960:1990, ukms.pred, lwd = 2, lty = 2)
title(xlab = "Year", cex.lab = 1.5, line = 2)
title(ylab = "Military Spending", cex.lab = 2, line = 2)
legend(list(x = 1960,y = 13), bty = "n", c("UK", "Synthetic UK"), lwd = 2, col = "black", lty = c(1, 2, -1), cex = 1.3, y.intersp=0.5)
text(1969, 6, "First year of Thatcher's tenure", cex = 1.5)
arrows(1976, 6, 1979, 6, length = .2)
text(1985, 12.5, bquote(paste("MSPE" == .(round(ukms.synth$loss.v, 3)))), cex = 1.5)
text(1985, 4, bquote(paste("p-value" == .(round(pval.joint.post_t, 3)))), cex = 1.5)
rug(ukms1$logsmilex_f2, side = 2, lwd=1, ticksize = .025)

#### INDIA  Military Spending with Women Political empowerment index #### ####
rm(list=ls())
library(haven)
library(Synth)
library(tictoc)

indms1 <- read_dta("India_ms56-84_725.dta")
indms1 <- indms1[indms1$year<1977,]
colnames(indms1)[4]=c("female")

vdem <- read_dta("vdem_gen_ind.dta")
indms1 <- merge(indms1, vdem, all.x = TRUE )

### 1. Prepping data
macro <- aggregate(female ~ cown, data = indms1, sum)
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
    tr.years <- indms1$year[(indms1$female==1 & indms1$cown==treated[i])]
    tr.start[i] <- min(tr.years)
    tr.end[i] <- max(tr.years)
}
indms1 <- merge(indms1, macro, by="cown")


predictors = c("polity2_20", "gdppc_pw", "gdpgrowth_pw", "v2x_gender", "women_legst", "cinc_c", "logtrade", "logthreat_f2", "allies", "issues", "rival")

### 2. Estimating Synthetic Counterfactual for INDIA
indms.dp <- dataprep(
    foo=indms1, 
    predictors = predictors,
    special.predictors = list(list("logsmilex_f2",1965, c("mean"))),
    predictors.op = "mean", 
    dependent = "logsmilex_f2", 
    unit.variable = "cown", 
    time.variable = "year",  
    treatment.identifier = 750, 
    controls.identifier = donor, 
    time.predictors.prior = c(1956:1965), 
    time.optimize.ssr = c(1956:1965), 
    unit.names.variable = "country", 
    time.plot = c(1956:tr.end[i]) 
)
indms.synth <- synth(indms.dp, "ALL")

synthms.tables <- synth.tab( dataprep.res = indms.dp, synth.res = indms.synth) 
indms.pred <- indms.dp$Y0plot %*% indms.synth$solution.w ####Y1s in Previous code
indms.alpha<- indms.dp$Y1plot - indms.pred


### 4. Placebo cases ##
### Run placebo for each donor unit

for (j in 1:length(donor))
{
    indms_pl.dp<- 
        dataprep( 
            foo = indms1, 
            predictors = predictors,
            predictors.op = "mean",
            special.predictors = list(list("logsmilex_f2",1965, c("mean"))),
            dependent = "logsmilex_f2", 
            unit.variable = "cown", 
            time.variable = "year",  
            treatment.identifier = donor[j], 
            controls.identifier = donor[-j], 
            time.predictors.prior = c(1956:1965), 
            time.optimize.ssr = c(1956:1965), 
            unit.names.variable = "country", 
            time.plot = c(1956:tr.end[i]) 
        ) 
    
    indms_synth_pl <- synth(indms_pl.dp)
    synthms_pl.tables <- synth.tab( dataprep.res = indms_pl.dp, synth.res = indms_synth_pl, round.digit = 5) 
    assign(paste("W",treated[i],donor[j], sep="."), (synthms_pl.tables$tab.w))
    assign(paste("V",treated[i],donor[j], sep="."), (indms_synth_pl$solution.v))
    assign(paste("Pred",treated[i],donor[j], sep="."), (synthms_pl.tables$tab.pred))
    
    assign(paste("Y1",treated[i],donor[j], sep="."), (indms_pl.dp$Y1))
    assign(paste("Y1s",treated[i],donor[j], sep="."), (indms_pl.dp$Y0plot%*%indms_synth_pl$solution.w))
    assign(paste("Alpha",treated[i],donor[j], sep="."), (indms_pl.dp$Y1-indms_pl.dp$Y0plot%*%indms_synth_pl$solution.w))
    #bad.fit <- numeric(33)
    #bad.fit[j] <- as.numeric(ifelse(indms_synth_pl$loss.v > indms.synth$loss.v*5, 1, 0))
    
}

### 5. Post Estimation results ###
## pre-RMSPE comparison ##
ind.prermspe <- sqrt(mean(indms.alpha[1:10]^2, na.rm=TRUE))
ind.postrmspe <- sqrt(mean(indms.alpha[11:21]^2, na.rm=TRUE))
ind.prermspe.d <- numeric(length(donor))
ind.postrmspe.d <- numeric(length(donor))

for (j in 1:length(donor))
{
    temp <- get(paste("Alpha", treated[i], donor[j], sep="."))
    ind.prermspe.d[j] <- sqrt(mean(temp[1:10]^2,na.rm=TRUE ))
    ind.postrmspe.d[j] <- sqrt(mean(temp[11:21]^2,na.rm=TRUE ))
}

# 5.1 e(avg_pre_rmspe_p) A measure of fit. Concerning if significant.
## The proportion of placebos that have a pre- treatment RMSPE at least as large as the average of the treated units. 
rmspe.fit <- length(which(ind.prermspe.d>ind.prermspe))/length(donor)

## 5.2  e(pval_joint_post_t)
## The proportion of placebos that have a ratio of post-treatment RMSPE over pre-treatment RMSPE at least as large as the average ratio for the treated units.
pval.joint.post_t <- length(which(ind.postrmspe.d/ind.prermspe.d>ind.postrmspe/ind.prermspe))/length(donor)

alphas <-indms.alpha
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
    pvals_t[t] <-  length(which(abs(alphas[t+10,-1]/ind.prermspe.d)>abs(alphas[t+10,1]/ind.prermspe)))/length(donor)
}

#country weights Table
synthms.tables$tab.w

# PLOT 1 (ind_ms_wp1)
plot(1956:1976, indms.dp$Y1plot, ylim = c(3, 13),  type = "n", lwd = 2, xlab = "", ylab = "")
box(); grid()
abline(v = 1966, lwd = 1)  
lines(1956:1976, indms.dp$Y1plot, lwd = 2)
lines(1956:1976, indms.pred, lwd = 2, lty = 2)
title(xlab = "Year", cex.lab = 1.5, line = 2)
title(ylab = "Military Spending", cex.lab = 2.5, line = 2)
legend(list(x = 1956,y = 13), bty = "n", c("India", "Synthetic India"), lwd = 2, col = "black", lty = c(1, 2, -1), cex = 1.3, y.intersp=0.7)
text(1960.5, 5, "First year of Gandhi's tenure", cex = 1.3)
arrows(1964.8, 5, 1966, 5, length = .2)
text(1973, 12.5, bquote(paste("MSPE" == .(round(indms.synth$loss.v, 3)))), cex = 1.5)
text(1973, 3.5, bquote(paste("p-value" == .(round(0.999, 3)))), cex = 1.5)
rug(indms1$logsmilex_f2, side = 2, lwd=1, ticksize = .025)



#### India  Military Spending with Women Civil Liberties Index #### ####
predictors = c( "polity2_20", "gdpgrowth_pw", "gdppc_pw", "v2x_gencs", "women_legst", "cinc_c", "logtrade", "logthreat", "allies", "issues", "rival")

### 2. Estimating Synthetic Counterfactual for INDIA
indms.dp <- dataprep(
    foo=indms1, 
    predictors = predictors,
    special.predictors = list(list("logsmilex_f2",1965, c("mean"))),
    predictors.op = "mean", 
    dependent = "logsmilex_f2", 
    unit.variable = "cown", 
    time.variable = "year",  
    treatment.identifier = 750, 
    controls.identifier = donor, 
    time.predictors.prior = c(1956:1965), 
    time.optimize.ssr = c(1956:1965), 
    unit.names.variable = "country", 
    time.plot = c(1956:tr.end[i]) 
)
indms.synth <- synth(indms.dp, "ALL")

synthms.tables <- synth.tab( dataprep.res = indms.dp, synth.res = indms.synth) 
indms.pred <- indms.dp$Y0plot %*% indms.synth$solution.w ####Y1s in Previous code
indms.alpha<- indms.dp$Y1plot - indms.pred

### 4. Placebo cases ##
### Run placebo for each donor unit
for (j in 1:length(donor))
{
    indms_pl.dp<- 
        dataprep( 
            foo = indms1, 
            predictors = predictors,
            predictors.op = "mean",
            special.predictors = list(list("logsmilex_f2",1965, c("mean"))),
            dependent = "logsmilex_f2", 
            unit.variable = "cown", 
            time.variable = "year",  
            treatment.identifier = donor[j], 
            controls.identifier = donor[-j], 
            time.predictors.prior = c(1956:1965), 
            time.optimize.ssr = c(1956:1965), 
            unit.names.variable = "country", 
            time.plot = c(1956:tr.end[i]) 
        ) 
    
    indms_synth_pl <- synth(indms_pl.dp)
    synthms_pl.tables <- synth.tab( dataprep.res = indms_pl.dp, synth.res = indms_synth_pl, round.digit = 5) 
    assign(paste("W",treated[i],donor[j], sep="."), (synthms_pl.tables$tab.w))
    assign(paste("V",treated[i],donor[j], sep="."), (indms_synth_pl$solution.v))
    assign(paste("Pred",treated[i],donor[j], sep="."), (synthms_pl.tables$tab.pred))
    
    assign(paste("Y1",treated[i],donor[j], sep="."), (indms_pl.dp$Y1))
    assign(paste("Y1s",treated[i],donor[j], sep="."), (indms_pl.dp$Y0plot%*%indms_synth_pl$solution.w))
    assign(paste("Alpha",treated[i],donor[j], sep="."), (indms_pl.dp$Y1-indms_pl.dp$Y0plot%*%indms_synth_pl$solution.w))
}

### 5. Post Estimation results ###
## pre-RMSPE comparison ##
ind.prermspe <- sqrt(mean(indms.alpha[1:10]^2, na.rm=TRUE))
ind.postrmspe <- sqrt(mean(indms.alpha[11:21]^2, na.rm=TRUE))
ind.prermspe.d <- numeric(length(donor))
ind.postrmspe.d <- numeric(length(donor))

for (j in 1:length(donor))
{
    temp <- get(paste("Alpha", treated[i], donor[j], sep="."))
    ind.prermspe.d[j] <- sqrt(mean(temp[1:10]^2,na.rm=TRUE ))
    ind.postrmspe.d[j] <- sqrt(mean(temp[11:21]^2,na.rm=TRUE ))
}

# 5.1 e(avg_pre_rmspe_p) A measure of fit. Concerning if significant.
## The proportion of placebos that have a pre- treatment RMSPE at least as large as the average of the treated units. 
rmspe.fit <- length(which(ind.prermspe.d>ind.prermspe))/length(donor)

## 5.2  e(pval_joint_post_t)
## The proportion of placebos that have a ratio of post-treatment RMSPE over pre-treatment RMSPE at least as large as the average ratio for the treated units.
pval.joint.post_t <- length(which(ind.postrmspe.d/ind.prermspe.d>ind.postrmspe/ind.prermspe))/length(donor)


alphas <-indms.alpha
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
    pvals_t[t] <-  length(which(abs(alphas[t+10,-1]/ind.prermspe.d)>abs(alphas[t+10,1]/ind.prermspe)))/length(donor)
}

#country weights Table
synthms.tables$tab.w

# PLOT 2 (ind_ms_wp2)
plot(1956:1976, indms.dp$Y1plot, ylim = c(3, 13),  type = "n", lwd = 2, xlab = "", ylab = "")
box(); grid()
abline(v = 1966, lwd = 1)  
lines(1956:1976, indms.dp$Y1plot, lwd = 2)
lines(1956:1976, indms.pred, lwd = 2, lty = 2)
title(xlab = "Year", cex.lab = 1.5, line = 2)
title(ylab = "Military Spending", cex.lab = 2.5, line = 2)
legend(list(x = 1956,y = 13), bty = "n", c("India", "Synthetic India"), lwd = 2, col = "black", lty = c(1, 2, -1), cex = 1.3, y.intersp=0.7)
text(1960.5, 5, "First year of Gandhi's tenure", cex = 1.3)
arrows(1964.8, 5, 1966, 5, length = .2)
text(1973, 12.5, bquote(paste("MSPE" == .(round(indms.synth$loss.v, 3)))), cex = 1.5)
text(1973, 3.5, bquote(paste("p-value" == .(round(0.999, 3)))), cex = 1.5)
rug(indms1$logsmilex_f2, side = 2, lwd=1, ticksize = .025)


######### 4.  MID ANALYSIS (Appendix D.3.4 )
#### ####ISRAEL  Military Spending with COUNT MID #### ####
rm(list=ls())
library(haven)
isrms1 <- read_dta("isr53-73fem_ms504r.dta")
library(Synth)
library(tictoc)

colnames(isrms1)[4]=c("female")

# Merging with MID data
mid <- read_dta("mid_isr.dta")
isrms1 <- merge(isrms1, mid, all.x = TRUE )

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

### 2. Estimating Synthetic Counterfactual for Israel with Count MID per year included 

predictors = c( "polity2_20", "gdpgrowth_pw", "gdppc_pw", "count_mid", "women_legst", "cinc_c", "logtrade", "logthreat", "allies", "issues", "rival")
isrms_r1.dp <- dataprep(
    foo=isrms1, 
    predictors = predictors,
    special.predictors = list(list("logsmilex_f2", 1968, c("mean"))),
    predictors.op = "mean", 
    dependent = "logsmilex_f2", 
    unit.variable = "cown", 
    time.variable = "year",  
    treatment.identifier = 666, 
    controls.identifier = donor, 
    time.predictors.prior = c(1956:1968), 
    time.optimize.ssr = c(1956:1968), 
    unit.names.variable = "country", 
    time.plot = c(1956:tr.end[i]) 
)

isrms_r1.synth <- synth(isrms_r1.dp)
synthms.tables <- synth.tab( dataprep.res = isrms_r1.dp, synth.res = isrms_r1.synth) 
isrms_r1.pred <- isrms_r1.dp$Y0plot %*% isrms_r1.synth$solution.w ####Y1s in Previous code
isrms_r1.alpha<- isrms_r1.dp$Y1plot - isrms_r1.pred

### 4. Placebo cases ##
### Run placebo for each donor unit
for (j in 1:length(donor))
{
    isrms_r1_pl.dp<- 
        dataprep( 
            foo = isrms1, 
            predictors = predictors,
            predictors.op = "mean",
            special.predictors = list(list("logsmilex_f2", 1968, c("mean"))),
            dependent = "logsmilex_f2", 
            unit.variable = "cown", 
            time.variable = "year",  
            treatment.identifier = donor[j], 
            controls.identifier = donor[-j], 
            time.predictors.prior = c(1956:1968), 
            time.optimize.ssr = c(1956:1968), 
            unit.names.variable = "country", 
            time.plot = c(1956:tr.end[i]) 
        ) 
    
    isrms_r1_synth_pl <- synth(isrms_r1_pl.dp)
    synthms2_pl.tables <- synth.tab( dataprep.res = isrms_r1_pl.dp, synth.res = isrms_r1_synth_pl, round.digit = 5) 
    assign(paste("W",treated[i],donor[j], sep="."), (synthms2_pl.tables$tab.w))
    assign(paste("V",treated[i],donor[j], sep="."), (isrms_r1_synth_pl$solution.v))
    assign(paste("Pred",treated[i],donor[j], sep="."), (synthms2_pl.tables$tab.pred))
    
    assign(paste("Y1",treated[i],donor[j], sep="."), (isrms_r1_pl.dp$Y1))
    assign(paste("Y1s",treated[i],donor[j], sep="."), (isrms_r1_pl.dp$Y0plot%*%isrms_r1_synth_pl$solution.w))
    assign(paste("Alpha",treated[i],donor[j], sep="."), (isrms_r1_pl.dp$Y1-isrms_r1_pl.dp$Y0plot%*%isrms_r1_synth_pl$solution.w))
}



# 5.1 e(avg_pre_rmspe_p) A measure of fit. Concerning if significant.
## The proportion of placebos that have a pre- treatment RMSPE at least as large as the average of the treated units. 

## pre-RMSPE comparison ##
isr.prermspe2 <- sqrt(mean(isrms_r1.alpha[1:13]^2, na.rm=TRUE))
isr.postrmspe2 <- sqrt(mean(isrms_r1.alpha[14:18]^2, na.rm=TRUE))
isr.prermspe.d2 <- numeric(length(donor))
isr.postrmspe.d2 <- numeric(length(donor))

for (j in 1:length(donor))
{
    temp <- get(paste("Alpha", treated[i], donor[j], sep="."))
    isr.prermspe.d2[j] <- sqrt(mean(temp[1:13]^2,na.rm=TRUE ))
    isr.postrmspe.d2[j] <- sqrt(mean(temp[14:18]^2,na.rm=TRUE ))
}

rmspe.fit2<- length(which(isr.prermspe.d2>isr.prermspe2))/length(donor)


## 5.2  e(pval_joint_post_t)
## The proportion of placebos that have a ratio of post-treatment RMSPE over pre-treatment RMSPE at least as large as the average ratio for the treated units.
pval.joint.post_t2 <- length(which(isr.postrmspe.d2/isr.prermspe.d2>isr.postrmspe2/isr.prermspe2))/length(donor)

alphas2 <-isrms_r1.alpha
for (j in 1: length(donor))
{
    temp12 <- get(paste("Alpha",treated[i],donor[j], sep="."))
    alphas2 <- cbind(alphas2, temp12)
    colnames(alphas2)[j+1] <- paste(treated[i],donor[j], sep=".")
}

## 5.3  e(pvals_t)
# A vector of the proportions of placebo pseudo t-statistics (unit’s effect divided by its pre-treatment RMSPE) that are at least as largeas the main pseudo t-statistic for each post-treatment period.
pvals_t2 <- numeric(5) 
for (t in 1:5)
{
    pvals_t2[t] <-  length(which(abs(alphas2[t+13,-1]/isr.prermspe.d2)>abs(alphas2[t+13,1]/isr.prermspe2)))/length(donor)
}

## 5.3.2 Plotting p values
plot(1:5, pvals_t2, ylim=c(0,1), axes = FALSE,xlab = "", ylab = "",  pch=20, bg="grey", lwd=5, cex=2)
box(); grid()
axis(1, at = 1:5, labels = seq(1, 5, 1), cex.axis = 1.25)
axis(2, at = seq(0, 1, 0.1), cex.axis = 1, las = 2)
title(xlab = "Number of years after treatment", cex.lab = 1.5, line = 2.4)
text(pvals_t2, labels=round(pvals_t2, 2), cex=0.9, font=2,  pos=3)
title(ylab = "Lead specific P-Values", cex.lab = 1.5, line = 2.6)



### 7. Getting country weights and balance tables
#country weights 
donorw <- synthms.tables$tab.w
#Balance Table 
synthms.tables$tab.pred

# PLOT 3 (isr_ms_mid3)
plot(1956:1973, isrms_r1.dp$Y1plot, ylim = c(3, 13),  type = "n", lwd = 2, xlab = "", ylab = "")
box(); grid()
abline(v = 1969, lwd = 1)  
lines(1956:1973, isrms_r1.dp$Y1plot, lwd = 2)
lines(1956:1973, isrms_r1.pred, lwd = 2, lty = 2)
title(xlab = "Year", cex.lab = 1.5, line = 2)
title(ylab = "Military Spending", cex.lab = 1, line = 2)
legend(list(x = 1956,y = 13.5), bty = "n", c("Israel", "Synthetic Israel"), lwd = 2, col = "black", lty = c(1, 2, -1), cex = 1.5, y.intersp=0.7)
text(1962.2, 5, "First year of Meir's tenure", cex = 1.7)
text(1971.3, 12.5, bquote(paste("MSPE" == .(round(isrms_r1.synth$loss.v, 3)))), cex = 1.5)
text(1971.4, 3.5, bquote(paste("p-value" == .(round(pval.joint.post_t2, 3)))), cex = 1.5)
arrows(1967, 5, 1969, 5, length = .1)
rug(isrms1$logsmilex_f2, side = 2, lwd=1, ticksize = .025)


### ISRAEL  Military Spending with previous_MID_avg   ####

### 2. Estimating Synthetic Counterfactual 

predictors = c( "polity2_20", "gdpgrowth_pw", "gdppc_pw", "previous_MID_avg", "women_legst", "cinc_c", "logtrade", "logthreat", "allies", "issues", "rival")
isrms_r1.dp <- dataprep(
    foo=isrms1, 
    predictors = predictors,
    special.predictors = list(list("logsmilex_f2", 1968, c("mean"))),
    predictors.op = "mean", 
    dependent = "logsmilex_f2", 
    unit.variable = "cown", 
    time.variable = "year",  
    treatment.identifier = 666, 
    controls.identifier = donor, 
    time.predictors.prior = c(1956:1968), 
    time.optimize.ssr = c(1956:1968), 
    unit.names.variable = "country", 
    time.plot = c(1956:tr.end[i]) 
)

isrms_r1.synth <- synth(isrms_r1.dp)
synthms.tables <- synth.tab( dataprep.res = isrms_r1.dp, synth.res = isrms_r1.synth) 
isrms_r1.pred <- isrms_r1.dp$Y0plot %*% isrms_r1.synth$solution.w ####Y1s in Previous code
isrms_r1.alpha<- isrms_r1.dp$Y1plot - isrms_r1.pred


### 4. Placebo cases ##
### Run placebo for each donor unit
for (j in 1:length(donor))
{
    isrms_r1_pl.dp<- 
        dataprep( 
            foo = isrms1, 
            predictors = predictors,
            predictors.op = "mean",
            special.predictors = list(list("logsmilex_f2", 1968, c("mean"))),
            dependent = "logsmilex_f2", 
            unit.variable = "cown", 
            time.variable = "year",  
            treatment.identifier = donor[j], 
            controls.identifier = donor[-j], 
            time.predictors.prior = c(1956:1968), 
            time.optimize.ssr = c(1956:1968), 
            unit.names.variable = "country", 
            time.plot = c(1956:tr.end[i]) 
        ) 
    
    isrms_r1_synth_pl <- synth(isrms_r1_pl.dp)
    synthms2_pl.tables <- synth.tab( dataprep.res = isrms_r1_pl.dp, synth.res = isrms_r1_synth_pl, round.digit = 5) 
    assign(paste("W",treated[i],donor[j], sep="."), (synthms2_pl.tables$tab.w))
    assign(paste("V",treated[i],donor[j], sep="."), (isrms_r1_synth_pl$solution.v))
    assign(paste("Pred",treated[i],donor[j], sep="."), (synthms2_pl.tables$tab.pred))
    
    assign(paste("Y1",treated[i],donor[j], sep="."), (isrms_r1_pl.dp$Y1))
    assign(paste("Y1s",treated[i],donor[j], sep="."), (isrms_r1_pl.dp$Y0plot%*%isrms_r1_synth_pl$solution.w))
    assign(paste("Alpha",treated[i],donor[j], sep="."), (isrms_r1_pl.dp$Y1-isrms_r1_pl.dp$Y0plot%*%isrms_r1_synth_pl$solution.w))
}



# 5.1 e(avg_pre_rmspe_p) A measure of fit. Concerning if significant.
## The proportion of placebos that have a pre- treatment RMSPE at least as large as the average of the treated units. 

## pre-RMSPE comparison ##
isr.prermspe2 <- sqrt(mean(isrms_r1.alpha[1:13]^2, na.rm=TRUE))
isr.postrmspe2 <- sqrt(mean(isrms_r1.alpha[14:18]^2, na.rm=TRUE))
isr.prermspe.d2 <- numeric(length(donor))
isr.postrmspe.d2 <- numeric(length(donor))

for (j in 1:length(donor))
{
    temp <- get(paste("Alpha", treated[i], donor[j], sep="."))
    isr.prermspe.d2[j] <- sqrt(mean(temp[1:13]^2,na.rm=TRUE ))
    isr.postrmspe.d2[j] <- sqrt(mean(temp[14:18]^2,na.rm=TRUE ))
}

rmspe.fit2<- length(which(isr.prermspe.d2>isr.prermspe2))/length(donor)


## 5.2  e(pval_joint_post_t)
## The proportion of placebos that have a ratio of post-treatment RMSPE over pre-treatment RMSPE at least as large as the average ratio for the treated units.
pval.joint.post_t2 <- length(which(isr.postrmspe.d2/isr.prermspe.d2>isr.postrmspe2/isr.prermspe2))/length(donor)

alphas2 <-isrms_r1.alpha
for (j in 1: length(donor))
{
    temp12 <- get(paste("Alpha",treated[i],donor[j], sep="."))
    alphas2 <- cbind(alphas2, temp12)
    colnames(alphas2)[j+1] <- paste(treated[i],donor[j], sep=".")
}

## 5.3  e(pvals_t)
# A vector of the proportions of placebo pseudo t-statistics (unit’s effect divided by its pre-treatment RMSPE) that are at least as largeas the main pseudo t-statistic for each post-treatment period.
pvals_t2 <- numeric(5) 
for (t in 1:5)
{
    pvals_t2[t] <-  length(which(abs(alphas2[t+13,-1]/isr.prermspe.d2)>abs(alphas2[t+13,1]/isr.prermspe2)))/length(donor)
}

## 5.3.2 Plotting p values
plot(1:5, pvals_t2, ylim=c(0,1), axes = FALSE,xlab = "", ylab = "",  pch=20, bg="grey", lwd=5, cex=2)
box(); grid()
axis(1, at = 1:5, labels = seq(1, 5, 1), cex.axis = 1.25)
axis(2, at = seq(0, 1, 0.1), cex.axis = 1, las = 2)
title(xlab = "Number of years after treatment", cex.lab = 1.5, line = 2.4)
text(pvals_t2, labels=round(pvals_t2, 2), cex=0.9, font=2,  pos=3)
title(ylab = "Lead specific P-Values", cex.lab = 1.5, line = 2.6)



### 7. Getting country weights and balance tables
#country weights 
donorw <- synthms.tables$tab.w
#Balance Table 
synthms.tables$tab.pred


# PLOT 2 (isr_ms_mid2)
plot(1956:1973, isrms_r1.dp$Y1plot, ylim = c(3, 13),  type = "n", lwd = 2, xlab = "", ylab = "")
box(); grid()
abline(v = 1969, lwd = 1)  
lines(1956:1973, isrms_r1.dp$Y1plot, lwd = 2)
lines(1956:1973, isrms_r1.pred, lwd = 2, lty = 2)
title(xlab = "Year", cex.lab = 1.5, line = 2)
title(ylab = "Military Spending", cex.lab = 1, line = 2)
legend(list(x = 1956,y = 13.5), bty = "n", c("Israel", "Synthetic Israel"), lwd = 2, col = "black", lty = c(1, 2, -1), cex = 1.5, y.intersp=0.7)
text(1962.2, 5, "First year of Meir's tenure", cex = 1.7)
text(1971.3, 12.5, bquote(paste("MSPE" == .(round(isrms_r1.synth$loss.v, 3)))), cex = 1.5)
text(1971.2, 3.5, bquote(paste("p-value" == .(round(pval.joint.post_t2, 3)))), cex = 1.5)
arrows(1967, 5, 1969, 5, length = .1)
rug(isrms1$logsmilex_f2, side = 2, lwd=1, ticksize = .025)

### ISRAEL  Military Spending with previous_fatal   ####

### 2. Estimating Synthetic Counterfactual 

predictors = c( "polity2_20", "gdpgrowth_pw", "gdppc_pw", "previous_fatal", "women_legst", "cinc_c", "logtrade", "logthreat", "allies", "issues", "rival")
isrms_r1.dp <- dataprep(
    foo=isrms1, 
    predictors = predictors,
    special.predictors = list(list("logsmilex_f2", 1968, c("mean"))),
    predictors.op = "mean", 
    dependent = "logsmilex_f2", 
    unit.variable = "cown", 
    time.variable = "year",  
    treatment.identifier = 666, 
    controls.identifier = donor, 
    time.predictors.prior = c(1956:1968), 
    time.optimize.ssr = c(1956:1968), 
    unit.names.variable = "country", 
    time.plot = c(1956:tr.end[i]) 
)

isrms_r1.synth <- synth(isrms_r1.dp)
synthms.tables <- synth.tab( dataprep.res = isrms_r1.dp, synth.res = isrms_r1.synth) 
isrms_r1.pred <- isrms_r1.dp$Y0plot %*% isrms_r1.synth$solution.w ####Y1s in Previous code
isrms_r1.alpha<- isrms_r1.dp$Y1plot - isrms_r1.pred



### 4. Placebo cases ##
### Run placebo for each donor unit
for (j in 1:length(donor))
{
    isrms_r1_pl.dp<- 
        dataprep( 
            foo = isrms1, 
            predictors = predictors,
            predictors.op = "mean",
            special.predictors = list(list("logsmilex_f2", 1968, c("mean"))),
            dependent = "logsmilex_f2", 
            unit.variable = "cown", 
            time.variable = "year",  
            treatment.identifier = donor[j], 
            controls.identifier = donor[-j], 
            time.predictors.prior = c(1956:1968), 
            time.optimize.ssr = c(1956:1968), 
            unit.names.variable = "country", 
            time.plot = c(1956:tr.end[i]) 
        ) 
    
    isrms_r1_synth_pl <- synth(isrms_r1_pl.dp)
    synthms2_pl.tables <- synth.tab( dataprep.res = isrms_r1_pl.dp, synth.res = isrms_r1_synth_pl, round.digit = 5) 
    assign(paste("W",treated[i],donor[j], sep="."), (synthms2_pl.tables$tab.w))
    assign(paste("V",treated[i],donor[j], sep="."), (isrms_r1_synth_pl$solution.v))
    assign(paste("Pred",treated[i],donor[j], sep="."), (synthms2_pl.tables$tab.pred))
    
    assign(paste("Y1",treated[i],donor[j], sep="."), (isrms_r1_pl.dp$Y1))
    assign(paste("Y1s",treated[i],donor[j], sep="."), (isrms_r1_pl.dp$Y0plot%*%isrms_r1_synth_pl$solution.w))
    assign(paste("Alpha",treated[i],donor[j], sep="."), (isrms_r1_pl.dp$Y1-isrms_r1_pl.dp$Y0plot%*%isrms_r1_synth_pl$solution.w))
}



# 5.1 e(avg_pre_rmspe_p) A measure of fit. Concerning if significant.
## The proportion of placebos that have a pre- treatment RMSPE at least as large as the average of the treated units. 

## pre-RMSPE comparison ##
isr.prermspe2 <- sqrt(mean(isrms_r1.alpha[1:13]^2, na.rm=TRUE))
isr.postrmspe2 <- sqrt(mean(isrms_r1.alpha[14:18]^2, na.rm=TRUE))
isr.prermspe.d2 <- numeric(length(donor))
isr.postrmspe.d2 <- numeric(length(donor))

for (j in 1:length(donor))
{
    temp <- get(paste("Alpha", treated[i], donor[j], sep="."))
    isr.prermspe.d2[j] <- sqrt(mean(temp[1:13]^2,na.rm=TRUE ))
    isr.postrmspe.d2[j] <- sqrt(mean(temp[14:18]^2,na.rm=TRUE ))
}

rmspe.fit2<- length(which(isr.prermspe.d2>isr.prermspe2))/length(donor)


## 5.2  e(pval_joint_post_t)
## The proportion of placebos that have a ratio of post-treatment RMSPE over pre-treatment RMSPE at least as large as the average ratio for the treated units.
pval.joint.post_t2 <- length(which(isr.postrmspe.d2/isr.prermspe.d2>isr.postrmspe2/isr.prermspe2))/length(donor)

alphas2 <-isrms_r1.alpha
for (j in 1: length(donor))
{
    temp12 <- get(paste("Alpha",treated[i],donor[j], sep="."))
    alphas2 <- cbind(alphas2, temp12)
    colnames(alphas2)[j+1] <- paste(treated[i],donor[j], sep=".")
}

## 5.3  e(pvals_t)
# A vector of the proportions of placebo pseudo t-statistics (unit’s effect divided by its pre-treatment RMSPE) that are at least as largeas the main pseudo t-statistic for each post-treatment period.
pvals_t2 <- numeric(5) 
for (t in 1:5)
{
    pvals_t2[t] <-  length(which(abs(alphas2[t+13,-1]/isr.prermspe.d2)>abs(alphas2[t+13,1]/isr.prermspe2)))/length(donor)
}

## 5.3.2 Plotting p values
plot(1:5, pvals_t2, ylim=c(0,1), axes = FALSE,xlab = "", ylab = "",  pch=20, bg="grey", lwd=5, cex=2)
box(); grid()
axis(1, at = 1:5, labels = seq(1, 5, 1), cex.axis = 1.25)
axis(2, at = seq(0, 1, 0.1), cex.axis = 1, las = 2)
title(xlab = "Number of years after treatment", cex.lab = 1.5, line = 2.4)
text(pvals_t2, labels=round(pvals_t2, 2), cex=0.9, font=2,  pos=3)
title(ylab = "Lead specific P-Values", cex.lab = 1.5, line = 2.6)


### 7. Getting country weights and balance tables
#country weights 
donorw <- synthms.tables$tab.w
#Balance Table 
synthms.tables$tab.pred

# PLOT 1 (isr_ms_mid1)
plot(1956:1973, isrms_r1.dp$Y1plot, ylim = c(3, 13),  type = "n", lwd = 2, xlab = "", ylab = "")
box(); grid()
abline(v = 1969, lwd = 1)  
lines(1956:1973, isrms_r1.dp$Y1plot, lwd = 2)
lines(1956:1973, isrms_r1.pred, lwd = 2, lty = 2)
title(xlab = "Year", cex.lab = 1.5, line = 2)
title(ylab = "Military Spending", cex.lab = 1, line = 2)
legend(list(x = 1956,y = 13.5), bty = "n", c("Israel", "Synthetic Israel"), lwd = 2, col = "black", lty = c(1, 2, -1), cex = 1.5, y.intersp=0.7)
text(1962.2, 5, "First year of Meir's tenure", cex = 1.7)
text(1971.3, 12.5, bquote(paste("MSPE" == .(round(isrms_r1.synth$loss.v, 3)))), cex = 1.5)
text(1971.2, 3.5, bquote(paste("p-value" == .(round(pval.joint.post_t2, 3)))), cex = 1.5)
arrows(1967, 5, 1969, 5, length = .1)
rug(isrms1$logsmilex_f2, side = 2, lwd=1, ticksize = .025)


#### #### UK  Military Spending with COUNT MID ####
rm(list=ls())
library(haven)
ukms1 <- read_dta("UK_ms.dta")
library(Synth)
library(tictoc)

colnames(ukms1)[4]=c("female")

# Merging with MID data
mid <- read_dta("mid_uk.dta")
ukms1 <- merge(ukms1, mid, all.x = TRUE )


### 1. Prepping data
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
ukms1 <- merge(ukms1, macro, by="cown")

predictors = c("polity2_20", "gdppc_pw", "gdpgrowth_pw", "count_mid", "women_legst", "cinc_c", "logtrade", "logthreat", "allies", "issues", "rival")

### 2. Estimating Synthetic Counterfactual for UK 
ukms.dp <- dataprep(
    foo=ukms1, 
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
    time.plot = c(1960:tr.end[i]) 
)
ukms.synth <- synth(ukms.dp)

synthms.tables <- synth.tab( dataprep.res = ukms.dp, synth.res = ukms.synth) 
ukms.pred <- ukms.dp$Y0plot %*% ukms.synth$solution.w ####Y1s in Previous code
ukms.alpha<- ukms.dp$Y1plot - ukms.pred


### 4. Placebo cases ##
### Run placebo for each donor unit
for (j in 1:length(donor))
{
    ukms_pl.dp<- 
        dataprep( 
            foo = ukms1, 
            predictors = predictors,
            predictors.op = "mean",
            special.predictors = list(list("logsmilex_f2",1978, c("mean"))),
            dependent = "logsmilex_f2", 
            unit.variable = "cown", 
            time.variable = "year",  
            treatment.identifier = donor[j], 
            controls.identifier = donor[-j], 
            time.predictors.prior = c(1960:1978), 
            time.optimize.ssr = c(1960:1978), 
            unit.names.variable = "country", 
            time.plot = c(1960:tr.end[i]) 
        ) 
    
    ukms_synth_pl <- synth(ukms_pl.dp)
    synthms_pl.tables <- synth.tab( dataprep.res = ukms_pl.dp, synth.res = ukms_synth_pl, round.digit = 5) 
    assign(paste("W",treated[i],donor[j], sep="."), (synthms_pl.tables$tab.w))
    assign(paste("V",treated[i],donor[j], sep="."), (ukms_synth_pl$solution.v))
    assign(paste("Pred",treated[i],donor[j], sep="."), (synthms_pl.tables$tab.pred))
    
    assign(paste("Y1",treated[i],donor[j], sep="."), (ukms_pl.dp$Y1))
    assign(paste("Y1s",treated[i],donor[j], sep="."), (ukms_pl.dp$Y0plot%*%ukms_synth_pl$solution.w))
    assign(paste("Alpha",treated[i],donor[j], sep="."), (ukms_pl.dp$Y1-ukms_pl.dp$Y0plot%*%ukms_synth_pl$solution.w))
}

### 5. Post Estimation results ###
## pre-RMSPE comparison ##
uk.prermspe <- sqrt(mean(ukms.alpha[1:19]^2, na.rm=TRUE))
uk.postrmspe <- sqrt(mean(ukms.alpha[20:31]^2, na.rm=TRUE))
uk.prermspe.d <- numeric(length(donor))
uk.postrmspe.d <- numeric(length(donor))

for (j in 1:length(donor))
{
    temp <- get(paste("Alpha", treated[i], donor[j], sep="."))
    uk.prermspe.d[j] <- sqrt(mean(temp[1:19]^2,na.rm=TRUE ))
    uk.postrmspe.d[j] <- sqrt(mean(temp[20:31]^2,na.rm=TRUE ))
}


# 5.1 e(avg_pre_rmspe_p) A measure of fit. Concerning if significant.
## The proportion of placebos that have a pre- treatment RMSPE at least as large as the average of the treated units. 
rmspe.fit <- length(which(uk.prermspe.d>uk.prermspe))/length(donor)


## 5.2  e(pval_joint_post_t)
## The proportion of placebos that have a ratio of post-treatment RMSPE over pre-treatment RMSPE at least as large as the average ratio for the treated units.
pval.joint.post_t <- length(which(uk.postrmspe.d/uk.prermspe.d>uk.postrmspe/uk.prermspe))/length(donor)

alphas <-ukms.alpha
for (j in 1: length(donor))
{
    temp1 <- get(paste("Alpha",treated[i],donor[j], sep="."))
    alphas <- cbind(alphas, temp1)
    colnames(alphas)[j+1] <- paste(treated[i],donor[j], sep=".")
}

## 5.3.1  e(pvals_t)
# A vector of the proportions of placebo pseudo t-statistics (unit’s effect divided by its pre-treatment RMSPE) that are at least as largeas the main pseudo t-statistic for each post-treatment period.
pvals_t <- numeric(12) 
for (t in 1:12)
{
    pvals_t[t] <-  length(which(abs(alphas[t+19,-1]/uk.prermspe.d)>abs(alphas[t+19,1]/uk.prermspe)))/length(donor)
}

#country weights =
donorw <- synthms.tables$tab.w
#Balance Table
synthms.tables$tab.pred

## PLOT 3. MID3 (uk_ms_mid3)
plot(1960:1990, ukms.dp$Y1plot, ylim = c(3, 13),  type = "n", lwd = 2, xlab = "", ylab = "")
box(); grid()
abline(v = 1979, lwd = 1)  
lines(1960:1990, ukms.dp$Y1plot, lwd = 2)
lines(1960:1990, ukms.pred, lwd = 2, lty = 2)
title(xlab = "Year", cex.lab = 1.5, line = 2)
title(ylab = "Military Spending", cex.lab = 2.5, line = 2)
legend(list(x = 1960,y = 13.5), bty = "n", c("UK", "Synthetic UK"), lwd = 2, col = "black", lty = c(1, 2, -1), cex = 1.3, y.intersp=0.5)
text(1969, 6, "First year of Thatcher's tenure", cex = 1.7)
arrows(1976, 6, 1979, 6, length = .2)
text(1985, 12.5, bquote(paste("MSPE" == .(round(ukms.synth$loss.v, 5)))), cex = 1.5)
text(1985, 3.5, bquote(paste("p-value" == .(round(pval.joint.post_t, 3)))), cex = 1.5)
rug(ukms1$logsmilex_f2, side = 2, lwd=1, ticksize = .025)




### UK  Military Spending with previous_MID_avg ####

predictors = c("polity2_20", "gdppc_pw", "gdpgrowth_pw", "previous_MID_avg", "women_legst", "cinc_c", "logtrade", "logthreat", "allies", "issues", "rival")

### 2. Estimating Synthetic Counterfactual for UK 
ukms.dp <- dataprep(
    foo=ukms1, 
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
    time.plot = c(1960:tr.end[i]) 
)
ukms.synth <- synth(ukms.dp)

synthms.tables <- synth.tab( dataprep.res = ukms.dp, synth.res = ukms.synth) 
ukms.pred <- ukms.dp$Y0plot %*% ukms.synth$solution.w ####Y1s in Previous code
ukms.alpha<- ukms.dp$Y1plot - ukms.pred



### 4. Placebo cases ##
### Run placebo for each donor unit
for (j in 1:length(donor))
{
    ukms_pl.dp<- 
        dataprep( 
            foo = ukms1, 
            predictors = predictors,
            predictors.op = "mean",
            special.predictors = list(list("logsmilex_f2",1978, c("mean"))),
            dependent = "logsmilex_f2", 
            unit.variable = "cown", 
            time.variable = "year",  
            treatment.identifier = donor[j], 
            controls.identifier = donor[-j], 
            time.predictors.prior = c(1960:1978), 
            time.optimize.ssr = c(1960:1978), 
            unit.names.variable = "country", 
            time.plot = c(1960:tr.end[i]) 
        ) 
    
    ukms_synth_pl <- synth(ukms_pl.dp)
    synthms_pl.tables <- synth.tab( dataprep.res = ukms_pl.dp, synth.res = ukms_synth_pl, round.digit = 5) 
    assign(paste("W",treated[i],donor[j], sep="."), (synthms_pl.tables$tab.w))
    assign(paste("V",treated[i],donor[j], sep="."), (ukms_synth_pl$solution.v))
    assign(paste("Pred",treated[i],donor[j], sep="."), (synthms_pl.tables$tab.pred))
    
    assign(paste("Y1",treated[i],donor[j], sep="."), (ukms_pl.dp$Y1))
    assign(paste("Y1s",treated[i],donor[j], sep="."), (ukms_pl.dp$Y0plot%*%ukms_synth_pl$solution.w))
    assign(paste("Alpha",treated[i],donor[j], sep="."), (ukms_pl.dp$Y1-ukms_pl.dp$Y0plot%*%ukms_synth_pl$solution.w))
}

### 5. Post Estimation results ###
## pre-RMSPE comparison ##
uk.prermspe <- sqrt(mean(ukms.alpha[1:19]^2, na.rm=TRUE))
uk.postrmspe <- sqrt(mean(ukms.alpha[20:31]^2, na.rm=TRUE))
uk.prermspe.d <- numeric(length(donor))
uk.postrmspe.d <- numeric(length(donor))

for (j in 1:length(donor))
{
    temp <- get(paste("Alpha", treated[i], donor[j], sep="."))
    uk.prermspe.d[j] <- sqrt(mean(temp[1:19]^2,na.rm=TRUE ))
    uk.postrmspe.d[j] <- sqrt(mean(temp[20:31]^2,na.rm=TRUE ))
}


# 5.1 e(avg_pre_rmspe_p) A measure of fit. Concerning if significant.
## The proportion of placebos that have a pre- treatment RMSPE at least as large as the average of the treated units. 
rmspe.fit <- length(which(uk.prermspe.d>uk.prermspe))/length(donor)


## 5.2  e(pval_joint_post_t)
## The proportion of placebos that have a ratio of post-treatment RMSPE over pre-treatment RMSPE at least as large as the average ratio for the treated units.
pval.joint.post_t <- length(which(uk.postrmspe.d/uk.prermspe.d>uk.postrmspe/uk.prermspe))/length(donor)

alphas <-ukms.alpha
for (j in 1: length(donor))
{
    temp1 <- get(paste("Alpha",treated[i],donor[j], sep="."))
    alphas <- cbind(alphas, temp1)
    colnames(alphas)[j+1] <- paste(treated[i],donor[j], sep=".")
}

## 5.3.1  e(pvals_t)
# A vector of the proportions of placebo pseudo t-statistics (unit’s effect divided by its pre-treatment RMSPE) that are at least as largeas the main pseudo t-statistic for each post-treatment period.
pvals_t <- numeric(12) 
for (t in 1:12)
{
    pvals_t[t] <-  length(which(abs(alphas[t+19,-1]/uk.prermspe.d)>abs(alphas[t+19,1]/uk.prermspe)))/length(donor)
}

#country weights =
donorw <- synthms.tables$tab.w
#Balance Table
synthms.tables$tab.pred

#PLOT 2 MID2 (uk_ms_mid2)
plot(1960:1990, ukms.dp$Y1plot, ylim = c(3, 13),  type = "n", lwd = 2, xlab = "", ylab = "")
box(); grid()
abline(v = 1979, lwd = 1)  
lines(1960:1990, ukms.dp$Y1plot, lwd = 2)
lines(1960:1990, ukms.pred, lwd = 2, lty = 2)
title(xlab = "Year", cex.lab = 1.5, line = 2)
title(ylab = "Military Spending", cex.lab = 2.5, line = 2)
legend(list(x = 1960,y = 13.5), bty = "n", c("UK", "Synthetic UK"), lwd = 2, col = "black", lty = c(1, 2, -1), cex = 1.3, y.intersp=0.5)
text(1969, 6, "First year of Thatcher's tenure", cex = 1.7)
arrows(1976, 6, 1979, 6, length = .2)
text(1985, 12.5, bquote(paste("MSPE" == .(round(ukms.synth$loss.v, 5)))), cex = 1.5)
text(1985, 3.5, bquote(paste("p-value" == .(round(pval.joint.post_t, 3)))), cex = 1.5)
rug(ukms1$logsmilex_f2, side = 2, lwd=1, ticksize = .025)


### UK  Military Spending with previous_fatal   ####

predictors = c("polity2_20", "gdppc_pw", "gdpgrowth_pw", "previous_fatal", "women_legst", "cinc_c", "logtrade", "logthreat", "allies", "issues", "rival")

### 2. Estimating Synthetic Counterfactual for UK 
ukms.dp <- dataprep(
    foo=ukms1, 
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
    time.plot = c(1960:tr.end[i]) 
)
ukms.synth <- synth(ukms.dp)

synthms.tables <- synth.tab( dataprep.res = ukms.dp, synth.res = ukms.synth) 
ukms.pred <- ukms.dp$Y0plot %*% ukms.synth$solution.w ####Y1s in Previous code
ukms.alpha<- ukms.dp$Y1plot - ukms.pred


### 4. Placebo cases ##
### Run placebo for each donor unit
for (j in 1:length(donor))
{
    
    ukms_pl.dp<- 
        dataprep( 
            foo = ukms1, 
            predictors = predictors,
            predictors.op = "mean",
            special.predictors = list(list("logsmilex_f2",1978, c("mean"))),
            dependent = "logsmilex_f2", 
            unit.variable = "cown", 
            time.variable = "year",  
            treatment.identifier = donor[j], 
            controls.identifier = donor[-j], 
            time.predictors.prior = c(1960:1978), 
            time.optimize.ssr = c(1960:1978), 
            unit.names.variable = "country", 
            time.plot = c(1960:tr.end[i]) 
        ) 
    
    ukms_synth_pl <- synth(ukms_pl.dp)
    synthms_pl.tables <- synth.tab( dataprep.res = ukms_pl.dp, synth.res = ukms_synth_pl, round.digit = 5) 
    assign(paste("W",treated[i],donor[j], sep="."), (synthms_pl.tables$tab.w))
    assign(paste("V",treated[i],donor[j], sep="."), (ukms_synth_pl$solution.v))
    assign(paste("Pred",treated[i],donor[j], sep="."), (synthms_pl.tables$tab.pred))
    
    assign(paste("Y1",treated[i],donor[j], sep="."), (ukms_pl.dp$Y1))
    assign(paste("Y1s",treated[i],donor[j], sep="."), (ukms_pl.dp$Y0plot%*%ukms_synth_pl$solution.w))
    assign(paste("Alpha",treated[i],donor[j], sep="."), (ukms_pl.dp$Y1-ukms_pl.dp$Y0plot%*%ukms_synth_pl$solution.w))
}


### 5. Post Estimation results ###
## pre-RMSPE comparison ##
uk.prermspe <- sqrt(mean(ukms.alpha[1:19]^2, na.rm=TRUE))
uk.postrmspe <- sqrt(mean(ukms.alpha[20:31]^2, na.rm=TRUE))
uk.prermspe.d <- numeric(length(donor))
uk.postrmspe.d <- numeric(length(donor))

for (j in 1:length(donor))
{
    temp <- get(paste("Alpha", treated[i], donor[j], sep="."))
    uk.prermspe.d[j] <- sqrt(mean(temp[1:19]^2,na.rm=TRUE ))
    uk.postrmspe.d[j] <- sqrt(mean(temp[20:31]^2,na.rm=TRUE ))
}


# 5.1 e(avg_pre_rmspe_p) A measure of fit. Concerning if significant.
## The proportion of placebos that have a pre- treatment RMSPE at least as large as the average of the treated units. 
rmspe.fit <- length(which(uk.prermspe.d>uk.prermspe))/length(donor)


## 5.2  e(pval_joint_post_t)
## The proportion of placebos that have a ratio of post-treatment RMSPE over pre-treatment RMSPE at least as large as the average ratio for the treated units.
pval.joint.post_t <- length(which(uk.postrmspe.d/uk.prermspe.d>uk.postrmspe/uk.prermspe))/length(donor)

alphas <-ukms.alpha
for (j in 1: length(donor))
{
    temp1 <- get(paste("Alpha",treated[i],donor[j], sep="."))
    alphas <- cbind(alphas, temp1)
    colnames(alphas)[j+1] <- paste(treated[i],donor[j], sep=".")
}

## 5.3.1  e(pvals_t)
# A vector of the proportions of placebo pseudo t-statistics (unit’s effect divided by its pre-treatment RMSPE) that are at least as largeas the main pseudo t-statistic for each post-treatment period.
pvals_t <- numeric(12) 
for (t in 1:12)
{
    pvals_t[t] <-  length(which(abs(alphas[t+19,-1]/uk.prermspe.d)>abs(alphas[t+19,1]/uk.prermspe)))/length(donor)
}

#country weights =
donorw <- synthms.tables$tab.w
#Balance Table
synthms.tables$tab.pred

## PLOT 1. MID1 (uk_ms_mid1)
plot(1960:1990, ukms.dp$Y1plot, ylim = c(3, 13),  type = "n", lwd = 2, xlab = "", ylab = "")
box(); grid()
abline(v = 1979, lwd = 1)  
lines(1960:1990, ukms.dp$Y1plot, lwd = 2)
lines(1960:1990, ukms.pred, lwd = 2, lty = 2)
title(xlab = "Year", cex.lab = 1.5, line = 2)
title(ylab = "Military Spending", cex.lab = 2.5, line = 2)
legend(list(x = 1960,y = 13.5), bty = "n", c("UK", "Synthetic UK"), lwd = 2, col = "black", lty = c(1, 2, -1), cex = 1.3, y.intersp=0.5)
text(1969, 6, "First year of Thatcher's tenure", cex = 1.7)
arrows(1976, 6, 1979, 6, length = .2)
text(1985, 12.5, bquote(paste("MSPE" == .(round(ukms.synth$loss.v, 5)))), cex = 1.5)
text(1985, 3.5, bquote(paste("p-value" == .(round(pval.joint.post_t, 3)))), cex = 1.5)
rug(ukms1$logsmilex_f2, side = 2, lwd=1, ticksize = .025)


#



#### #### INDIA  Military Spending with COUNT MID ####
rm(list=ls())
library(haven)
library(Synth)
library(tictoc)

indms1 <- read_dta("India_ms56-84_725.dta")
indms1 <- indms1[indms1$year<1977,]
colnames(indms1)[4]=c("female")

mid <- read_dta("mid_ind.dta")
indms1 <- merge(indms1, mid, all.x = TRUE )

### 1. Prepping data
macro <- aggregate(female ~ cown, data = indms1, sum)
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
    tr.years <- indms1$year[(indms1$female==1 & indms1$cown==treated[i])]
    tr.start[i] <- min(tr.years)
    tr.end[i] <- max(tr.years)
}
indms1 <- merge(indms1, macro, by="cown")

predictors = c("polity2_20", "gdppc_pw", "gdpgrowth_pw", "count_mid", "women_legst", "cinc_c", "logtrade", "logthreat_f2", "allies", "issues", "rival")

### 2. Estimating Synthetic Counterfactual for INDIA
indms.dp <- dataprep(
    foo=indms1, 
    predictors = predictors,
    special.predictors = list(list("logsmilex_f2",1965, c("mean"))),
    predictors.op = "mean", 
    dependent = "logsmilex_f2", 
    unit.variable = "cown", 
    time.variable = "year",  
    treatment.identifier = 750, 
    controls.identifier = donor, 
    time.predictors.prior = c(1956:1965), 
    time.optimize.ssr = c(1956:1965), 
    unit.names.variable = "country", 
    time.plot = c(1956:tr.end[i]) 
)
indms.synth <- synth(indms.dp, "ALL")

synthms.tables <- synth.tab( dataprep.res = indms.dp, synth.res = indms.synth) 
indms.pred <- indms.dp$Y0plot %*% indms.synth$solution.w ####Y1s in Previous code
indms.alpha<- indms.dp$Y1plot - indms.pred


### 4. Placebo cases ##
### Run placebo for each donor unit

for (j in 1:length(donor))
{
    indms_pl.dp<- 
        dataprep( 
            foo = indms1, 
            predictors = predictors,
            predictors.op = "mean",
            special.predictors = list(list("logsmilex_f2",1965, c("mean"))),
            dependent = "logsmilex_f2", 
            unit.variable = "cown", 
            time.variable = "year",  
            treatment.identifier = donor[j], 
            controls.identifier = donor[-j], 
            time.predictors.prior = c(1956:1965), 
            time.optimize.ssr = c(1956:1965), 
            unit.names.variable = "country", 
            time.plot = c(1956:tr.end[i]) 
        ) 
    
    indms_synth_pl <- synth(indms_pl.dp)
    synthms_pl.tables <- synth.tab( dataprep.res = indms_pl.dp, synth.res = indms_synth_pl, round.digit = 5) 
    assign(paste("W",treated[i],donor[j], sep="."), (synthms_pl.tables$tab.w))
    assign(paste("V",treated[i],donor[j], sep="."), (indms_synth_pl$solution.v))
    assign(paste("Pred",treated[i],donor[j], sep="."), (synthms_pl.tables$tab.pred))
    
    assign(paste("Y1",treated[i],donor[j], sep="."), (indms_pl.dp$Y1))
    assign(paste("Y1s",treated[i],donor[j], sep="."), (indms_pl.dp$Y0plot%*%indms_synth_pl$solution.w))
    assign(paste("Alpha",treated[i],donor[j], sep="."), (indms_pl.dp$Y1-indms_pl.dp$Y0plot%*%indms_synth_pl$solution.w))
    #bad.fit <- numeric(33)
    #bad.fit[j] <- as.numeric(ifelse(indms_synth_pl$loss.v > indms.synth$loss.v*5, 1, 0))
    
}


### 5. Post Estimation results ###
## pre-RMSPE comparison ##
ind.prermspe <- sqrt(mean(indms.alpha[1:10]^2, na.rm=TRUE))
ind.postrmspe <- sqrt(mean(indms.alpha[11:21]^2, na.rm=TRUE))
ind.prermspe.d <- numeric(length(donor))
ind.postrmspe.d <- numeric(length(donor))

for (j in 1:length(donor))
{
    temp <- get(paste("Alpha", treated[i], donor[j], sep="."))
    ind.prermspe.d[j] <- sqrt(mean(temp[1:10]^2,na.rm=TRUE ))
    ind.postrmspe.d[j] <- sqrt(mean(temp[11:21]^2,na.rm=TRUE ))
}

# 5.1 e(avg_pre_rmspe_p) A measure of fit. Concerning if significant.
## The proportion of placebos that have a pre- treatment RMSPE at least as large as the average of the treated units. 
rmspe.fit <- length(which(ind.prermspe.d>ind.prermspe))/length(donor)

## 5.2  e(pval_joint_post_t)
## The proportion of placebos that have a ratio of post-treatment RMSPE over pre-treatment RMSPE at least as large as the average ratio for the treated units.
pval.joint.post_t <- length(which(ind.postrmspe.d/ind.prermspe.d>ind.postrmspe/ind.prermspe))/length(donor)

alphas <-indms.alpha
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
    pvals_t[t] <-  length(which(abs(alphas[t+10,-1]/ind.prermspe.d)>abs(alphas[t+10,1]/ind.prermspe)))/length(donor)
}

#country weights Table
synthms.tables$tab.w
#Balance Table
synthms.tables$tab.pred

###. PLOT 3 (ind_ms_mid3)
plot(1956:1976, indms.dp$Y1plot, ylim = c(3, 13),  type = "n", lwd = 2, xlab = "", ylab = "")
box(); grid()
abline(v = 1966, lwd = 1)  
lines(1956:1976, indms.dp$Y1plot, lwd = 2)
lines(1956:1976, indms.pred, lwd = 2, lty = 2)
title(xlab = "Year", cex.lab = 1.5, line = 2)
title(ylab = "Military Spending", cex.lab = 1.5, line = 2)
legend(list(x = 1956,y = 13), bty = "n", c("India", "Synthetic India"), lwd = 2, col = "black", lty = c(1, 2, -1), cex = 1.3, y.intersp=0.6)
text(1960.5, 5, "First year of Gandhi's tenure", cex = 1.3)
arrows(1964.8, 5, 1966, 5, length = .2)
text(1973, 12.5, bquote(paste("MSPE" == .(round(indms.synth$loss.v, 3)))), cex = 1.5)
text(1973, 4.5, bquote(paste("p-value" == .(round(pval.joint.post_t, 3)))), cex = 1.5)
rug(indms1$logsmilex_f2, side = 2, lwd=1, ticksize = .025)



### India  Military Spending with previous_MID_avg ####

predictors = c("polity2_20", "gdppc_pw", "gdpgrowth_pw", "previous_MID_avg", "women_legst", "cinc_c", "logtrade", "logthreat_f2", "allies", "issues", "rival")

### 2. Estimating Synthetic Counterfactual for INDIA
indms.dp <- dataprep(
    foo=indms1, 
    predictors = predictors,
    special.predictors = list(list("logsmilex_f2",1965, c("mean"))),
    predictors.op = "mean", 
    dependent = "logsmilex_f2", 
    unit.variable = "cown", 
    time.variable = "year",  
    treatment.identifier = 750, 
    controls.identifier = donor, 
    time.predictors.prior = c(1956:1965), 
    time.optimize.ssr = c(1956:1965), 
    unit.names.variable = "country", 
    time.plot = c(1956:tr.end[i]) 
)
indms.synth <- synth(indms.dp, "ALL")

synthms.tables <- synth.tab( dataprep.res = indms.dp, synth.res = indms.synth) 
indms.pred <- indms.dp$Y0plot %*% indms.synth$solution.w ####Y1s in Previous code
indms.alpha<- indms.dp$Y1plot - indms.pred


### 4. Placebo cases ##
### Run placebo for each donor unit

for (j in 1:length(donor))
{
    indms_pl.dp<- 
        dataprep( 
            foo = indms1, 
            predictors = predictors,
            predictors.op = "mean",
            special.predictors = list(list("logsmilex_f2",1965, c("mean"))),
            dependent = "logsmilex_f2", 
            unit.variable = "cown", 
            time.variable = "year",  
            treatment.identifier = donor[j], 
            controls.identifier = donor[-j], 
            time.predictors.prior = c(1956:1965), 
            time.optimize.ssr = c(1956:1965), 
            unit.names.variable = "country", 
            time.plot = c(1956:tr.end[i]) 
        ) 
    
    indms_synth_pl <- synth(indms_pl.dp)
    synthms_pl.tables <- synth.tab( dataprep.res = indms_pl.dp, synth.res = indms_synth_pl, round.digit = 5) 
    assign(paste("W",treated[i],donor[j], sep="."), (synthms_pl.tables$tab.w))
    assign(paste("V",treated[i],donor[j], sep="."), (indms_synth_pl$solution.v))
    assign(paste("Pred",treated[i],donor[j], sep="."), (synthms_pl.tables$tab.pred))
    
    assign(paste("Y1",treated[i],donor[j], sep="."), (indms_pl.dp$Y1))
    assign(paste("Y1s",treated[i],donor[j], sep="."), (indms_pl.dp$Y0plot%*%indms_synth_pl$solution.w))
    assign(paste("Alpha",treated[i],donor[j], sep="."), (indms_pl.dp$Y1-indms_pl.dp$Y0plot%*%indms_synth_pl$solution.w))
    #bad.fit <- numeric(33)
    #bad.fit[j] <- as.numeric(ifelse(indms_synth_pl$loss.v > indms.synth$loss.v*5, 1, 0))
    
}

### 5. Post Estimation results ###
## pre-RMSPE comparison ##
ind.prermspe <- sqrt(mean(indms.alpha[1:10]^2, na.rm=TRUE))
ind.postrmspe <- sqrt(mean(indms.alpha[11:21]^2, na.rm=TRUE))
ind.prermspe.d <- numeric(length(donor))
ind.postrmspe.d <- numeric(length(donor))

for (j in 1:length(donor))
{
    temp <- get(paste("Alpha", treated[i], donor[j], sep="."))
    ind.prermspe.d[j] <- sqrt(mean(temp[1:10]^2,na.rm=TRUE ))
    ind.postrmspe.d[j] <- sqrt(mean(temp[11:21]^2,na.rm=TRUE ))
}

# 5.1 e(avg_pre_rmspe_p) A measure of fit. Concerning if significant.
## The proportion of placebos that have a pre- treatment RMSPE at least as large as the average of the treated units. 
rmspe.fit <- length(which(ind.prermspe.d>ind.prermspe))/length(donor)


## 5.2  e(pval_joint_post_t)
## The proportion of placebos that have a ratio of post-treatment RMSPE over pre-treatment RMSPE at least as large as the average ratio for the treated units.
pval.joint.post_t <- length(which(ind.postrmspe.d/ind.prermspe.d>ind.postrmspe/ind.prermspe))/length(donor)

alphas <-indms.alpha
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
    pvals_t[t] <-  length(which(abs(alphas[t+10,-1]/ind.prermspe.d)>abs(alphas[t+10,1]/ind.prermspe)))/length(donor)
}

#country weights Table
synthms.tables$tab.w
#Balance Table
synthms.tables$tab.pred

###. PLOT 2 (ind_ms_mid2)
plot(1956:1976, indms.dp$Y1plot, ylim = c(3, 13),  type = "n", lwd = 2, xlab = "", ylab = "")
box(); grid()
abline(v = 1966, lwd = 1)  
lines(1956:1976, indms.dp$Y1plot, lwd = 2)
lines(1956:1976, indms.pred, lwd = 2, lty = 2)
title(xlab = "Year", cex.lab = 1.5, line = 2)
title(ylab = "Military Spending", cex.lab = 1.5, line = 2)
legend(list(x = 1956,y = 13), bty = "n", c("India", "Synthetic India"), lwd = 2, col = "black", lty = c(1, 2, -1), cex = 1.3, y.intersp=0.6)
text(1960.5, 5, "First year of Gandhi's tenure", cex = 1.3)
arrows(1964.8, 5, 1966, 5, length = .2)
text(1973, 12.5, bquote(paste("MSPE" == .(round(indms.synth$loss.v, 3)))), cex = 1.5)
text(1973, 4.5, bquote(paste("p-value" == .(round(pval.joint.post_t, 3)))), cex = 1.5)
rug(indms1$logsmilex_f2, side = 2, lwd=1, ticksize = .025)

### India  Military Spending with previous_fatal ####  

predictors = c("polity2_20", "gdppc_pw", "gdpgrowth_pw", "previous_fatal", "women_legst", "cinc_c", "logtrade", "logthreat_f2", "allies", "issues", "rival")

### 2. Estimating Synthetic Counterfactual for INDIA
indms.dp <- dataprep(
    foo=indms1, 
    predictors = predictors,
    special.predictors = list(list("logsmilex_f2",1965, c("mean"))),
    predictors.op = "mean", 
    dependent = "logsmilex_f2", 
    unit.variable = "cown", 
    time.variable = "year",  
    treatment.identifier = 750, 
    controls.identifier = donor, 
    time.predictors.prior = c(1956:1965), 
    time.optimize.ssr = c(1956:1965), 
    unit.names.variable = "country", 
    time.plot = c(1956:tr.end[i]) 
)
indms.synth <- synth(indms.dp, "ALL")

synthms.tables <- synth.tab( dataprep.res = indms.dp, synth.res = indms.synth) 
indms.pred <- indms.dp$Y0plot %*% indms.synth$solution.w ####Y1s in Previous code
indms.alpha<- indms.dp$Y1plot - indms.pred


### 4. Placebo cases ##
### Run placebo for each donor unit

for (j in 1:length(donor))
{
    indms_pl.dp<- 
        dataprep( 
            foo = indms1, 
            predictors = predictors,
            predictors.op = "mean",
            special.predictors = list(list("logsmilex_f2",1965, c("mean"))),
            dependent = "logsmilex_f2", 
            unit.variable = "cown", 
            time.variable = "year",  
            treatment.identifier = donor[j], 
            controls.identifier = donor[-j], 
            time.predictors.prior = c(1956:1965), 
            time.optimize.ssr = c(1956:1965), 
            unit.names.variable = "country", 
            time.plot = c(1956:tr.end[i]) 
        ) 
    
    indms_synth_pl <- synth(indms_pl.dp)
    synthms_pl.tables <- synth.tab( dataprep.res = indms_pl.dp, synth.res = indms_synth_pl, round.digit = 5) 
    assign(paste("W",treated[i],donor[j], sep="."), (synthms_pl.tables$tab.w))
    assign(paste("V",treated[i],donor[j], sep="."), (indms_synth_pl$solution.v))
    assign(paste("Pred",treated[i],donor[j], sep="."), (synthms_pl.tables$tab.pred))
    
    assign(paste("Y1",treated[i],donor[j], sep="."), (indms_pl.dp$Y1))
    assign(paste("Y1s",treated[i],donor[j], sep="."), (indms_pl.dp$Y0plot%*%indms_synth_pl$solution.w))
    assign(paste("Alpha",treated[i],donor[j], sep="."), (indms_pl.dp$Y1-indms_pl.dp$Y0plot%*%indms_synth_pl$solution.w))
    #bad.fit <- numeric(33)
    #bad.fit[j] <- as.numeric(ifelse(indms_synth_pl$loss.v > indms.synth$loss.v*5, 1, 0))
    
}

### 5. Post Estimation results ###
## pre-RMSPE comparison ##
ind.prermspe <- sqrt(mean(indms.alpha[1:10]^2, na.rm=TRUE))
ind.postrmspe <- sqrt(mean(indms.alpha[11:21]^2, na.rm=TRUE))
ind.prermspe.d <- numeric(length(donor))
ind.postrmspe.d <- numeric(length(donor))

for (j in 1:length(donor))
{
    temp <- get(paste("Alpha", treated[i], donor[j], sep="."))
    ind.prermspe.d[j] <- sqrt(mean(temp[1:10]^2,na.rm=TRUE ))
    ind.postrmspe.d[j] <- sqrt(mean(temp[11:21]^2,na.rm=TRUE ))
}

# 5.1 e(avg_pre_rmspe_p) A measure of fit. Concerning if significant.
## The proportion of placebos that have a pre- treatment RMSPE at least as large as the average of the treated units. 
rmspe.fit <- length(which(ind.prermspe.d>ind.prermspe))/length(donor)


## 5.2  e(pval_joint_post_t)
## The proportion of placebos that have a ratio of post-treatment RMSPE over pre-treatment RMSPE at least as large as the average ratio for the treated units.
pval.joint.post_t <- length(which(ind.postrmspe.d/ind.prermspe.d>ind.postrmspe/ind.prermspe))/length(donor)

alphas <-indms.alpha
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
    pvals_t[t] <-  length(which(abs(alphas[t+10,-1]/ind.prermspe.d)>abs(alphas[t+10,1]/ind.prermspe)))/length(donor)
}

#country weights Table
synthms.tables$tab.w
#Balance Table
synthms.tables$tab.pred

###. PLOT 1 (ind_ms_mid1)
plot(1956:1976, indms.dp$Y1plot, ylim = c(3, 13),  type = "n", lwd = 2, xlab = "", ylab = "")
box(); grid()
abline(v = 1966, lwd = 1)  
lines(1956:1976, indms.dp$Y1plot, lwd = 2)
lines(1956:1976, indms.pred, lwd = 2, lty = 2)
title(xlab = "Year", cex.lab = 1.5, line = 2)
title(ylab = "Military Spending", cex.lab = 1.5, line = 2)
legend(list(x = 1956,y = 13), bty = "n", c("India", "Synthetic India"), lwd = 2, col = "black", lty = c(1, 2, -1), cex = 1.3, y.intersp=0.6)
text(1960.5, 5, "First year of Gandhi's tenure", cex = 1.3)
arrows(1964.8, 5, 1966, 5, length = .2)
text(1973, 12.5, bquote(paste("MSPE" == .(round(indms.synth$loss.v, 3)))), cex = 1.5)
text(1973, 4.5, bquote(paste("p-value" == .(round(0.999, 5)))), cex = 1.5)
rug(indms1$logsmilex_f2, side = 2, lwd=1, ticksize = .025)

