#### ULKAR IMAMVERDI AND PATRICK SHEA    ###
###    Female Leaders and Foreign Policy ###
### Part 1: SINGLE TREATMENT ANALYSIS    ###
############################################

########SYNTHETIC CASE STUDIES: ISRAEL  Military Spending ########
rm(list=ls())
library(haven)
isrms1 <- read_dta("isr53-73fem_ms504r.dta")
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

### 2. Estimating Synthetic Counterfactual for Israel 

predictors2 = c( "polity2_20", "gdpgrowth_pw", "gdppc_pw", "women_legst", "cinc_c", "logtrade", "logthreat", "allies", "issues", "rival")
isrms2.dp <- dataprep(
    foo=isrms1, 
    predictors = predictors2,
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
isrms2.synth <- synth(isrms2.dp)


synthms2.tables <- synth.tab( dataprep.res = isrms2.dp, synth.res = isrms2.synth) 
isrms2.pred <- isrms2.dp$Y0plot %*% isrms2.synth$solution.w ####Y1s in Previous code
isrms2.alpha<- isrms2.dp$Y1plot - isrms2.pred

## 3. Plot results
plot(1956:1973, isrms2.dp$Y1plot, ylim = c(3, 13),  type = "n", lwd = 2, xlab = "", ylab = "")
box(); grid()
abline(v = 1969, lwd = 1)  
lines(1956:1973, isrms2.dp$Y1plot, lwd = 2)
lines(1956:1973, isrms2.pred, lwd = 2, lty = 2)
title(xlab = "Year", cex.lab = 1.5, line = 2)
title(ylab = "Military Spending", cex.lab = 1, line = 2)
legend(list(x = 1956,y = 15.5), bty = "n", c("Israel", "Synthetic Israel"), lwd = 2, col = "black", lty = c(1, 2, -1), cex = 1.8, y.intersp=0.3)
text(1962.8, 5, "First year of Meir's tenure", cex = 2.3)
text(1971.3, 12.5, bquote(paste("MSPE" == .(round(isrms2.synth$loss.v, 5)))), cex = 1.9)
arrows(1967, 5, 1969, 5, length = .1)
rug(isrms1$logsmilex_f2, side = 2, lwd=1, ticksize = .025)

### 4. Placebo cases ##
### Run placebo for each donor unit
for (j in 1:length(donor))
{
    isrms2_pl.dp<- 
        dataprep( 
            foo = isrms1, 
            predictors = predictors2,
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
    
    isrms2_synth_pl <- synth(isrms2_pl.dp)
    synthms2_pl.tables <- synth.tab( dataprep.res = isrms2_pl.dp, synth.res = isrms2_synth_pl, round.digit = 5) 
    assign(paste("W",treated[i],donor[j], sep="."), (synthms2_pl.tables$tab.w))
    assign(paste("V",treated[i],donor[j], sep="."), (isrms2_synth_pl$solution.v))
    assign(paste("Pred",treated[i],donor[j], sep="."), (synthms2_pl.tables$tab.pred))
    
    assign(paste("Y1",treated[i],donor[j], sep="."), (isrms2_pl.dp$Y1))
    assign(paste("Y1s",treated[i],donor[j], sep="."), (isrms2_pl.dp$Y0plot%*%isrms2_synth_pl$solution.w))
    assign(paste("Alpha",treated[i],donor[j], sep="."), (isrms2_pl.dp$Y1-isrms2_pl.dp$Y0plot%*%isrms2_synth_pl$solution.w))
}


## pre-RMSPE comparison ##
isr.prermspe2 <- sqrt(mean(isrms2.alpha[1:13]^2, na.rm=TRUE))
isr.postrmspe2 <- sqrt(mean(isrms2.alpha[14:18]^2, na.rm=TRUE))
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

alphas2 <-isrms2.alpha
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
text(pvals_t2, labels = round(pvals_t2, 2), pos = 3, cex = 1.5)
title(xlab = "Number of years after treatment", cex.lab = 1.5, line = 2.4)
title(ylab = "Lead specific P-Values", cex.lab = 2.2, line = 2.6)



### PLot Plabeco and Synthetic Israel
plot(1956:1973, isrms2.alpha, type = "n", ylim = c(-3, 3), lwd = 2, xlab = "", ylab = "", axes = FALSE) 
box(); grid()
#for(j in 1: 24){
#    lines(1956:1973, alphas2[, j], type = "l", col = "gray80", lwd = 2)
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
text(1964.5, -2.5, "Treatment time", cex = 1.5)
arrows(1967, -2.5, 1969, -2.5, length = .1)
legend(list(x = 1956,y = 3.5), ncol = 1, bty = "n", c("Israel", "Placebo Countries"), lty=c(1,1), lwd = 2, col = c("black", "gray80"), cex = 1.25, y.intersp=0.5)

# 7. Tables 
#country weights Table B14
synthms2.tables$tab.w
#country weights Table B8
synthms2.tables$tab.pred



########SYNTHETIC CASE STUDIES: ISRAEL FP US diff ######
rm(list=ls())
library(haven)
setwd("~/Dropbox/Women Leaders and Conflict/Empirical Analysis/Rdata and scripts")
isrfp1 <- read_dta("isr50-73fem_fp936r.dta")
library(Synth)

colnames(isrfp1)[4]=c("female")
#isrfp1[isrfp1$cown == 666, c("logsmilex_f2", "lagdv")]

### 1. Prepping data
macro <- aggregate(female ~ cown, data = isrfp1, sum)
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
    tr.years <- isrfp1$year[(isrfp1$female==1 & isrfp1$cown==treated[i])]
    tr.start[i] <- min(tr.years)
    tr.end[i] <- max(tr.years)
}
isrfp1 <- merge(isrfp1, macro, by="cown")

isrfp1$fpdiff<-numeric(960)
isrfp1$fpdiff=isrfp1$idealpoint2_v-isrfp1$idealpoint2_v[isrfp1$cown==2]

### 2. Estimating Synthetic Counterfactual for Israel 
isrfp.dp <- dataprep(
    foo=isrfp1, 
    predictors = c("polity2_20", "logrgdp_f2", "gdpgrowth_pw", "lntenure_a", "cinc_c" ),
    predictors.op = "mean", 
    special.predictors = list(list("fpdiff", 1968, c("mean"))),
    dependent = "fpdiff", 
    unit.variable = "cown", 
    time.variable = "year",  
    treatment.identifier = 666, 
    controls.identifier = donor, 
    time.predictors.prior = c(1964:1968), 
    time.optimize.ssr = c(1964:1968), 
    unit.names.variable = "country", 
    time.plot = c(1964:tr.end[i]) 
)
isrfp.synth <- synth(isrfp.dp)
synthfp.tables <- synth.tab( dataprep.res = isrfp.dp, synth.res = isrfp.synth) 
isrfp.pred <- isrfp.dp$Y0plot %*% isrfp.synth$solution.w ####Y1s in previous code
isrfp.alpha<- isrfp.dp$Y1plot - isrfp.pred


### 3. Plot results

plot(1964:1973, isrfp.dp$Y1plot, ylim = c(-4, 2),  type = "n", lwd = 2, xlab = "", ylab = "")
box(); grid()
abline(v = 1969, lwd = 1)  
lines(1964:1973, isrfp.dp$Y1plot, lwd = 2)
lines(1964:1973, isrfp.pred, lwd = 2, lty = 2)
title(xlab = "Year", cex.lab = 1.5, line = 2)
title(ylab = "Ideal Point Estimates", cex.lab = 3, line = 2)
legend(list(x = 1964,y = -1.5), bty = "n", c("Israel", "Synthetic Israel"), lwd = 2, col = "black", lty = c(1, 2, -1), cex = 1.5, y.intersp=0.25)

#text(1963, -1, "First year of Meir's tenure", cex = 1.6)
text(1971.5, 1.5, bquote(paste("MSPE" == .(round(isrfp.synth$loss.v, 5)))), cex = 1.5)
#arrows(1967, -1, 1969, -1, length = .1)
rug(isrfp1$fpdiff, side = 2, lwd=1, ticksize = .035)


### 4. Placebo cases ##
### Run placebo for each donor unit
for (j in 1:length(donor))
{
    isrfp_pl.dp<- 
        dataprep( 
            foo = isrfp1, 
            predictors = c("polity2_20", "logrgdp_f2", "gdpgrowth_pw", "lntenure_a", "cinc_c" ),
            special.predictors = list(list("fpdiff", 1968, c("mean"))),
            predictors.op = "mean", 
            dependent = "fpdiff", 
            unit.variable = "cown", 
            time.variable = "year",  
            treatment.identifier = donor[j], 
            controls.identifier = donor[-j], 
            time.predictors.prior = c(1964:1968), 
            time.optimize.ssr = c(1964:1968), 
            unit.names.variable = "country", 
            time.plot = c(1964:tr.end[i]) 
        ) 
    
    isrfp_synth_pl <- synth(isrfp_pl.dp)
    synthfp_pl.tables <- synth.tab( dataprep.res = isrfp_pl.dp, synth.res = isrfp_synth_pl, round.digit = 5) 
    assign(paste("W",treated[i],donor[j], sep="."), (synthfp_pl.tables$tab.w))  # country weights
    assign(paste("V",treated[i],donor[j], sep="."), (isrfp_synth_pl$solution.v)) # covariate weights
    assign(paste("Pred",treated[i],donor[j], sep="."), (synthfp_pl.tables$tab.pred)) #sum table for synth, original and sample mean
    
    assign(paste("Y1",treated[i],donor[j], sep="."), (isrfp_pl.dp$Y1))
    assign(paste("Y1s",treated[i],donor[j], sep="."), (isrfp_pl.dp$Y0plot%*%isrfp_synth_pl$solution.w))
    assign(paste("Alpha",treated[i],donor[j], sep="."), (isrfp_pl.dp$Y1-isrfp_pl.dp$Y0plot%*%isrfp_synth_pl$solution.w))
}


### 5. Post Estimation results ###

## pre-RMSPE comparison ##
isr.prermspe <- sqrt(mean(isrfp.alpha[1:5]^2, na.rm=TRUE))
isr.postrmspe <- sqrt(mean(isrfp.alpha[6:10]^2, na.rm=TRUE))
isr.prermspe.d <- numeric(length(donor))
isr.postrmspe.d <- numeric(length(donor))

for (j in 1:length(donor))
{
    temp <- get(paste("Alpha", treated[i], donor[j], sep="."))
    isr.prermspe.d[j] <- sqrt(mean(temp[1:5]^2,na.rm=TRUE ))
    isr.postrmspe.d[j] <- sqrt(mean(temp[6:10]^2,na.rm=TRUE ))
}

# 5.1 e(avg_pre_rmspe_p) A measure of fit. Concerning if significant.
## The proportion of placebos that have a pre- treatment RMSPE at least as large as the average of the treated units. 
rmspe.fit <- length(which(isr.prermspe.d>isr.prermspe))/length(donor)


## 5.2  e(pval_joint_post_t)
## The proportion of placebos that have a ratio of post-treatment RMSPE over pre-treatment RMSPE at least as large as the average ratio for the treated units.
pval.joint.post_t <- length(which(isr.postrmspe.d/isr.prermspe.d>isr.postrmspe/isr.prermspe))/length(donor)

alphas <-isrfp.alpha
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
    pvals_t[t] <-  length(which(abs(alphas[t+5,-1]/isr.prermspe.d)>abs(alphas[t+5,1]/isr.prermspe)))/length(donor)
}

plot(1:5, pvals_t, ylim=c(0,1), axes = FALSE,xlab = "", ylab = "",  pch=20, bg="grey", lwd=5, cex=2.5)
box(); grid()
axis(1, at = 1:5, labels = seq(1, 5, 1), cex.axis = 1.25)
axis(2, at = seq(0, 1, 0.1), cex.axis = 1, las = 2)
title(xlab = "Number of years after treatment", cex.lab = 1.5, line = 2.4)
title(ylab = "Lead specific P-Values", cex.lab = 1.9, line = 2.6)


### 6. Plotting placebo cases and Israel

plot(1964:1973, isrfp.alpha, type = "n", ylim = c(-4, 2), lwd = 2, xlab = "", ylab = "", axes = FALSE) 
box(); grid()
for(j in 1: 39){
    lines(1964:1973, alphas[, j], type = "l", col = "gray80", lwd = 2)
} 
lines(1964:1973, alphas[,1 ], lwd = 2) #Israel 
abline(v=1969, lwd = 1)
abline(h=0, lwd = 1)
axis(1, at = 1964:1973, labels = seq(1964, 1973, 1), cex.axis = 1.25)
axis(2, at = seq(-3, 3, 1), cex.axis = 1.25, las = 2)
title(xlab = "Year", cex.lab = 1.5, line = 2)
title(ylab = "Difference in Ideal-Point estimates", cex.lab = 1.5, line = 2)
text(1966, -3, "Treatment time", cex = 1.2)
arrows(1967, -3, 1969, -3, length = .1)
legend(list(x = 1969.5,y = -1.8), ncol = 1, bty = "n", c("Israel", "Placebo Countries"), lty=c(1,1), lwd = 2, col = c("black", "gray80"), cex = 1.25, y.intersp = 0.25)

#country weights Table B15
synthfp.tables$tab.w
#Balance Table B9
synthfp.tables$tab.pred


########SYNTHETIC CASE STUDIES: UK Military Spending ####

rm(list=ls())
setwd("~/Dropbox/Women Leaders and Conflict/Empirical Analysis/Rdata and scripts")

library(haven)
ukms1 <- read_dta("UK_ms.dta")

library(Synth)
library(foreign)
library(tictoc)
library(plm)
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
ukms.synth <- synth(ukms.dp, "ALL")

synthms.tables <- synth.tab( dataprep.res = ukms.dp, synth.res = ukms.synth) 
ukms.pred <- ukms.dp$Y0plot %*% ukms.synth$solution.w ####Y1s in Previous code
ukms.alpha<- ukms.dp$Y1plot - ukms.pred

plot(1960:1990, ukms.dp$Y1plot, ylim = c(3, 13),  type = "n", lwd = 2, xlab = "", ylab = "")
box(); grid()
abline(v = 1979, lwd = 1)  
lines(1960:1990, ukms.dp$Y1plot, lwd = 2)
lines(1960:1990, ukms.pred, lwd = 2, lty = 2)
title(xlab = "Year", cex.lab = 1.5, line = 2)
title(ylab = "Military Spending", cex.lab = 2.5, line = 2)
legend(list(x = 1960,y = 14.8), bty = "n", c("UK", "Synthetic UK"), lwd = 2, col = "black", lty = c(1, 2, -1), cex = 1.3, y.intersp=0.3)
text(1969.5, 6, "First year of Thatcher's tenure", cex = 1.5)
arrows(1976, 6, 1979, 6, length = .2)
text(1985, 12.5, bquote(paste("MSPE" == .(round(ukms.synth$loss.v, 5)))), cex = 1.5)
rug(ukms1$logsmilex_f2, side = 2, lwd=1, ticksize = .025)

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
    #bad.fit <- numeric(33)
    #bad.fit[j] <- as.numeric(ifelse(ukms_synth_pl$loss.v > ukms.synth$loss.v*5, 1, 0))
     
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

## 5.3.2 Plotting p values
plot(1:12, pvals_t, ylim=c(0,1), axes = FALSE,xlab = "", ylab = "",  pch=20, bg="grey", lwd=5, cex=2.5)
box(); grid()
axis(1, at = 1:10, labels = seq(1, 10, 1), cex.axis = 1.25)
axis(2, at = seq(0, 1, 0.1), cex.axis = 1, las = 2)
title(xlab = "Number of years after treatment", cex.lab = 1.5, line = 2.4)
title(ylab = "Lead specific P-Values", cex.lab = 1.9, line = 2.6)


### 6. Plotting placebo cases and UK
plot(1960:1990, ukms.alpha, type = "n", ylim = c(-3, 3), lwd = 2, xlab = "", ylab = "", axes = FALSE) 
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
text(1972.3, -2.48, "Treatment time", cex = 1.4)
arrows(1976.5, -2.5, 1979, -2.5, length = .1)
legend(list(x = 1958,y = 3), ncol = 1, bty = "n", c("UK", "Placebo Countries"), lty=c(1,1), lwd = 2, col = c("black", "gray80"), cex = 1.25, y.intersp=0.6)

#country weights Table B10
synthms.tables$tab.w
#Balance Table B4
synthms.tables$tab.pred


########SYNTHETIC CASE STUDIES: UK Foreign Policy US DIFF ####
rm(list=ls())
setwd("~/Dropbox/Women Leaders and Conflict/Empirical Analysis/Rdata and scripts")

library(haven)
ukfp1 <- read_dta("UK_fp.dta")

library(Synth)
library(foreign)
library(tictoc)

colnames(ukfp1)[4]=c("female")

ukfp1$fpdiff<-numeric(1054)
ukfp1$fpdiff=ukfp1$idealpoint2_v-ukfp1$idealpoint2_v[ukfp1$cown==2]

### 1. Prepping data
macro <- aggregate(female ~ cown, data = ukfp1, sum)
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
    tr.years <- ukfp1$year[(ukfp1$female==1 & ukfp1$cown==treated[i])]
    tr.start[i] <- min(tr.years)
    tr.end[i] <- max(tr.years)
}
ukfp1 <- merge(ukfp1, macro, by="cown")


### 2. Estimating Synthetic Counterfactual for UK
### 2.1 Setting data
ukfp.dp <- dataprep(
    foo=ukfp1, 
    predictors = c("polity2_20", "logrgdp_f2", "gdpgrowth_pw", "lntenure", "cinc_c" ),
    special.predictors = list(list("fpdiff", 1978, c("mean"))),
    predictors.op = "mean", 
    dependent = "fpdiff",
    unit.variable = "cown", 
    time.variable = "year",  
    treatment.identifier = 200, 
    controls.identifier = donor[1:33], 
    time.predictors.prior = c(1974:1978), 
    time.optimize.ssr = c(1974:1978), 
    unit.names.variable = "country", 
    time.plot = c(1974:tr.end[i]) 
)
### 2.2 Running synth
ukfp.synth <- synth(ukfp.dp, "ALL")

synthfp.tables <- synth.tab( dataprep.res = ukfp.dp, synth.res = ukfp.synth) 
ukfp.pred <- ukfp.dp$Y0plot %*% ukfp.synth$solution.w ####Y1s in Previous code
ukfp.alpha<- ukfp.dp$Y1plot - ukfp.pred

plot(1974:1990, ukfp.dp$Y1plot, ylim = c(-4, 2),  type = "n", lwd = 2, xlab = "", ylab = "")
box(); grid()
abline(v = 1979, lwd = 1)  
lines(1974:1990, ukfp.dp$Y1plot, lwd = 2)
lines(1974:1990, ukfp.pred, lwd = 2, lty = 2)
title(xlab = "Year", cex.lab = 1.5, line = 2)
title(ylab = "Ideal-Point Estimates", cex.lab = 2.5, line = 2)
legend(list(x = 1973,y = 3), bty = "n", c("UK", "Synthetic UK"), lwd = 2, col = "black", lty = c(1, 2, -1), cex = 1.3, y.intersp=0.3)
text(1976, -2, "First year of ", cex = 1.1)
text(1976.3, -2.4, "Thatcher's tenure", cex = 1.1)
arrows(1978, -2.3, 1979, -2.3, length = .2)
text(1987, 1.5, bquote(paste("MSPE" == .(round(ukfp.synth$loss.v, 5)))), cex = 1.5)
rug(ukfp1$fpdiff, side = 2, lwd=1, ticksize = .025)

donor2<-c(20,  42,  70,  90,  94, 100, 130, 135, 140, 150, 155, 205, 210, 211, 212, 220, 230, 235, 305, 325,
          350, 375, 380, 390, 560, 600, 630, 740, 800, 820, 900, 920)
### 4. Placebo cases ##
### Run placebo for each donor unit
for (j in 1:length(donor))
{
    ukfp_pl.dp<- 
        dataprep( 
            foo = ukfp1, 
            predictors = c("polity2_20", "logrgdp_f2", "gdpgrowth_pw", "lntenure", "cinc_c" ),
            special.predictors = list(list("fpdiff", 1978, c("mean"))),
            predictors.op = "mean", 
            dependent = "fpdiff",
            unit.variable = "cown", 
            time.variable = "year",  
            treatment.identifier = donor[j], 
            controls.identifier = donor[-j], 
            time.predictors.prior = c(1974:1978), 
            time.optimize.ssr = c(1974:1978), 
            unit.names.variable = "country", 
            time.plot = c(1974:tr.end[i]) 
        ) 
    
    ukfp_synth_pl <- synth(ukfp_pl.dp)
    synthfp_pl.tables <- synth.tab( dataprep.res = ukfp_pl.dp, synth.res = ukfp_synth_pl, round.digit = 5) 
    assign(paste("W",treated[i],donor[j], sep="."), (synthfp_pl.tables$tab.w))
    assign(paste("V",treated[i],donor[j], sep="."), (ukfp_synth_pl$solution.v))
    assign(paste("Pred",treated[i],donor[j], sep="."), (synthfp_pl.tables$tab.pred))
    
    assign(paste("Y1",treated[i],donor[j], sep="."), (ukfp_pl.dp$Y1))
    assign(paste("Y1s",treated[i],donor[j], sep="."), (ukfp_pl.dp$Y0plot%*%ukfp_synth_pl$solution.w))
    assign(paste("Alpha",treated[i],donor[j], sep="."), (ukfp_pl.dp$Y1-ukfp_pl.dp$Y0plot%*%ukfp_synth_pl$solution.w))
}

### 5. Post Estimation results ###
## pre-RMSPE comparison ##
uk.prermspe <- sqrt(mean(ukfp.alpha[1:5]^2, na.rm=TRUE))
uk.postrmspe <- sqrt(mean(ukfp.alpha[6:17]^2, na.rm=TRUE))
uk.prermspe.d <- numeric(length(donor))
uk.postrmspe.d <- numeric(length(donor))

for (j in 1:length(donor))
{
    temp <- get(paste("Alpha", treated[i], donor[j], sep="."))
    uk.prermspe.d[j] <- sqrt(mean(temp[1:5]^2,na.rm=TRUE ))
    uk.postrmspe.d[j] <- sqrt(mean(temp[6:17]^2,na.rm=TRUE ))
}

# 5.1 e(avg_pre_rmspe_p) A measure of fit. Concerning if significant.
## The proportion of placebos that have a pre- treatment RMSPE at least as large as the average of the treated units. 
rmspe.fit <- length(which(uk.prermspe.d>uk.prermspe))/length(donor)


## 5.2  e(pval_joint_post_t)
## The proportion of placebos that have a ratio of post-treatment RMSPE over pre-treatment RMSPE at least as large as the average ratio for the treated units.
pval.joint.post_t <- length(which(uk.postrmspe.d/uk.prermspe.d>uk.postrmspe/uk.prermspe))/length(donor)

alphas <-ukfp.alpha
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
    pvals_t[t] <-  length(which(abs(alphas[t+5,-1]/uk.prermspe.d)>abs(alphas[t+5,1]/uk.prermspe)))/length(donor)
}

## 5.3.2 Plotting p values
plot(1:12, pvals_t, ylim=c(0,1), axes = FALSE,xlab = "", ylab = "",  pch=20, bg="grey", lwd=5, cex=2.5)
box(); grid()
axis(1, at = 1:10, labels = seq(1, 10, 1), cex.axis = 1.25)
axis(2, at = seq(0, 1, 0.1), cex.axis = 1, las = 2)
title(xlab = "Number of years after treatment", cex.lab = 1.5, line = 2.4)
title(ylab = "Lead specific P-Values", cex.lab = 1.9, line = 2.6)


### 6. Plotting placebo cases and UK
plot(1974:1990, ukfp.alpha, type = "n", ylim = c(-3, 3), lwd = 2, xlab = "", ylab = "", axes = FALSE) 
box(); grid()
for(j in 1: 33){
    lines(1974:1990, alphas[, j], type = "l", col = "gray80", lwd = 2)
} 

lines(1974:1990, alphas[,1 ], lwd = 2) #UK
abline(v=1979, lwd = 1)
abline(h=0, lwd = 1)
axis(1, at = 1974:1990, labels = seq(1974, 1990, 1), cex.axis = 1.25)
axis(2, at = seq(-3, 3, 1), cex.axis = 1.25, las = 2)
title(xlab = "Year", cex.lab = 1.5, line = 2)
title(ylab = "Difference in Ideal-Point estimates", cex.lab = 1.5, line = 2)
text(1975, -2.5, "Treatment time", cex = 1.2)
arrows(1978, -2.5, 1979, -2.5, length = .1)
legend(list(x = 1973,y = 4), ncol = 1, bty = "n", c("UK", "Placebo Countries"), lty=c(1,1), lwd = 2, col = c("black", "gray80"), cex = 1.25, y.intersp=0.3)

#country weights Table B11
synthfp.tables$tab.w
#Balance Table B5
synthfp.tables$tab.pred

########SYNTHETIC CASE STUDIES: India Military Spending ####
rm(list=ls())
setwd("~/Dropbox/Women Leaders and Conflict/Empirical Analysis/Rdata and scripts")

library(haven)
indms1 <- read_dta("India_ms56-84_725.dta")
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

###. 3 Plotting the results INDIA
plot(1956:1976, indms.dp$Y1plot, ylim = c(3, 13),  type = "n", lwd = 2, xlab = "", ylab = "")
box(); grid()
abline(v = 1966, lwd = 1)  
lines(1956:1976, indms.dp$Y1plot, lwd = 2)
lines(1956:1976, indms.pred, lwd = 2, lty = 2)
title(xlab = "Year", cex.lab = 1.5, line = 2)
title(ylab = "Military Spending", cex.lab = 2.5, line = 2)
legend(list(x = 1956,y = 15), bty = "n", c("India", "Synthetic India"), lwd = 2, col = "black", lty = c(1, 2, -1), cex = 1.3, y.intersp=0.3)
text(1960.5, 5, "First year of Gandhi's tenure", cex = 1.5)
arrows(1964.8, 5, 1966, 5, length = .2)
text(1973, 12.5, bquote(paste("MSPE" == .(round(indms.synth$loss.v, 5)))), cex = 1.5)
rug(indms1$logsmilex_f2, side = 2, lwd=1, ticksize = .025)

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

plot(1:11, pvals_t, ylim=c(0,1), axes = FALSE,xlab = "", ylab = "",  pch=20, bg="grey", lwd=5, cex=2.5)
box(); grid()
axis(1, at = 1:11, labels = seq(1, 11, 1), cex.axis = 1.25)
axis(2, at = seq(0, 1, 0.1), cex.axis = 1, las = 2)
title(xlab = "Number of years after treatment", cex.lab = 1.5, line = 2.4)
title(ylab = "Lead specific P-Values", cex.lab = 1.9, line = 2.6)


### 6. Plotting placebo cases and India
plot(1956:1976, indms.alpha, type = "n", ylim = c(-3, 3), lwd = 2, xlab = "", ylab = "", axes = FALSE) 
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
legend(list(x = 1954,y = 4), ncol = 1, bty = "n", c("India", "Placebo Countries"), lty=c(1,1), lwd = 2, col = c("black", "gray80"), cex = 1.25, y.intersp=0.3)

#country weights Table B12 
synthms.tables$tab.w
#Balance Table B6
synthms.tables$tab.pred



########SYNTHETIC CASE STUDIES: India FP US DIFF ####
rm(list=ls())
setwd("~/Dropbox/Women Leaders and Conflict/Empirical Analysis/Rdata and scripts")

library(haven)
indfp1 <- read_dta("ind_fp.usdiff.dta")
#indfp1 <- indfp1[indfp1$year<1977,]
#indfp1 <- indfp1[indfp1$year>1960,]

library(Synth)
library(foreign)

#colnames(indfp1)[4]=c("female")

### 1. Prepping data
macro <- aggregate(female ~ cown, data = indfp1, sum)
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
    tr.years <- indfp1$year[(indfp1$female==1 & indfp1$cown==treated[i])]
    tr.start[i] <- min(tr.years)
    tr.end[i] <- max(tr.years)
}
indfp1 <- merge(indfp1, macro, by="cown")

indfp1$fpdiff<-numeric(768)
indfp1$fpdiff=indfp1$idealpoint2_v-indfp1$idealpoint2_v[indfp1$cown==2]
predictors = c("polity2_20", "logrgdp_f2", "gdpgrowth_pw", "lntenure", "cinc_c" )

### 2. Estimating Synthetic Counterfactual for INDIA
indfp.dp <- dataprep(
    foo=indfp1, 
    predictors = predictors,
    special.predictors = list(list("fpdiff", 1965, c("mean"))),
    predictors.op = "mean", 
    dependent = "fpdiff", 
    unit.variable = "cown", 
    time.variable = "year",  
    treatment.identifier = 750, 
    controls.identifier = donor, 
    time.predictors.prior = c(1961:1965), 
    time.optimize.ssr = c(1961:1965), 
    unit.names.variable = "country", 
    time.plot = c(1961:tr.end[i]) 
)

indfp.synth <- synth(indfp.dp, "ALL")

synthfp.tables <- synth.tab( dataprep.res = indfp.dp, synth.res = indfp.synth) 
indfp.pred <- indfp.dp$Y0plot %*% indfp.synth$solution.w ####Y1s in Previous code
indfp.alpha<- indfp.dp$Y1plot - indfp.pred

###. 3 Plotting the results INDIA
plot(1961:1976, indfp.dp$Y1plot, ylim = c(-4, 2),  type = "n", lwd = 2, xlab = "", ylab = "")
box(); grid()
abline(v = 1966, lwd = 1)  
lines(1961:1976, indfp.dp$Y1plot, lwd = 2)
lines(1961:1976, indfp.pred, lwd = 2, lty = 2)
title(xlab = "Year", cex.lab = 1.5, line = 2)
title(ylab = "Ideal-Point Estimates", cex.lab = 2.5, line = 2)
rug(indfp1$fpdiff, side = 2, lwd=1, ticksize = .025)
legend(list(x = 1961,y = 2), bty = "n", c("India", "Synthetic India"), lwd = 2, col = "black", lty = c(1, 2, -1), cex = 1, y.intersp=0.3)
text(1963, -3, "First year of Gandhi's tenure", cex = 1)
arrows(1965, -3, 1966, -3, length = .2)
text(1974, 1.8, bquote(paste("MSPE" == .(round(indfp.synth$loss.v, 5)))), cex = 1.5)


### 4. Placebo cases ##
### Run placebo for each donor unit

for (j in 1:length(donor))
{
    indfp_pl.dp<- 
        dataprep( 
            foo = indfp1, 
            predictors = predictors,
            special.predictors = list(list("fpdiff", 1965, c("mean"))),
            predictors.op = "mean", 
            dependent = "fpdiff", 
            unit.variable = "cown", 
            time.variable = "year",  
            treatment.identifier = donor[j], 
            controls.identifier = donor[-j], 
            time.predictors.prior = c(1961:1965), 
            time.optimize.ssr = c(1961:1965), 
            unit.names.variable = "country", 
            time.plot = c(1961:tr.end[i]) 
        ) 
    indfp_synth_pl <- synth(indfp_pl.dp, "ALL")
    synthfp_pl.tables <- synth.tab( dataprep.res = indfp_pl.dp, synth.res = indfp_synth_pl, round.digit = 5) 
    assign(paste("W",treated[i],donor[j], sep="."), (synthfp_pl.tables$tab.w))
    assign(paste("V",treated[i],donor[j], sep="."), (indfp_synth_pl$solution.v))
    assign(paste("Pred",treated[i],donor[j], sep="."), (synthfp_pl.tables$tab.pred))
    
    assign(paste("Y1",treated[i],donor[j], sep="."), (indfp_pl.dp$Y1))
    assign(paste("Y1s",treated[i],donor[j], sep="."), (indfp_pl.dp$Y0plot%*%indfp_synth_pl$solution.w))
    assign(paste("Alpha",treated[i],donor[j], sep="."), (indfp_pl.dp$Y1-indfp_pl.dp$Y0plot%*%indfp_synth_pl$solution.w))
}

### 5. Post Estimation results ###
## pre-RMSPE comparison ##
ind.prermspe <- sqrt(mean(indfp.alpha[1:5]^2, na.rm=TRUE))
ind.postrmspe <- sqrt(mean(indfp.alpha[6:16]^2, na.rm=TRUE))
ind.prermspe.d <- numeric(length(donor))
ind.postrmspe.d <- numeric(length(donor))

for (j in 1:47)
{
    temp <- get(paste("Alpha", treated[i], donor[j], sep="."))
    ind.prermspe.d[j] <- sqrt(mean(temp[1:5]^2,na.rm=TRUE ))
    ind.postrmspe.d[j] <- sqrt(mean(temp[6:16]^2,na.rm=TRUE ))
}

# 5.1 e(avg_pre_rmspe_p) A measure of fit. Concerning if significant.
## The proportion of placebos that have a pre- treatment RMSPE at least as large as the average of the treated units. 
rmspe.fit <- length(which(ind.prermspe.d>ind.prermspe))/length(donor)


## 5.2  e(pval_joint_post_t)
## The proportion of placebos that have a ratio of post-treatment RMSPE over pre-treatment RMSPE at least as large as the average ratio for the treated units.
pval.joint.post_t <- length(which(ind.postrmspe.d/ind.prermspe.d>ind.postrmspe/ind.prermspe))/length(donor)

alphas <-indfp.alpha
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
    pvals_t[t] <-  length(which(abs(alphas[t+5,-1]/ind.prermspe.d)>abs(alphas[t+5,1]/ind.prermspe)))/length(donor)
}

plot(1:11, pvals_t, ylim=c(0,1), axes = FALSE,xlab = "", ylab = "",  pch=20, bg="grey", lwd=5, cex=2.5)
box(); grid()
axis(1, at = 1:11, labels = seq(1, 11, 1), cex.axis = 1.25)
axis(2, at = seq(0, 1, 0.1), cex.axis = 1, las = 2)
title(xlab = "Number of years after treatment", cex.lab = 1.5, line = 2.4)
title(ylab = "Lead specific P-Values", cex.lab = 1.9, line = 2.6)



### 6. Plotting placebo cases and India
plot(1961:1976, indfp.alpha, type = "n", ylim = c(-3, 3), lwd = 2, xlab = "", ylab = "", axes = FALSE) 
box(); grid()

for(j in 1: 47)
{lines(1961:1976, alphas[, j], type = "l", col = "gray80", lwd = 2)}
lines(1961:1976, alphas[,1 ], lwd = 2) #INDIA
abline(v=1966, lwd = 1)
abline(h=0, lwd = 1)
axis(1, at = 1961:1976, labels = seq(1961, 1976, 1), cex.axis = 1.25)
axis(2, at = seq(-3, 3, 1), cex.axis = 1.25, las = 2)
title(xlab = "Year", cex.lab = 1.8, line = 2)
title(ylab = "Difference in Ideal-Point Estimates", cex.lab = 1.5, line = 2)
text(1963, -2.5, "Treatment time", cex = 1.6)
arrows(1965, -2.5, 1966, -2.5, length = .1)
legend(list(x = 1959.2,y = 3), ncol = 1, bty = "n", c("India", "Placebo Countries"), lty=c(1,1), lwd = 2, col = c("black", "gray80"), cex = 1.35, y.intersp=0.4)

#country weights Table B13 
synthfp.tables$tab.w
#Balance Table B7
synthfp.tables$tab.pred




