library(ggplot2)
library(sciplot)
library(magrittr) 
library(dplyr)    
library(glmmTMB)
library(plyr)
library(knitr)
library(PerformanceAnalytics)
library(psych)
library(rcompanion)
library(MASS)

dat = read.csv("blossoms.csv")
head(dat)
names(dat)
View(dat)
str(dat)

populations = as.factor(dat$pop)
patches = as.factor(dat$patch)

###############################################
# summarize data
###############################################

tapply(dat$UBL, dat$pop, mean, na.rm=T)
tapply(dat$LBL, dat$pop, mean, na.rm=T)

popstats2 = ddply(dat, .(pop), summarize,
                  UBLm = mean(UBL, na.rm=T),
                  UBLsd = sd(UBL, na.rm=T),
                  UBLse = se(UBL, na.rm=T),
                  LBLm = mean(LBL, na.rm=T),
                  LBLsd = sd(LBL, na.rm=T),
                  LBLse = se(LBL, na.rm=T))
popstats2[,-1] = round(popstats2[,-1], 2)
kable(popstats2)

#################################################################################################################

#################################################################################
# After exploring and summarizing the data, fit some linear models to estimate 
# the slopes of one trait on
# another. Interpret the results. Do the analysis on both arithmetic 
# and log scale. Choose traits that belong
# to the same vs. different functional groups, can you detect any patterns?
##################################################################################

# do regression on UBL and LBL for population S1
# notice slope slightly smaller than 1
# so usually lowr is bigger in a Plant
# then check through ancova how the various populations vary
oldpar = par(no.readonly = TRUE)

subdat = dat[dat$pop=="S1",]
y = subdat$UBL
x = subdat$LBL
reg = lm(y~x)
cf = summary(reg)$coef
predvals = cf[1,1] + cf[2,1]*x
eq = paste0("y = ", round(cf[2,1],1), "*x ", "+ ", round(cf[1,1],1))
par(mfrow=c(1,2))
newx = seq(12.09, 28.15, length.out =length(x))
predy = cf[1,1] + cf[2,1]*newx
plot(x, y, las=1,
     ylab="UBL data for population S1", xlab="LBL data for population S1", main=eq)
lines(newx, predy, col="red", lwd = 2 )
abline(a=0, b=1, col="blue", lwd = 2)
segments(x, y, x, predvals)
z = residuals(reg)
plotNormalHistogram(z, prob = FALSE, col="black", border="green",
                    main = "Normal Distribution Overlay on Residuals Histogram",
                    linecol="red", lwd=3 )
summary(reg)

fitdistr(residuals(reg), "normal")
coefs = summary(reg)$coef
(coefs[2,1]*(mean(x) + sd(x))) - (coefs[2,1]*mean(x))
sd(x)
##############################################################
#
# for +1 mm change in LBL, the UBL will change by + 0.9 mm (rounded) 
#
##############################################################




##############################################################################
##############################################################################
#
# Doing boxplots to visualize data of UBL and LBL
#
##############################################################################
##############################################################################
par(mfrow=c(1,1))
par(bg = "ivory")
boxplot(UBL~populations, data = dat, xlab="Populations", ylim = c(12, 28),
        ylab="UBL",boxwex=0.35, main = "Boxplot of UBL data for populations", lwd = 2, border= "slategrey", # colour of the box borders
        col = "slategray2", # colour of the inside of the boxes
        col.axis = 'grey20', # colour of the axis numbers 
        col.lab = 'grey20', # colour of the axis labels
        frame = F)
stripchart(UBL~populations,
           data = dat,
           method = "jitter",  main = "method = 'jitter', jitter = 0.2",
           pch = 16, # specify the type of point to use
           cex = 1,
           col = c(1,2,3,4,5,6,7,9),
           vertical = TRUE, 
           add = TRUE)
# at = c(0.75,1.75)
par(mfrow=c(1,1))
boxplot(LBL~populations, data = dat, xlab="Populations", ylim = c(12, 28),
        ylab="LBL",boxwex=0.35, main = "Boxplot of LBL data for populations", lwd = 2, border= "slategrey", # colour of the box borders
        col = "slategray2", # colour of the inside of the boxes
        col.axis = 'grey20', # colour of the axis numbers 
        col.lab = 'grey20', # colour of the axis labels
        frame = F)
stripchart(LBL~populations,
           data = dat,
           method = "jitter",  main = "method = 'jitter', jitter = 0.2",
           pch = 16, # specify the type of point to use
           cex = 1,
           col = c(1,2,3,4,5,6,7,9),
           vertical = TRUE, 
           add = TRUE)

par(oldpar)


##############################################################################
##############################################################################
#
# end of boxplots to visualize data
#
##############################################################################
##############################################################################


#################################################################################################################


#################################################################################
#
#
# now I try to do one way ANOVA to check differences between populations when it comes to UBL and LBL
#
#
##################################################################################



