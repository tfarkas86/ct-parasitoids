setwd("~/Dropbox/Projects/CT Parasitoids/RData")
paraData <- read.csv("paraDataGLMM.csv", header=T)

library(nlme)
library(MASS)
library(lme4)
library(glmmML)
library(pscl)

paraData$fYear <- as.factor(paraData$year)
paraData$fWeek <- as.factor(paraData$week)
paraData$fTree <- as.factor(paraData$tree)

an1 <- glmmPQL(paraH ~ hostsp + densNat + fYear + fWeek + site, random = (~1 | fTree), 
               family=binomial, data=paraData)
summary(an1)



#an2 <- lmer(para ~ hostsp + densNat + fYear + fWeek + site + (1 | fTree),
#            family=binomial, data=paraData) 
# AIC 3705.6, all models with interactions increase AIC

#####

#an3 <- lmer(paraT ~ hostsp + densNat + fYear + fWeek + site + (1 | fTree),
           # family=binomial, data=paraData)
# AIC 2624
cZero <- cbind(mVC=c(0, 0, 1), hVC=c(0, 1, 0))
hZero <- cbind(mVC=c(1, 0, 0), hVC=c(0, 0, 1))
mZero <- cbind(mVC=c(1, 0, 0), hVC=c(0, 1, 0))

contrasts(paraData$site) <- mZero
an3 <- lmer(paraT ~ hostsp + densNat*site + fYear + fWeek + site + (1 | fTree),
            family=binomial, data=paraData)
summary(an3)
# AIC 2623, interaction significant, but density insignificant for each site individually

######
an4 <- lmer(paraH ~ hostsp + densNat + fYear*site + fWeek + site + (1 | fTree),
          family=binomial, data=paraData)
summary(an4)
## AIC = 2181.9

###
drop1(an4, test="Chi")

an4 <- glm(paraH ~ hostsp + densNat + fYear + fWeek + site,
           family=binomial, data=paraData)
summary(an3)
drop1(an3, test="Chi")

## quantitative data only

qData <- subset(paraData, paraData$isQ. == 1)

an6 <- glmmPQL(paraT ~ hostsp + densNat + fWeek + fYear + site, random = (~1 | fTree), 
              family=binomial, data=qData)
summary(an6)

0.0004331161^2/(.0004331161^2+1.008807^2)

an5 <- lmer(paraH ~ hostsp + densNat + fWeek + fYear + site + (1 | fTree), 
            family=binomial, data=paraData)
summary(an5, test="Chi")
drop1(an5)

## randome effects models

an1 <- glmmPQL(para ~ 1, random = (~1 | fTree), family=binomial, data=paraData)
summary(an1)
an1 <- glmmPQL(paraH ~ hostsp, random = (~1 | fTree), family=binomial, data=paraData)
summary(an1)

1.48969^2/(1.48969^2+.7219514^2)