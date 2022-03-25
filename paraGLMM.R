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

##

an2 <- lmer(para ~ hostsp + densNat + fYear + fWeek + site + (1 | fTree),
            family=binomial, data=paraData) 
an3 <- lmer(para ~ hostsp + densNat + fYear + fWeek + site + (1 | fTree),
            family=binomial, data=paraData)
anova(an2, an3)
summary(an2)
drop1(an3, test="Chi")
Anova(an3, type="III")

# MCMC

anM1 <- MCMCglmm(para ~ hostsp + densNat, random=~fTree, family="categorical", data=paraData, verbose=F)
# AIC 3705.6, all models with interactions increase AIC

cSite <- cbind(h=c(0, 0, 1), m=c(0, 1, 0)) # b = -1.186E+3, z = -1.419, p = 0.155761
hSite <- cbind(c=c(1, 0, 0), m=c(0, 0, 1)) # b = -4.885E+2, z = -0.600, p = 0.548611
mSite <- cbind(c=c(1, 0, 0), h=c(0, 1, 0)) # b = 6.125E+2, z = 1.333, p = 0.182484
siteCon <- cbind(c(-2,1,1), c(0,-1, 1))
  
contrasts(paraData$site) <- siteCon
an3 <- lmer(paraT ~ hostsp + densNat*site + fYear + fWeek + site + (1 | fTree),
            family=binomial, data=paraData)
an3a <- lmer(paraT ~ hostsp + densNat + fYear + fWeek + site + (1|as.factor(rep("A", length(paraT)))),
             family=binomial, data=paraData)
an3b <- lmer(paraT ~ densNat*site + (1|fTree), family=binomial, data=paraData)
drop1(an3a, test="Chi")
Anova(an3a, type="III")
anova(an3b)
summary(an3a)
# AIC 2623, interaction significant, but density insignificant for each site individually

an4 <- lmer(paraH ~ hostsp + densNat + fYear*site + fWeek + site + (1 | fTree),
          family=binomial, data=paraData)
summary(an4)
drop1(an4, test="Chi")
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

#################### FOCAL SPECIES ANALYSES #############################

# load data

focData <- read.csv("paraFocalGLMM.csv", header=T)

# create factors from numeric variables

focData$fYear <- as.factor(focData$year)
focData$fWeek <- as.factor(focData$week)
focData$fTree <- as.factor(focData$tree)
focData$densFoc <- focData$DensFoc

### analysis

# simple logistic regression

an1 <- glm(para ~ host + cat + densIFoc + fYear + fWeek + site, 
           family=binomial, data=focData) # AIC = 2282.9
drop1(an1, test="Chi")

an2 <- glm(paraT ~ host + cat + densIFoc + fYear + fWeek + site, 
           family=binomial, data=focData) # AIC = 1711.2
drop1(an2, test="Chi")

an3 <- glm(paraH ~ host + cat + densIFoc + fYear + fWeek + site, 
           family=binomial, data=focData) # AIC = 1292.4
drop1(an3, test="Chi")

# binomial GLMM

# all parasitoids

an4 <- lmer(para ~ host + cat + densIFoc + fYear + fWeek + site + (1 | fTree), 
           family=binomial, data=focData) # AIC = 2284.1
an5 <- lmer(para ~ host + cat + densIFoc*site + fYear + fWeek + site + (1 | fTree), 
            family=binomial, data=focData) # AIC 2282.7
an6 <- lmer(para ~ host + cat + densIFoc*site + fYear*site + fWeek + (1 | fTree), 
            family=binomial, data=focData) # AIC 2272.0
summary(an6)
drop1(an6, test="Chi")

summary(glmmPQL(para ~ host + cat + densIFoc*site + fYear*site + fWeek, 
                random = (~1 | fTree), family=binomial, data=focData))
# 0.7916537^2/(0.7916537^2+.8816011^2) # = ICC = 0.4463989

# flies only

an7 <- lmer(paraT ~ host + cat + densIFoc + fYear + fWeek + site + (1 | fTree), 
            family=binomial, data=focData) # AIC = 1713.076
an8 <- lmer(paraT ~ host + cat + densIFoc*fYear + fYear + fWeek + site + (1 | fTree), 
            family=binomial, data=focData) # AIC = 1709.8
an9 <- lmer(paraT ~ host + cat + densIFoc*fYear + densIFoc*site + fYear + fWeek + site + (1 | fTree), 
            family=binomial, data=focData) # AIC = 1709.393
an10 <- lmer(paraT ~ host + cat + densIFoc*fYear + densIFoc*site + fWeek + (1 | fTree), 
            family=binomial, data=focData) # AIC = 1709.393
AIC(an10)
drop1(an10, test='Chi')

summary(glmmPQL(paraT ~ host + cat + densIFoc*fYear + densIFoc*site + fWeek, 
                random = (~1 | fTree), family=binomial, data=focData))
# 1.730617^2/(1.730617^2+.6480582^2) # = ICC = 0.8770197

# contrasts for year and site

cSite <- cbind(h=c(0,1,0), m=c(0, 0, 1)) # b = -3.000E+3, z = -1.166, p = 0.243800
hSite <- cbind(c=c(1, 0, 0), m=c(0, 0, 1)) # b = -6.566E+3, z = -1.500, p = 0.133644
mSite <- cbind(c=c(1, 0, 0), h=c(0, 1, 0)) # b = 6.610E+2, z = 0.66, p = 0.508088

year4 <- cbind(o5=c(0,1,0,0,0), o6=c(0,0,1,0,0), o7=c(0,0,0,1,0),o8=c(0,0,0,0,1)) # b = -7.944E+3, z = -1.805, p = 0.07089
year5 <- cbind(o4=c(1,0,0,0,0), o6=c(0,0,1,0,0), o7=c(0,0,0,1,0),o8=c(0,0,0,0,1)) # b = -1.590E+3, z = - 0.488, p = 0.626403
year6 <- cbind(o4=c(1,0,0,0,0), o5=c(0,1,0,0,0), o7=c(0,0,0,1,0),o8=c(0,0,0,0,1)) # b = 6.148E+3, z = 2.506, p = 0.012212
year7 <- cbind(o4=c(1,0,0,0,0), o5=c(0,1,0,0,0), o6=c(0,0,1,0,0),o8=c(0,0,0,0,1)) # b = -3.543E+2, z = -0.290, p = 0.772145
year8 <- cbind(o4=c(1,0,0,0,0), o5=c(0,1,0,0,0), o6=c(0,0,1,0,0),o7=c(0,0,0,1,0)) # b = 5.555E+3, z = 1.020, p = 0.307764

# analyses

contrasts(focData$site) <- mSite
an10a <- lmer(paraT ~ host + cat + fYear + densIFoc*site + fWeek + (1 | fTree), 
             family=binomial, data=focData) # AIC = 1709.393
summary(an10a)

contrasts(focData$fYear) <- year8

an10b <- lmer(paraT ~ host + cat + densIFoc*fYear + site + fWeek + (1 | fTree), 
             family=binomial, data=focData) # AIC = 1709.393
summary(an10b)
# waps only

an11 <- lmer(paraH ~ host + cat + densIFoc + fYear + fWeek + site + (1 | fTree), 
            family=binomial, data=focData) # AIC = 1284.674
an12 <- lmer(paraH ~ host + cat*site + densIFoc + fYear + fWeek + site + (1 | fTree), 
             family=binomial, data=focData) # AIC = 1281.877
an13 <- lmer(paraH ~ host + cat*densIFoc + cat*site + fYear + fWeek  + (1 | fTree), 
             family=binomial, data=focData) # AIC = 1281.538
summary(an13)
AIC(an13)
drop1(an13, test="Chi")

summary(glmmPQL(paraH ~ host + cat*densIFoc + cat*site + fYear + fWeek, 
                random = (~1 | fTree), family=binomial, data=focData))
 # 2.544015^2/(2.544015^2+.4804167^2) # = ICC = 0.9655667

# create contrasts for individual caterpillar comparisons

adCat <- rbind(ad=c(0, 0, 0, 0, 0, 0), es=c(1, 0, 0, 0, 0, 0), hi=c(0, 1, 0, 0, 0, 0),
                mc=c(0, 0, 1, 0, 0, 0), md=c(0, 0, 0, 1, 0, 0), or=c(0, 0, 0, 0, 1, 0), 
                pt=c(0, 0, 0, 0, 0, 1)) # b = -2.540E+3, z = -0.116, p = 0.9074
esCat <- rbind(ad=c(1, 0, 0, 0, 0, 0), es=c(0, 0, 0, 0, 0, 0), hi=c(0, 1, 0, 0, 0, 0), 
                mc=c(0, 0, 1, 0, 0, 0), md=c(0, 0, 0, 1, 0, 0), or=c(0, 0, 0, 0, 1, 0), 
                pt=c(0, 0, 0, 0, 0, 1)) # b = -1.807E+4, z = -0.458, p = 0.6573
hiCat <- rbind(ad=c(1, 0, 0, 0, 0, 0), es=c(0, 1, 0, 0, 0, 0), hi=c(0, 0, 0, 0, 0, 0), 
                mc=c(0, 0, 1, 0, 0, 0), md=c(0, 0, 0, 1, 0, 0), or=c(0, 0, 0, 0, 1, 0), 
                pt=c(0, 0, 0, 0, 0, 1)) # b = 7.769E+3, z = 0.659, p = 0.5096
mcCat <- rbind(ad=c(1, 0, 0, 0, 0, 0), es=c(0, 1, 0, 0, 0, 0), hi=c(0, 0, 1, 0, 0, 0), 
                mc=c(0, 0, 0, 0, 0, 0), md=c(0, 0, 0, 1, 0, 0), or=c(0, 0, 0, 0, 1, 0), 
                pt=c(0, 0, 0, 0, 0, 1)) # b = 8.386E+2, z = 0.946, p = 0.34414
mdCat <- rbind(ad=c(1, 0, 0, 0, 0, 0), es=c(0, 1, 0, 0, 0, 0), hi=c(0, 0, 1, 0, 0, 0), 
                mc=c(0, 0, 0, 1, 0, 0), md=c(0, 0, 0, 0, 0, 0), or=c(0, 0, 0, 0, 1, 0), 
                pt=c(0, 0, 0, 0, 0, 1)) # b = 4.499E+4, z = 1.658, p = 0.09734
orCat <- rbind(ad=c(1, 0, 0, 0, 0, 0), es=c(0, 1, 0, 0, 0, 0), hi=c(0, 0, 1, 0, 0, 0), 
                mc=c(0, 0, 0, 1, 0, 0), md=c(0, 0, 0, 0, 1, 0), or=c(0, 0, 0, 0, 0, 0), 
                pt=c(0, 0, 0, 0, 0, 1)) # b = 9.523E+3, z = 1.683, p = 0.0923
ptCat <- rbind(ad=c(1, 0, 0, 0, 0, 0), es=c(0, 1, 0, 0, 0, 0), hi=c(0, 0, 1, 0, 0, 0), 
                mc=c(0, 0, 0, 1, 0, 0), md=c(0, 0, 0, 0, 1, 0), or=c(0, 0, 0, 0, 0, 1), 
                pt=c(0, 0, 0, 0, 0, 0)) # b = -7.958E+3, z = -1.196, p = 0.231536

contrasts(focData$cat) <- ptCat

an13 <- lmer(paraH ~ host + cat*densIFoc + cat*site + fYear + fWeek  + (1 | fTree), 
             family=binomial, data=focData) # AIC = 1281.53
summary(an13)