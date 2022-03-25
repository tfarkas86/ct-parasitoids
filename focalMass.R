setwd("~/Dropbox/Projects/CT Parasitoids/RData")
massData <- read.table("focalMass.txt", header=T)

massData$fYear <- factor(massData$year)
names(massData)

host <- c("BC", "BE", "BI", "HI", "RM", "RO", "WH", "WO")

library(nlme)

adData <- subset(massData, species=="Ad")
anAd <- lme(mass ~ host + sexCon, random = ~1 | fYear/site, data = adData)
summary(anAd)

intAd <- 0.0619
adMasses <- c(intAd, intAd-.0209, intAd-.0017, intAd-.0192, intAd-.0169, intAd-.0252, intAd-.0089, intAd-.0186)


esData <- subset(massData, species=="Es")

anEs <- lme(mass ~ host + sexCon, random = ~1 | fYear/site, data = esData)
summary(anEs)

intEs <- 0.02138
esMasses <- c(NA, intEs, intEs-.0033, intEs+.0116, intEs+.0031, intEs+.0022, intEs+.0039, intEs+.0003)

hiData <- subset(massData, species=="Hi")

anHi <- lme(mass ~ host + sexCon, random = ~1 | fYear/site, data = hiData)
summary(anHi)

intHi <- 0.0478
hiMasses <- c(intHi, intHi-.0191, intHi-.0144, intHi+.0104, intHi+.0119, intHi+.0142, intHi+.0085, intEs+.0227)

mcData <- subset(massData, species=="Mc")

anMc <- lme(mass ~ host + sexCon, random = ~1 | fYear/site, data = mcData)
summary(anMc)

intMc <- 0.0275
mcMasses <- c(intMc, NA, intMc-.0080, intMc+.0048, intMc+.0085, intMc+.0128, intMc+.0031, intMc+.0108)


mdData <- subset(massData, species=="Md")

anMd <- lme(mass ~ host + sexCon, random = ~1 | fYear/site, data = mdData)
summary(anMd)

intMd <- 0.0526
mdMasses <- c(intMd, intMd-.0011, intMd-.0167, intMd+.0026, intMd+.0030, intMd+.0075, intEs+.0013, intEs+.0114)

orData <- subset(massData, species=="Or")

anOr <- lme(mass ~ host + sexCon, random = ~1 | fYear/site, data = orData)
summary(anOr)

intOr <- 0.0623
orMasses <- c(intOr, intOr-.0249, intOr-.0110, intOr+.0114, intOr+.0147, intOr+.0192, intOr+.0131, intOr+.0233)

ptData <- subset(massData, species=="Pt")

anPt <- lme(mass ~ host + sexCon, random = ~1 | fYear/site, data = ptData)
summary(anPt)

intPt <- 0.0414
ptMasses <- c(intPt, intPt-.0256, intPt-.0132, intPt+.0009, intPt+.0075, intPt+.0080, intPt-.0012, intPt+.0047)

focMasses <- data.frame(host, adMasses, esMasses, hiMasses, mcMasses, mdMasses, orMasses, ptMasses)
massString <- c(adMasses, esMasses, hiMasses, mcMasses, mdMasses, orMasses, ptMasses)

paraData <- read.table("focalParaByHost.txt", header=T)
paraData$mass <- massString

adY <- cbind(AdParaData$pPara, AdParaData$pParaT, AdParaData$pParaH)
esY <- cbind(EsParaData$pPara, EsParaData$pParaT, EsParaData$pParaH)
hiY <- cbind(HiParaData$pPara, HiParaData$pParaT, HiParaData$pParaH)
mcY <- cbind(McParaData$pPara, McParaData$pParaT, McParaData$pParaH)
mdY <- cbind(MdParaData$pPara, MdParaData$pParaT, MdParaData$pParaH)
orY <- cbind(OrParaData$pPara, OrParaData$pParaT, OrParaData$pParaH)
ptY <- cbind(PtParaData$pPara, PtParaData$pParaT, PtParaData$pParaH)

AdParaData <- subset(paraData, catSp=="Ad")
EsParaData <- subset(paraData, catSp=="Es")
HiParaData <- subset(paraData, catSp=="Hi")
McParaData <- subset(paraData, catSp=="Mc")
MdParaData <- subset(paraData, catSp=="Md")
OrParaData <- subset(paraData, catSp=="Or")
PtParaData <- subset(paraData, catSp=="Pt")

adMassModA <- lm(pPara ~ mass, data=AdParaData)
esMassModA <- lm(pPara ~ mass, data=EsParaData)
hiMassModA <- lm(pPara ~ mass, data=HiParaData)
mcMassModA <- lm(pPara ~ mass, data=McParaData)
mdMassModA <- lm(pPara ~ mass, data=MdParaData)
orMassModA <- lm(pPara ~ mass, data=OrParaData)
ptMassModA <- lm(pPara ~ mass, data=PtParaData)

adMassModT <- lm(pParaT ~ mass, data=AdParaData)
esMassModT <- lm(pParaT ~ mass, data=EsParaData)
hiMassModT <- lm(pParaT ~ mass, data=HiParaData)
mcMassModT <- lm(pParaT ~ mass, data=McParaData)
mdMassModT <- lm(pParaT ~ mass, data=MdParaData)
orMassModT <- lm(pParaT ~ mass, data=OrParaData)
ptMassModT <- lm(pParaT ~ mass, data=PtParaData)

adMassModH <- lm(pParaH ~ mass, data=AdParaData)
esMassModH <- lm(pParaH ~ mass, data=EsParaData)
hiMassModH <- lm(pParaH ~ mass, data=HiParaData)
mcMassModH <- lm(pParaH ~ mass, data=McParaData)
mdMassModH <- lm(pParaH ~ mass, data=MdParaData)
orMassModH <- lm(pParaH ~ mass, data=OrParaData)
ptMassModH <- lm(pParaH ~ mass, data=PtParaData)

summary(adMassModA)
summary(esMassModA)
summary(hiMassModA)
summary(mcMassModA)
summary(mdMassModA)
summary(orMassModA)
summary(ptMassModA)

summary(adMassModT)
summary(esMassModT)
summary(hiMassModT)
summary(mcMassModT)
summary(mdMassModT)
summary(orMassModT)
summary(ptMassModT)

summary(adMassModH)
summary(esMassModH)
summary(hiMassModH)
summary(mcMassModH)
summary(mdMassModH)
summary(orMassModH)
summary(ptMassModH)

adMassMan <- manova(adY ~ mass, data=AdParaData)
esMassMan <- manova(esY ~ mass, data=EsParaData)
hiMassMan <- manova(hiY ~ mass, data=HiParaData)
mcMassMan <- manova(mcY ~ mass, data=McParaData)
mdMassMan <- manova(mdY ~ mass, data=MdParaData)
orMassMan <- manova(orY ~ mass, data=OrParaData)
ptMassMan <- manova(ptY ~ mass, data=PtParaData)

summary(adMassMan)
summary(esMassMan)
summary(hiMassMan)
summary(mcMassMan)
summary(mdMassMan)
summary(orMassMan)
summary(ptMassMan)


