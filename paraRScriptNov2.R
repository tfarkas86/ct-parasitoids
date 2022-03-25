# Analysis on line 164

# load parasitism data

setwd("~/Dropbox/Projects/CT Parasitoids/RData")
classes <- c("factor", "factor", "character", "factor", "factor", 
             "character", "factor", "factor", "logical", "integer", 
             "numeric", "character", "factor", "factor", 
             "factor", "factor", "integer")
data <- read.table("allCatsRawParaData.csv", header=T, sep=",", colClasses=classes)
data$tree <- as.factor(paste(data$hostID, data$date, data$site, sep="_"))
data$weekHost <- as.factor(paste(data$year, data$week,data$host, sep="_"))
data$hostType <- ifelse(is.na(data$noLvs), "-", "Q")
data$paraT <- ifelse(data$paraType == "T", T, F)
data$paraH <- ifelse(data$paraType == "H", T, F)
data$para <- ifelse(data$para == "Y", T,F)

# create caterpillar groups in data data.frame

data$isQCat <- data$catType == "Q" & data$isCat
data$isNative <- data$isCat & (data$catSp != "Lymantria dispar")
data$isGen <- data$isCat & (data$breadth == "G")
data$isFocal <- (data$catSp == "Phigalia titea") | (data$catSp == "Ennomos subsignaria") | 
  (data$catSp == "Malacosoma disstria") | (data$catSp == "Melanolophia canadaria") |
  (data$catSp == "Himella intractata") | (data$catSp == "Achatia distincta") |
  (data$catSp == "Orthosia rubescens")

# create tree-level data frame

hostData <- data.frame(year=data$year, week=data$week, date=data$date, site=data$site, 
                       host=data$host, hostID=data$hostID,
                       noLvs=data$noLvs, avgLA=data$avgLA)
treeData <- unique(hostData)
treeData$tree <- as.factor(paste(treeData$hostID, treeData$date, treeData$site, sep="_"))
treeData$totLA <- ifelse(is.na(treeData$noLvs), 0, treeData$noLvs * treeData$avgLA)
treeData$weekHost <- as.factor(paste(treeData$year, treeData$week, treeData$host, sep="_"))

# add caterpillars to treeData data frame

treeData$noCats <- sapply(treeData$tree, function(x) sum((x == data$tree) & data$isCat))
treeData$noQCats <- sapply(treeData$tree, function(x) sum((x == data$tree) & data$isQCat))
treeData$noQNats <-  sapply(treeData$tree, function(x) sum((x == data$tree) & data$isQCat & data$isNative))
treeData$noQGens <- sapply(treeData$tree, function(x) sum((x == data$tree) & data$isQCat & data$isGen))
treeData$noQGenNats <- sapply(treeData$tree, function(x) sum((x == data$tree) & data$isQCat & data$isGen & data$isNative))
treeData$noQFocs <- sapply(treeData$tree, function(x) sum((x == data$tree) & data$isQCat & data$isFocal))
treeData$noRows <- sapply(treeData$tree, function(x) sum(x == data$tree))

treeLA <- vector()
for (i in 1:length(treeData[,1])) {
  treeLA <- c(treeLA, rep(treeData$totLA[i], treeData$noRows[i]))
}
data$treeLA <- treeLA

treeNoCats <- vector()
for (i in 1:length(treeData[,1])) {
  treeNoCats <- c(treeNoCats, rep(treeData$noCats[i], treeData$noRows[i]))
}
data$treeNoCats <- treeNoCats

treeNoQCats <- vector()
for (i in 1:length(treeData[,1])) {
  treeNoQCats <- c(treeNoQCats, rep(treeData$noQCats[i], treeData$noRows[i]))
}
data$treeNoQCats <- treeNoQCats

treeNoQNats <- vector()
for (i in 1:length(treeData[,1])) {
  treeNoQNats <- c(treeNoQNats, rep(treeData$noQNats[i], treeData$noRows[i]))
}
data$treeNoQNats <- treeNoQNats

treeNoQGens <- vector()
for (i in 1:length(treeData[,1])) {
  treeNoQGens <- c(treeNoQGens, rep(treeData$noQGens[i], treeData$noRows[i]))
}
data$treeNoQGens <- treeNoQGens

treeNoQGenNats <- vector()
for (i in 1:length(treeData[,1])) {
  treeNoQGenNats <- c(treeNoQGenNats, rep(treeData$noQGenNats[i], treeData$noRows[i]))
}
data$treeNoQGenNats <- treeNoQGenNats

treeNoQFocs <- vector()
for (i in 1:length(treeData[,1])) {
  treeNoQFocs <- c(treeNoQFocs, rep(treeData$noQFocs[i], treeData$noRows[i]))
}
data$treeNoQFocs <- treeNoQFocs

data$treeDens <- data$treeNoQCats / data$treeLA
data$treeDensGen <- data$treeNoQGens / data$treeLA
data$treeDensNat <- data$treeNoQNats / data$treeLA 
data$treeDensGenNats <- data$treeNoQGenNats / data$treeLA
data$treeDensFocs <- data$treeNoQFocs / data$treeLA

# create week-level data frame

weekHostData <- data.frame(year=data$year, week=data$week, 
                           host=data$host)
weekData <- unique(weekHostData)
weekData$weekHost <- as.factor(paste(weekData$year, weekData$week, weekData$host, sep="_"))

# add caterpillars and leaf area to weekData data frame

weekData$totLA <-sapply(weekData$weekHost, function(x) sum((x == treeData$weekHost) * treeData$totLA))
weekData$noRows <- sapply(weekData$weekHost, function(x) sum(x == data$weekHost))
weekData$noCats <- sapply(weekData$weekHost, function(x) sum((x == data$weekHost) & data$isCat))
weekData$noQCats <- sapply(weekData$weekHost, function(x) sum((x == data$weekHost) & data$isQCat))
weekData$noQNats <-  sapply(weekData$weekHost, function(x) sum((x == data$weekHost) & data$isQCat & data$isNative))
weekData$noQGens <- sapply(weekData$weekHost, function(x) sum((x == data$weekHost) & data$isQCat & data$isGen))
weekData$noQGenNats <- sapply(weekData$weekHost, function(x) sum((x == data$weekHost) & data$isQCat & data$isGen & data$isNative))
weekData$noQFocs <- sapply(weekData$weekHost, function(x) sum((x == data$weekHost) & data$isQCat & data$isFocal))

# add week-level density information to data data.frame

weekLA <- vector()
for (i in 1:length(weekData[,1])) {
  weekLA <- c(weekLA, rep(weekData$totLA[i], weekData$noRows[i]))
}
data$weekLA <- weekLA

weekNoCats <- vector()
for (i in 1:length(weekData[,1])) {
  weekNoCats <- c(weekNoCats, rep(weekData$noCats[i], weekData$noRows[i]))
}
data$weekNoCats <- weekNoCats

weekNoQCats <- vector()
for (i in 1:length(weekData[,1])) {
  weekNoQCats <- c(weekNoQCats, rep(weekData$noQCats[i], weekData$noRows[i]))
}
data$weekNoQCats <- weekNoQCats

weekNoQNats <- vector()
for (i in 1:length(weekData[,1])) {
  weekNoQNats <- c(weekNoQNats, rep(weekData$noQNats[i], weekData$noRows[i]))
}
data$weekNoQNats <- weekNoQNats

weekNoQGens <- vector()
for (i in 1:length(weekData[,1])) {
  weekNoQGens <- c(weekNoQGens, rep(weekData$noQGens[i], weekData$noRows[i]))
}
data$weekNoQGens <- weekNoQGens

weekNoQGenNats <- vector()
for (i in 1:length(weekData[,1])) {
  weekNoQGenNats <- c(weekNoQGenNats, rep(weekData$noQGenNats[i], weekData$noRows[i]))
}
data$weekNoQGenNats <- weekNoQGenNats

weekNoQFocs <- vector()
for (i in 1:length(weekData[,1])) {
  weekNoQFocs <- c(weekNoQFocs, rep(weekData$noQFocs[i], weekData$noRows[i]))
}
data$weekNoQFocs <- weekNoQFocs

data$weekDens <- data$weekNoQCats / data$weekLA
data$weekDensGen <- data$weekNoQGens / data$weekLA
data$weekDensNat <- data$weekNoQNats / data$weekLA 
data$weekDensGenNats <- data$weekNoQGenNats / data$weekLA
data$weekDensFocs <- data$weekNoQFocs / data$weekLA

data2 <- data[-901:-942,]
data3 <- data2[c(-2377,-2378),]
data4 <- data3[-4396,]
data5 <- data4[-5953,]
data <- data5

####### GLMM ANALYSIS #########

library(lme4)
library(car)

# all generalists (excl. L. dispar), tree-level-Native-generalist density, model selection by forward minimization of AIC

an1 <- lmer(para ~ year + week + site + host + treeDensGenNats + (1|tree), subset=(isNative & isGen & (hostType=="Q")), data=data, family=binomial)
an2 <- lmer(paraT ~ year + site + week + host + treeDensGenNats + (1|tree),  subset=(isNative & isGen & (hostType=="Q")), data=data, family=binomial)
an3 <- lmer(paraH ~ year + week + site + host + treeDensGenNats + (1|tree), subset=(isNative & isGen & (hostType=="Q")), data=data, family=binomial)

# all generalists (excl. L. dispar), tree-level-all-generalist density, model selection by forward min of AIC

an4 <- lmer(para ~ year + week + site + host + treeDensGen + (1|tree), subset=(isNative & isGen & (hostType=="Q")), data=data, family=binomial)
an5 <- lmer(paraT ~ year + week + site + host + treeDensGen + (1|tree), subset=(isNative & isGen & (hostType=="Q")), data=data, family=binomial)
an6 <- lmer(paraH ~ year + week + site + host + treeDensGen  + (1|tree), subset=(isNative & isGen & (hostType=="Q")), data=data, family=binomial)

# all generalists (excl. L. dispar), week-level-Native-generalist density, model selection by forward min of AIC

an7 <- lmer(para ~ year + week + site + host + weekDensGenNats + year*site + (1|tree), subset=(isNative & isGen), data=data4, family=binomial)
an7a <- lmer(para ~ year + week + site + host + weekDensGenNats + (1|tree), subset=(isNative & isGen), data=data4, family=binomial)
an8 <- lmer(paraT ~ year + week + site + host + weekDensGenNats + (1|tree), subset=(isNative & isGen), data=data, family=binomial)
an9 <- lmer(paraH ~ year + week + site + host + weekDensGenNats + year*site + (1|tree), subset=(isNative & isGen), data=data, family=binomial)
an9a <- lmer(paraH ~ year + week + site + host + weekDensGenNats + (1|tree), subset=(isNative & isGen), data=data, family=binomial)

# all generalists (excl. L. dispar), week-level-all-generalist density, model selection by forward min of AIC

an10 <- lmer(para ~ year + week + site + host + weekDensGen + year*site + (1|tree), subset=(isNative & isGen), data=data, family=binomial)
an10a <- lmer(para ~ year + week + site + host + weekDensGen + (1|tree), subset=(isNative & isGen), data=data, family=binomial) 
an11 <- lmer(paraT ~ year + week + site + host + weekDensGen + (1|tree), subset=(isNative & isGen), data=data, family=binomial)
an12 <- lmer(paraH ~ year + week + site + host + weekDensGen + year*site + (1|tree), subset=(isNative & isGen), data=data, family=binomial)
an12a <- lmer(paraH ~ year + week + site + host + weekDensGen + (1|tree), subset=(isNative & isGen), data=data, family=binomial)

# focal caterpillars, tree-level-Native-generalist density, model selection by forward min of AIC

an13 <- lmer(para ~ year + week + site + host + catSp + treeDensGenNats + (1|tree), subset=(isFocal & hostType=="Q"), data=data, family=binomial)
an14 <- lmer(paraT ~ year + week + site + host + catSp + treeDensGenNats + (1|tree), subset=(isFocal & hostType=="Q"), data=data, family=binomial)
an15 <- lmer(paraH ~ year + week + site + host + catSp + treeDensGenNats + year*site + (1|tree), subset=(isFocal & hostType=="Q"), data=data, family=binomial)
an15a <- lmer(paraH ~ year + week + site + host + catSp + treeDensGenNats + (1|tree), subset=(isFocal & hostType=="Q"), data=data, family=binomial)

# focal caterpillars, tree-level-all-generalist density, model selection by forward min of AIC

an16 <- lmer(para ~ year + week + site + host + catSp + treeDensGen + (1|tree), subset=(isFocal & hostType=="Q"), data=data, family=binomial)
an17 <- lmer(paraT ~ year + week + site + host + catSp + treeDensGen + (1|tree), subset=(isFocal & hostType=="Q"), data=data, family=binomial)
an18 <- lmer(paraH ~ year + week + site + host + catSp + treeDensGen + year*site + (1|tree), subset=(isFocal & hostType=="Q"), data=data, family=binomial)
an18a <- lmer(paraH ~ year + week + site + host + catSp + treeDensGen + (1|tree), subset=(isFocal & hostType=="Q"), data=data, family=binomial)

# focal caterpillars, week-level-Native-generalist density

an19 <- lmer(para ~ year + week + site + host + catSp + weekDensGenNats + (1|tree), subset=(isFocal), data=data, family=binomial)
an20 <- lmer(paraT ~ year + week + site + host + catSp + weekDensGenNats + year*site + (1|tree), subset=(isFocal), data=data, family=binomial)
an20a <- lmer(paraT ~ year + week + site + host + catSp + weekDensGenNats + (1|tree), subset=(isFocal), data=data, family=binomial)
an21 <- lmer(paraH ~ year + week + site + host + catSp + weekDensGenNats + year*site + (1|tree), subset=(isFocal), data=data, family=binomial)
an21a <- lmer(paraH ~ year + week + site + host + catSp + weekDensGenNats  + (1|tree), subset=(isFocal), data=data, family=binomial)
an21b <- lmer(paraH ~ host + catSp + year + site + weekDensGenNats + (1|tree), subset=(isFocal), data=data, family=binomial)

# focal caterpillars, week-level-all-generalist density

an22 <- lmer(para ~ year + week + site + host + catSp + weekDensGen + year*site + catSp*weekDensGen + (1|tree) , subset=(isFocal), data=data, family=binomial)
an22a <- lmer(para ~ year + week + site + host + catSp + weekDensGen + (1|tree) , subset=(isFocal), data=data, family=binomial)
an23 <- lmer(paraT ~ year + week + site + host + catSp + weekDensGen  + (1|tree), subset=(isFocal), data=data, family=binomial)
an24 <- lmer(paraH ~ year + week + site + host + catSp + weekDensGen + year*site + (1|tree), subset=(isFocal), data=data, family=binomial)
an24a <- lmer(paraH ~ year + week + site + host + catSp + weekDensGen + (1|tree), subset=(isFocal), data=data, family=binomial)
##### Get AIC and p-values #####

d1.1 <- drop1(an1, test="Chi")
d1.2 <- drop1(an2, test="Chi")
d1.3 <- drop1(an3, test="Chi")
d1.4 <- drop1(an4, test="Chi")
d1.5 <- drop1(an5, test="Chi")
d1.6 <- drop1(an6, test="Chi")
d1.7 <- drop1(an7, test="Chi")
d1.7a <- drop1(an7a, test="Chi")
d1.7b <- rbind(d1.7a, d1.7[5,])
d1.8 <- drop1(an8, test="Chi")
d1.9 <- drop1(an9, test="Chi")
d1.9a <- drop1(an9a, test="Chi")
d1.9b <- rbind(d1.9a, d1.9[5,])
d1.10 <- drop1(an10, test="Chi")
d1.10a <- drop1(an10a, test="Chi")
d1.10b <- rbind(d1.10a, d1.10[5,])
d1.11 <- drop1(an11, test="Chi")
d1.12 <- drop1(an12, test="Chi")
d1.12a <- drop1(an12a, test="Chi")
d1.12b <- rbind(d1.12a, d1.12[5,])
d1.13 <- drop1(an13, test="Chi")
d1.14 <- drop1(an14, test="Chi")
d1.15 <- drop1(an15, test="Chi")
d1.15a <- drop1(an15a, test="Chi")
d1.15b <- rbind(d1.15a, d1.15[6,])
d1.16 <- drop1(an16, test="Chi")
d1.17 <- drop1(an17, test="Chi")
d1.18 <- drop1(an18, test="Chi")
d1.18a <- drop1(an18a, test="Chi")
d1.18b <- rbind(d1.18a, d1.18[6,])
d1.19 <- drop1(an19, test="Chi")
d1.20 <- drop1(an20, test="Chi")
d1.20a <- drop1(an20a, test="Chi")
d1.20b <- rbind(d1.20a, d1.20[6,])
d1.21 <- drop1(an21, test="Chi")
d1.21a <- drop1(an21a, test="Chi")
d1.21b <- rbind(d1.21a, d1.21[6,])
d1.22 <- drop1(an22, test="Chi") # 2 interactions
d1.22a <- drop1(an22a, test="Chi")
d1.22b <- rbind(d1.22a, d1.22[4:5,])
d1.23 <- drop1(an23, test="Chi")
d1.24 <- drop1(an24, test="Chi")
d1.24a <- drop1(an24a, test="Chi")
d1.24b <- rbind(d1.24a, d1.24[6,])

d1.1f <- cbind(d1.1, dAIC=d1.1$AIC[1:6]-d1.1$AIC[1])
d1.2f <- cbind(d1.2, dAIC=d1.2$AIC[1:6]-d1.2$AIC[1])
d1.3f <- cbind(d1.3, dAIC=d1.3$AIC[1:6]-d1.3$AIC[1])
d1.4f <- cbind(d1.4, dAIC=d1.4$AIC[1:6]-d1.4$AIC[1])
d1.5f <- cbind(d1.5, dAIC=d1.5$AIC[1:6]-d1.5$AIC[1])
d1.6f <- cbind(d1.6, dAIC=d1.6$AIC[1:6]-d1.6$AIC[1])
d1.7af <- cbind(d1.7a, dAIC=d1.7a$AIC[1:6]-d1.7a$AIC[1])
d1.7f <- rbind(d1.7af, "year:site"=c(d1.7[5,], dAIC=d1.7$AIC[5]-d1.7$AIC[1]))
d1.8f <- cbind(d1.8, dAIC=d1.8$AIC[1:6]-d1.8$AIC[1])
d1.9af <- cbind(d1.9a, dAIC=d1.9a$AIC[1:6]-d1.9a$AIC[1])
d1.9f <- rbind(d1.9af, "year:site"=c(d1.9[5,], dAIC=d1.9$AIC[5]-d1.9$AIC[1]))
d1.10af <- cbind(d1.10a, dAIC=d1.10a$AIC[1:6]-d1.10a$AIC[1])
d1.10f <- rbind(d1.10af, "year:site"=c(d1.10[5,], dAIC=d1.10$AIC[5]-d1.10$AIC[1]))
d1.11f <- cbind(d1.11, dAIC=d1.11$AIC[1:6]-d1.11$AIC[1])
d1.12af <- cbind(d1.12a, dAIC=d1.12a$AIC[1:6]-d1.12a$AIC[1])
d1.12f <- rbind(d1.12af, "year:site"=c(d1.12[5,], dAIC=d1.12$AIC[5]-d1.12$AIC[1]))
d1.13f <- cbind(d1.13, dAIC=d1.13$AIC[1:7]-d1.13$AIC[1])
d1.14f <- cbind(d1.14, dAIC=d1.14$AIC[1:7]-d1.14$AIC[1])
d1.15af <- cbind(d1.15a, dAIC=d1.15a$AIC[1:7]-d1.15a$AIC[1])
d1.15f <- rbind(d1.15af, "year:site"=c(d1.15[6,], dAIC=d1.15$AIC[6]-d1.15$AIC[1]))
d1.16f <- cbind(d1.16, dAIC=d1.16$AIC[1:7]-d1.16$AIC[1])
d1.17f <- cbind(d1.17, dAIC=d1.17$AIC[1:7]-d1.17$AIC[1])
d1.18af <- cbind(d1.18a, dAIC=d1.18a$AIC[1:7]-d1.18a$AIC[1])
d1.18f <- rbind(d1.18af, "year:site"=c(d1.18[6,], dAIC=d1.18$AIC[6]-d1.18$AIC[1]))
d1.19f <- cbind(d1.19, dAIC=d1.19$AIC[1:7]-d1.19$AIC[1])
d1.20af <- cbind(d1.20a, dAIC=d1.20a$AIC[1:7]-d1.20a$AIC[1])
d1.20f <- rbind(d1.20af, "year:site"=c(d1.20[6,], dAIC=d1.20$AIC[6]-d1.20$AIC[1]))
d1.21af <- cbind(d1.21a, dAIC=d1.21a$AIC[1:7]-d1.21a$AIC[1])
d1.21f <- rbind(d1.21af, "year:site"=c(d1.21[6,], dAIC=d1.21$AIC[6]-d1.21$AIC[1]))
d1.22af <- cbind(d1.22a, dAIC=d1.22a$AIC[1:7]-d1.22a$AIC[1])
d1.22f <- rbind(d1.22af, cbind(d1.22[4:5,], dAIC=d1.22$AIC[4:5]-d1.22$AIC[1]))
d1.23f <- cbind(d1.23, dAIC=d1.23$AIC[1:7]-d1.23$AIC[1])
d1.24af <- cbind(d1.24a, dAIC=d1.24a$AIC[1:7]-d1.24a$AIC[1])
d1.24f <- rbind(d1.24af, "year:site"=c(d1.24[6,], dAIC=d1.24$AIC[6]-d1.24$AIC[1]))

paraResults <- data.frame(round(rbind(d1.1f, d1.2f, d1.3f, d1.4f, d1.5f, d1.6f, d1.7f,d1.8f,
                     d1.9f, d1.10f, d1.11f, d1.12f, d1.13f, d1.14f, d1.15f, d1.16f, 
                     d1.17f, d1.18f, d1.19f, d1.20f, d1.21f, d1.22f, d1.23f, d1.24f),4))

setwd("~/Dropbox/Projects/CT Parasitoids/")
write.table(paraResults, "paraResultsNov9.csv", sep=",")

## Figures with Error Bars ##
load("~/Dropbox/Projects/CT Parasitoids/paraAnalysesWorkspace.RData")
library(lme4)
library(effects)

## generalist parasitism

ef1 <- effect("host", an7)
ef2 <- effect("host", an8)
ef3 <- effect("host", an9)

odds1 <- exp(ef1$fit)
odds2 <- exp(ef2$fit)
odds3 <- exp(ef3$fit)
prob1a <- odds1 / (odds1 + 1)
prob1b <- prob1a[-c(7:8)]
prob1 <- c(prob1a[7:8], prob1b)
prob2a <- odds2 / (odds2 + 1)
prob2b <- prob2a[-c(7:8)]
prob2 <- c(prob2a[7:8], prob2b)
prob3a <- odds3 / (odds3 + 1)
prob3b <- prob3a[-c(7:8)]
prob3 <- c(prob3a[7:8], prob3b)

uprOdds1 <- exp(ef1$upper)
uprOdds2 <- exp(ef2$upper)
uprOdds3 <- exp(ef3$upper)
uprProb1a <- uprOdds1 / (uprOdds1 + 1)
uprProb1b <- uprProb1a[-c(7:8)]
uprProb1 <- c(uprProb1a[7:8], uprProb1b)
uprProb2a <- uprOdds2 / (uprOdds2 + 1)
uprProb2b <- uprProb2a[-c(7:8)]
uprProb2 <- c(uprProb2a[7:8], uprProb2b)
uprProb3a <- uprOdds3 / (uprOdds3 + 1)
uprProb3b <- uprProb3a[-c(7:8)]
uprProb3 <- c(uprProb3a[7:8], uprProb3b)
uprProb3c <- uprProb3 + prob2

lwrOdds1 <- exp(ef1$lower)
lwrOdds2 <- exp(ef2$lower)
lwrOdds3 <- exp(ef3$lower)
lwrProb1a <- lwrOdds1 / (lwrOdds1 + 1)
lwrProb1b <- lwrProb1a[-c(7:8)]
lwrProb1 <- c(lwrProb1a[7:8], lwrProb1b)
lwrProb2a <- lwrOdds2 / (lwrOdds2 + 1)
lwrProb2b <- lwrProb2a[-c(7:8)]
lwrProb2 <- c(lwrProb2a[7:8], lwrProb2b)
lwrProb3a <- lwrOdds3 / (lwrOdds3 + 1)
lwrProb3b <- lwrProb3a[-c(7:8)]
lwrProb3 <- c(lwrProb3a[7:8], lwrProb3b)
lwrProb3c <- lwrProb3 + prob2

genMeans <- t(cbind(prob3, prob2, prob1))
genUpr <- t(cbind(uprProb3, uprProb2, uprProb1))
genLwr <- t(cbind(lwrProb3, lwrProb2, lwrProb1))

genMeans2 <- t(cbind(prob2, prob3))
genUpr2 <- t(cbind(uprProb2, uprProb3c))
genLwr2 <- t(cbind(lwrProb2, lwrProb3c))

names <- c("WH", "WO", "BC", "BE", "BI", "HI", "RM", "RO")

setwd("~/Dropbox/Projects/CT Parasitoids/Figures")
tiff("genParaBeside.tiff", width=5, height=5, unit="in", res=500)
mids <- barplot(genMeans, beside=T, names.arg=names, ylab = "adjusted probability of parasitism", 
                xlab = "host-plant species", ylim=c(0,.4), col=c("grey50", "grey70", "grey90"), las=1)
arrows(mids, genLwr, mids, genUpr, code=3, angle=90, length=0.03)
legend (x=0, y=.36, legend=c("wasps only", "flies only", "all parasitoids"), 
        fill=c("grey50", "grey70", "grey90"), border="grey100", bty="n")
text(x=1.7, y=.38, labels="a", cex=1.5)
dev.off()

#png("genParaStack.png")
#mids2 <- barplot(genMeans2, beside=F, names.arg=names, ylab = "adjusted probability of parasitism", ylim=c(0,.4))
#arrows(mids2-.1, lwrProb2, mids2-.1, uprProb2, code=3, angle=90, length=0.05)
#arrows(mids2+.1, lwrProb3c, mids2+.1, uprProb3c, code=3, angle=90, length=0.05)
#dev.off()

# focal parasitism

library(effects)

ef1 <- effect("host", an19)
ef2 <- effect("host", an20)


ef3 <- effect("host", an21b)

odds1 <- exp(ef1$fit)
odds2 <- exp(ef2$fit)
odds3 <- exp(ef3$fit)
prob1a <- odds1 / (odds1 + 1)
prob1b <- prob1a[-c(7:8)]
prob1 <- c(prob1a[7:8], prob1b)
prob2a <- odds2 / (odds2 + 1)
prob2b <- prob2a[-c(7:8)]
prob2 <- c(prob2a[7:8], prob2b)
prob3a <- odds3 / (odds3 + 1)
prob3b <- prob3a[-c(7:8)]
prob3 <- c(prob3a[7:8], prob3b)

uprOdds1 <- exp(ef1$upper)
uprOdds2 <- exp(ef2$upper)
uprOdds3 <- exp(ef3$upper)
uprProb1a <- uprOdds1 / (uprOdds1 + 1)
uprProb1b <- uprProb1a[-c(7:8)]
uprProb1 <- c(uprProb1a[7:8], uprProb1b)
uprProb2a <- uprOdds2 / (uprOdds2 + 1)
uprProb2b <- uprProb2a[-c(7:8)]
uprProb2 <- c(uprProb2a[7:8], uprProb2b)
uprProb3a <- uprOdds3 / (uprOdds3 + 1)
uprProb3b <- uprProb3a[-c(7:8)]
uprProb3 <- c(uprProb3a[7:8], uprProb3b)
uprProb3c <- uprProb3 + prob2

lwrOdds1 <- exp(ef1$lower)
lwrOdds2 <- exp(ef2$lower)
lwrOdds3 <- exp(ef3$lower)
lwrProb1a <- lwrOdds1 / (lwrOdds1 + 1)
lwrProb1b <- lwrProb1a[-c(7:8)]
lwrProb1 <- c(lwrProb1a[7:8], lwrProb1b)
lwrProb2a <- lwrOdds2 / (lwrOdds2 + 1)
lwrProb2b <- lwrProb2a[-c(7:8)]
lwrProb2 <- c(lwrProb2a[7:8], lwrProb2b)
lwrProb3a <- lwrOdds3 / (lwrOdds3 + 1)
lwrProb3b <- lwrProb3a[-c(7:8)]
lwrProb3 <- c(lwrProb3a[7:8], lwrProb3b)
lwrProb3c <- lwrProb3 + prob2

focMeans2 <- cbind(prob3, prob2, prob1)
focUpr2 <- cbind(uprProb3, uprProb2, uprProb1)
focLwr2<- cbind(lwrProb3, lwrProb2, lwrProb1)
datFrame <- data.frame(focMeans2, focUpr2, focLwr2)
datFrameSort <- datFrame[with(datFrame, order(prob1)),]

focMeans <- t(datFrameSort[,1:3])
focUpr <- t(datFrameSort[,4:6])
focLwr <- t(datFrameSort[,7:9])
focNames <- c("WH", "BI", "WO", "BC", "BE", "HI", "RO", "RM")
setwd("~/Dropbox/Projects/CT Parasitoids/Figures")
tiff("focParaBeside2.tiff", width=5, height=5, unit="in", res=500)
mids <- barplot(focMeans, beside=T, names.arg=focNames, ylab = NULL, 
                xlab = "host-plant species", ylim=c(0,.40), col=c("grey50", "grey70", "grey90"), 
                las=1)
arrows(mids, focLwr, mids, focUpr, code=3, angle=90, length=0.03)
#legend (x=0, y=.31, legend=c("wasps only", "flies only", "all parasitoids"), 
 #       fill=c("grey50", "grey70", "grey90"), border="grey100", bty="n")
text(x=1.7, y=.38, labels="b", cex=1.5)
dev.off()

# caterpillar species parasitism

ef1 <- effect("catSp", an19)
ef2 <- effect("catSp", an20)
ef3 <- effect("catSp", an21b)

odds1 <- exp(ef1$fit)
odds2 <- exp(ef2$fit)
odds3 <- exp(ef3$fit)
prob1a <- odds1 / (odds1 + 1)
prob1b <- prob1a[-c(7:8)]
prob1 <- c(prob1a[7:8], prob1b)
prob2a <- odds2 / (odds2 + 1)
prob2b <- prob2a[-c(7:8)]
prob2 <- c(prob2a[7:8], prob2b)
prob3a <- odds3 / (odds3 + 1)
prob3b <- prob3a[-c(7:8)]
prob3 <- c(prob3a[7:8], prob3b)

uprOdds1 <- exp(ef1$upper)
uprOdds2 <- exp(ef2$upper)
uprOdds3 <- exp(ef3$upper)
uprProb1a <- uprOdds1 / (uprOdds1 + 1)
uprProb1b <- uprProb1a[-c(7:8)]
uprProb1 <- c(uprProb1a[7:8], uprProb1b)
uprProb2a <- uprOdds2 / (uprOdds2 + 1)
uprProb2b <- uprProb2a[-c(7:8)]
uprProb2 <- c(uprProb2a[7:8], uprProb2b)
uprProb3a <- uprOdds3 / (uprOdds3 + 1)
uprProb3b <- uprProb3a[-c(7:8)]
uprProb3 <- c(uprProb3a[7:8], uprProb3b)
uprProb3c <- uprProb3 + prob2

lwrOdds1 <- exp(ef1$lower)
lwrOdds2 <- exp(ef2$lower)
lwrOdds3 <- exp(ef3$lower)
lwrProb1a <- lwrOdds1 / (lwrOdds1 + 1)
lwrProb1b <- lwrProb1a[-c(7:8)]
lwrProb1 <- c(lwrProb1a[7:8], lwrProb1b)
lwrProb2a <- lwrOdds2 / (lwrOdds2 + 1)
lwrProb2b <- lwrProb2a[-c(7:8)]
lwrProb2 <- c(lwrProb2a[7:8], lwrProb2b)
lwrProb3a <- lwrOdds3 / (lwrOdds3 + 1)
lwrProb3b <- lwrProb3a[-c(7:8)]
lwrProb3 <- c(lwrProb3a[7:8], lwrProb3b)
lwrProb3c <- lwrProb3 + prob2

catMeans2 <- cbind(prob3, prob2, prob1)
catUpr2 <- cbind(uprProb3, uprProb2, uprProb1)
catLwr2<- cbind(lwrProb3, lwrProb2, lwrProb1)
catFrame <- data.frame(focMeans2, focUpr2, focLwr2)
catFrame <- catFrame[-2,]
catFrameSort <- catFrame[with(catFrame, order(prob1)),]

catMeans <- t(catFrameSort[,1:3])
catUpr <- t(catFrameSort[,4:6])
catLwr <- t(catFrameSort[,7:9])
catNames <- rep("", 7)
setwd("~/Dropbox/Projects/CT Parasitoids/Figures")
tiff("catParaBeside.tiff", width=5, height=5, unit="in", res=500)

mids <- barplot(catMeans, beside=T, names.arg=catNames, ylab = "adjusted probability of parasitism", 
                ylim=c(0,.6), col=c("grey50", "grey70", "grey90"), 
                las=1)
arrows(mids, catLwr, mids, catUpr, code=3, angle=90, length=0.03)
legend (x=0, y=.6, legend=c("wasps only", "flies only", "all parasitoids"), 
       fill=c("grey50", "grey70", "grey90"), border="grey100", bty="n")
#legend (x=0, y=.31, legend=c("wasps only", "flies only", "all parasitoids"), 
#       fill=c("grey50", "grey70", "grey90"), border="grey100", bty="n")
dev.off()
