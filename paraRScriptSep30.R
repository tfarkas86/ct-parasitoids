# load parasitism data

setwd("~/Dropbox/Projects/CT Parasitoids/RData")
             "character", "factor", "factor", "logical", "integer", 
             "numeric", "character", "factor", "factor", 
             "factor", "factor", "integer")
data <- read.table("allCatsRawParaData.csv", header=T, sep=",", colClasses=classes)
data$tree <- as.factor(paste(data$hostID, data$date, data$site, sep="_"))
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

# create branch-level data frame

hostData <- data.frame(year=data$year, week=data$week, date=data$date, site=data$site, 
                       host=data$host, hostID=data$hostID,
                       noLvs=data$noLvs, avgLA=data$avgLA)
treeData <- unique(hostData)
treeData$tree <- as.factor(paste(treeData$hostID, treeData$date, treeData$site, sep="_"))
treeData$totLA <- treeData$noLvs * treeData$avgLA

# add caterpillars to treeData data frame

treeData$noCats <- sapply(treeData$tree, function(x) sum((x == data$tree) & data$isCat))
treeData$noQCats <- sapply(treeData$tree, function(x) sum((x == data$tree) & data$isQCat))
treeData$noQNats <-  sapply(treeData$tree, function(x) sum((x == data$tree) & data$isQCat & data$isNative))
treeData$noQGens <- sapply(treeData$tree, function(x) sum((x == data$tree) & data$isQCat & data$isGen))
treeData$noQGenNats <- sapply(treeData$tree, function(x) sum((x == data$tree) & data$isQCat & data$isGen & data$isNative))
treeData$noQFocs <- sapply(treeData$tree, function(x) sum((x == data$tree) & data$isQCat & data$isFocal))
treeData$noRows <- sapply(treeData$tree, function(x) sum(x == data$tree))

# add density information to data data.frame

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

data$dens <- data$treeNoQCats / data$treeLA
data$densGen <- data$treeNoQGens / data$treeLA
data$densNat <- data$treeNoQNats / data$treeLA 
data$densGenNats <- data$treeNoQGenNats / data$treeLA
data$densFocs <- data$treeNoQFocs / data$treeLA

####### ANALYSIS #########

# all generalist, native caterpillars on trees with density data, density = all caterpillars

an1 <- lmer(para ~ year + week + site + host + densGenNats + (1|tree), subset=(isNative & isGen & (hostType=="Q")), data=data, family=binomial)
Anova(an1, type="III")
an2 <- lmer(paraT ~ year + week + site + host + densGenNats + (1|tree), subset=(isNative & isGen & (hostType=="Q")), data=data, family=binomial)
Anova(an2, type="III")
an3 <- lmer(paraH ~ year + week + site + host + densGenNats + (1|tree), subset=(isNative & isGen & (hostType=="Q")), data=data, family=binomial)
Anova(an3, type="III")

# focal caterpillars only, density = all caterpillars

an1 <- lmer(para ~ year + week + site + host + catSp + densGenNats + year*site
            + (1|tree), subset=(isNative & isGen & (hostType=="Q") & isFocal), data=data, family=binomial)
Anova(an1, type="III")
an2 <- lmer(paraT ~ year + week + site + host + catSp + (1|tree), subset=(isFocal), data=data, family=binomial)
Anova(an2, type="III")
an3 <- lmer(paraH ~ year + week + site + host + densGenNats + (1|tree), subset=(isNative & isGen & (hostType=="Q") & isFocal), data=data, family=binomial)
Anova(an3, type="III")

# logistic regression

an1 <- glm(para ~ host + year + week + site + densGenNats, family=binomial, data=data, subset=(catSp == "Phigalia titea" & (hostType=="Q")))
Anova(an1, type=3)

q# ficin

noLvs <- as.list(data$noLvs)

substr(data$hostID[1], 3,3)
          