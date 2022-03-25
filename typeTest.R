rm(list=ls())
sampData <- read.table("sampTest.txt", header=T)

an1 <- lm(sampData$catDens ~ sampData$type + sampData$timing + sampData$type*sampData$timing)
summary(an1)

