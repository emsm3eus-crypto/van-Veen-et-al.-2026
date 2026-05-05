
library(reshape)
library(ggplot2)
library(reshape)
library(car)
library(nlme)
require(multcomp)
library(rstudioapi)

dataFile <- selectDirectory(
  caption = "Select Directory",
  label = "Select",
  path = getActiveProject()
)
setwd(dataFile)

a <- list.files(pattern='exp', ignore.case=T) #this should be any common string among all of your sholl files
keyFiles <- list.files(pattern = 'Key', ignore.case=T)
c <- lapply(a, read.csv2)

for(i in 1:length(c)){
  tmp <- c[[i]]
  names(tmp)[1] <- 'Radius'
  c[[i]] <- tmp
}

h <- melt(c, id='Radius') 

h <- h[!h$variable == 'X',] 

f <- unique(h$variable)
ff <- data.frame(f, c(1:length(f)))
names(ff) <- c('name', 'num')

h$cellNum <- c()
k<-c()
for (i in 1:length(f)){
  tmp1 <- as.character(ff[ff$num == i, 1][1])
  tmp <- h[h$variable == tmp1, ]
  y <- rep(i, times=nrow(tmp))
  k <- c(k,y)
}
h$cellNum <- k

key <- read.csv(keyFiles[1])
print(key) 
print(a)

h <- merge(h, key, 'L1')
names(h) <- c('List', 'Radius', 'ImageName', 'Intersections', 'cellNum', 'Condition', 'FileName', 'rep')


head(h)
range(h$Radius)

h$Condition <- as.factor(h$Condition) 

cond_test_corAR <- nlme::lme(
  Intersections ~  1 + Condition,
  data = h,
  random = ~ 1 | cellNum / Radius, 
  control = lmeControl(opt = "optim")
)

cond_test_corAR
summary(cond_test_corAR)
print(anova(cond_test_corAR))
AOV <- anova(cond_test_corAR)
AOV[,"p-value"]
