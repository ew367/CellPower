

#### bin SDs and plot power

setwd("/lustre/projects/Research_Project-MRC190311/DNAm/powerCalcs/update2023/cmdArgTesting/")

library(pwr)
library(dplyr)
library(ggplot2)

load("SDs.rdat")

#dim(allSDs)
#[1] 866529      6

allSDs <- allSDs[complete.cases(allSDs),]
#dim(allSDs)
#[1] 865918      6

nSamples <- 100
meanDiff <- 0.01
nBins <- seq(20,800,20)

#start_time <- Sys.time()

allProps <- c()
for(cell in colnames(allSDs)){
  
  test <- as.data.frame(allSDs) %>% select(any_of(cell))
  
  props <-c()
for(i in nBins){
  
  test <- test %>% mutate(points_bin = ntile(test, n=i))
  
  meanSDs <- c()
  
  for(j in unique(test$points_bin)){
    meanSD <- mean(test[which(test$points_bin == j),1])
    meanSDs[j] <- meanSD 
  
  }
  
  # calc power for each meanSD
  pwrOut <- pwr.t.test(n=nSamples, d = meanDiff/(meanSDs), sig.level = 9e-8,type="two.sample",alternative="two.sided")$power
  props <- c(props, sum(pwrOut > 0.8)/length(pwrOut))
  
}
  allProps <- c(allProps, props)
  
}

#end_time <- Sys.time()
#print(end_time - start_time)




allProps <- as.data.frame(matrix(allProps, nrow=40, byrow = F))
colnames(allProps) <- colnames(allSDs)

check <- allProps
check$binSize <- nBins

plotdf <- reshape2::melt(check, id = "binSize")

p <- ggplot(plotdf, aes(x = binSize, y = value, colour = variable))+
  geom_line()+
  ylab("Proportion > 0.8")+
  ggtitle(paste0(nSamples, " samples, ", meanDiff*100, "% difference"))

p

plotFile <- paste0("plots/binSize", nSamples, "samples_", meanDiff, "meanDiff.pdf")
pdf("plots/binSizePlot.pdf")
print(p)
dev.off()

diffs <- c()
for(i in colnames(allProps)){
  d <- diff(allProps[,i])
  diffs <- c(diffs, d)
}

diffs <- as.data.frame(matrix(diffs, nrow=39, byrow = F))
colnames(diffs) <- colnames(allProps)
binSize <- seq(20, 780, 20)
diffs <- cbind(binSize, diffs)

save(diffs, file="binSizDiffs.rdat")

print(diffs[diffs$binSize > 460,])

# find rows where absolute value is < 0.01
diffsAbs <- abs(diffs)
diffsAbs[apply(diffsAbs[,2:ncol(diffsAbs)], 1, function(x){all(x<0.001)}),]




# test how long it take to run using 500 bins and 1000 bins

#load("SDs.rdat")
#source("Package/BrainPower/R/calcDiffFunctionV2.r")

#setwd("/lustre/projects/Research_Project-MRC190311/DNAm/powerCalcs/update2023/cmdArgTesting")

# load betas
#betas <- fread("data/celltypeNormbeta.csv", showProgress = T, data.table = F)
#betasMatrix <- betasMatrix[complete.cases(betasMatrix),]
#rownames(betasMatrix) <- betasMatrix$probeID

# load pheno
#QCmetrics <- read.csv("data/QCmetrics.csv", stringsAsFactors = F)

# run function
#start_time <- Sys.time()

#calcDiff(betas, QCmetrics, 200, betasMatrix, 500)

#end_time <- Sys.time()

#print(end_time - start_time)


