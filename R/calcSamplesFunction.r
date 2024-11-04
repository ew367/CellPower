#' Calculate the number of samples needed in each group to have power to detect the specified mean difference
#'
#' @param betasOrSDs - the name of an r matrix object containing either:
# 1) a normalised betas matrix of cell type specific DNA methylation data (sites as rows and samples as columns), or
# 2) a matrix of standard deviations derived where each column represents a different cell type or condition (CpGs as rows and columns as cell types)
#' @param pheno -  the name of an r dataframe object containing the phenotype data that matches the normalised betas matrix, with a column called 'Cell.type' specifying the cell type of each sample. MUST BE SUPPLIED IN ORDER TO USE A NORMALISED BETAS MATRIX.
#' @param meanDiff - Numeric value of the percentage difference in methylation you want to detect (default = 5)
#' @param dataType - Character string of the data type input provided to the betasOrSDs argument. One of "betaMatrix" or "SDs"
#' @param binSize - Numeric value for the number of bins to use. To improve package efficiency and speed, CpGs with similar standard deviations are binned together. Increasing this number may improve the accuracy of the output, but will increase the processing time (default = 500)
#'
#' @return - Returns a list of matrices with one matrix per cell type/condition. Each matrix contains the power for each site (columns) and different sample sizes (rows).
#' @export
#' @import pwr
#' @import doParallel
#' @import foreach
#' @importFrom dplyr mutate
#' @importFrom dplyr left_join
#'
#' @examples
#' allSamples <- calcSamples(allSDs, dataType = "SDs")
#' allSamples <- calcSamples(betas, phenoFile, meanDiff = 5, dataType = "betaMatrix")


calcSamples <- function(betasOrSDs, pheno=NULL, meanDiff=5, dataType, binSize = 500){

  if(dataType == "betaMatrix" & is.null(pheno)){
    stop("a phenotype file must be supplied alongside a betas matrix")
  }

  meanDiffDec <- as.numeric(meanDiff)/100 # calculate as percentage
  print(paste0("level of difference to be tested: ", meanDiff, "%"))

  # if input is matrix of normalised beta values, calculate standard deviations for each cell type

  if(dataType == "SDs"){

    allSDs <- betasOrSDs

  } else if (dataType == "betaMatrix"){

    allSDs <- c()

    for(i in unique(pheno$Cell.type)){

      # subset to cell type
      betasCell <- as.matrix(betasOrSDs[,colnames(betasOrSDs) %in% pheno[which(pheno$Cell.type == i), "Basename"]])
      print(paste0(ncol(betasCell), " ", i, " samples found."))

      # calculate SD
      betasSD <- apply(as.matrix(betasCell),1,sd)

      # return SDs
      betasSD <- as.matrix(betasSD, labels = T)
      colnames(betasSD) <- i
      allSDs <- cbind(allSDs, betasSD)
    }

  } else {
    stop("dataType must be one of betaMatrix or SDs")
  }


  allSDs <- allSDs[complete.cases(allSDs),] # remove NA values from matrix

  # calculate max number of samples needed
  maxSamples <- c()

  for(i in 1:ncol(allSDs)){

    cellSDs <- allSDs[,i]

    nSamples = 100
    pwr = 0.1
    while(sum(pwr > 0.8)/length(pwr) < 0.85){
      pwr <- pwr.t.test(n=nSamples, d = meanDiffDec/(cellSDs), sig.level = 9e-8,type="two.sample",alternative="two.sided")$power
      nSamples = nSamples+100
    }

    maxSamples <- c(maxSamples, nSamples)
  }

  nSamples <- max(maxSamples)


  # remove variables that are no longer needed
  rm(list=setdiff(ls(), c("allSDs", "nSamples", "pheno", "meanDiff", "meanDiffDec", "binSize")))


  #### Create output tables for each cell type ====

  #define function to calculate power for each sample size
  powerCpG <- function(betasSD){
    require(pwr)
    power1 <- c()
    Samples <- seq(nSamples/100, nSamples, nSamples/100)
    meandiffdec <- meanDiffDec

    # bin data
    betasSD <- as.data.frame(betasSD) %>% mutate(points_bin = ntile(betasSD, n=binSize))

    meanSDs <- c()
    for(j in unique(betasSD$points_bin)){
      meanSD <- mean(betasSD[which(betasSD$points_bin == j),1])
      meanSDs[j] <- meanSD
    }


    for(sample in Samples){

      power1 <- c(power1, pwr.t.test(n=sample, d = meandiffdec/(meanSDs), sig.level = 9e-8,type="two.sample",alternative="two.sided")$power)

      }

    power1 <- as.data.frame(matrix(power1, nrow=100, byrow = T))

    returnList <- list(betasSD, power1)
    return(returnList)
  }


  print("Setting up parallel processors...")

    no_cores <- detectCores() - 1
    cl <- makeForkCluster(no_cores)
    registerDoParallel(cl)
    clusterExport(cl, c("allSDs", "nSamples", "meanDiffDec", "binSize"), envir=environment())


    samples <- seq(nSamples/100, nSamples, nSamples/100)
    results <- foreach(i = colnames(allSDs)) %dopar% {powerCpG(allSDs[,i])}
    stopCluster(cl)
    names(results) <- colnames(allSDs)


    # match back to full cpg list

    allCellRes <- list()

    for(i in 1:length(names(results))){

    betaSDs <- as.data.frame(results[[i]][1])
    powerOut <- as.data.frame(t(as.data.frame(results[[i]][2])))
    powerOut$points_bin <- rep(1:binSize)
    outRes <- dplyr::left_join(betaSDs, powerOut, by="points_bin")
    rownames(outRes) <- rownames(betaSDs)
    outRes <- outRes[,3:ncol(outRes)]
    outRes <- cbind(samples, t(outRes))

    allCellRes[[i]] <- outRes

  }

    names(allCellRes) <- names(results)
    return(allCellRes)

}
