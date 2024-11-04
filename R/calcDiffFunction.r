#' Calculate the percentage difference in methylation detectable from the specified number of samples
#'
#' @param betasOrSDs - the name of an r matrix object containing either:
# 1) a normalised betas matrix of cell type specific DNA methylation data (sites as rows and samples as columns), or
# 2) a matrix of standard deviations derived where each column represents a different cell type or condition (CpGs as rows and columns as cell types)
#' @param pheno -  the name of an r dataframe object containing the phenotype data that matches the normalised betas matrix, with a column called 'Cell.type' specifying the cell type of each sample. MUST BE SUPPLIED IN ORDER TO USE A NORMALISED BETAS MATRIX.
#' @param nSamples - Numeric value of the number of samples in each group (default = 100))
#' @param dataType - Character string of the data type input into the function. One of "betaMatrix" or "SDs"
#' @param binSize - Numeric value for the number of bins to use. To improve package efficiency and speed, CpGs with similar standard deviations are binned together. Increasing this number may slightly improve the accuracy of the output, but will increase the processing time (default = 500)
#'
#' @return - Returns a list of matrices with one matrix per cell type/condition. Each matrix contains the power for each site (columns) and different mean differences (rows).
#' @export
#' @import pwr
#' @import doParallel
#' @import foreach
#' @importFrom dplyr mutate
#' @importFrom dplyr left_join
#' @importFrom dplyr ntile
#' @importFrom stats complete.cases
#' @importFrom parallel detectCores
#' @importFrom parallel makeForkCluster
#' @importFrom parallel clusterExport
#'
#' @examples
#' allSamples <- calcDiff(allSDs, dataType = "SDs")

# or for beta values (not supplied in this example): allSamples <- calcDiff(betas, phenoFile, nSamples = 100, dataType = "betaMatrix")



calcDiff <- function(betasOrSDs, pheno=NULL, nSamples=100, dataType, binSize = 500){

  if(dataType == "betaMatrix" & is.null(pheno)){
    stop("a phenotype file must be supplied to use a betas matrix")
  }


  print(paste0("Number of samples to be tested: ", nSamples))

  if(dataType == "SDs"){

    allSDs <- betasOrSDs

  } else if (dataType == "betaMatrix"){

    allSDs <- c() # vector to hold allSDs from loop

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


  rm(list=setdiff(ls(), c("allSDs", "pheno", "nSamples", "binSize"))) # remove variables that are no longer needed


  #define function to calculate power for each sample size
  powerCpG <- function(betasSD){
    require(pwr)
    power1 <- c()
    meanDiff <- seq(0.001, 0.1, 0.001)


    # bin data
    betasSD <- as.data.frame(betasSD) %>% dplyr::mutate(points_bin = dplyr::ntile(betasSD, n=binSize))

    meanSDs <- c()
    for(j in unique(betasSD$points_bin)){
      meanSD <- mean(betasSD[which(betasSD$points_bin == j),1])
      meanSDs[j] <- meanSD
    }

    for(diff in meanDiff){
      power1 <- c(power1, pwr::pwr.t.test(d=diff/meanSDs,n=nSamples,sig.level=9e-8,type="two.sample",alternative="two.sided")$power)
    }
    power1 <- as.data.frame(matrix(power1, nrow=length(meanDiff), byrow = T))

    returnList <- list(betasSD, power1)
    return(returnList)

  }


  print("Setting up parallel processors...")

    require(dplyr)
    require(foreach)
    require(doParallel)

chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")

if (nzchar(chk) && chk == "TRUE") {
    no_cores <- 2L
} else {
    no_cores <- parallel::detectCores() - 1
}

  cl <- parallel::makeForkCluster(no_cores)
  doParallel::registerDoParallel(cl)
  parallel::clusterExport(cl, c("allSDs", "nSamples", "binSize"), envir=environment())

  meanDiff <- seq(0.001, 0.1, 0.001)
  results <- foreach(i = colnames(allSDs)) %dopar% {powerCpG(allSDs[,i])}
  stopCluster(cl)
  names(results) <- colnames(allSDs)

  allCellRes <- list()

  for(i in 1:length(names(results))){

    betaSDs <- as.data.frame(results[[i]][1])
    powerOut <- as.data.frame(t(as.data.frame(results[[i]][2])))
    powerOut$points_bin <- seq(1:binSize)
    outRes <- dplyr::left_join(betaSDs, powerOut, by="points_bin")
    rownames(outRes) <- rownames(betaSDs)
    outRes <- outRes[,3:ncol(outRes)]
    outRes <- cbind(meanDiff, t(outRes))

    allCellRes[[i]] <- outRes

  }

  names(allCellRes) <- names(results)
  return(allCellRes)

}

