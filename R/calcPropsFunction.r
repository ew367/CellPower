
#' Calculates the proportion of probes (CpGs) with power > 0.8
#'
#' @param allCellRes - a list of matrices with one matrix per cell type returned from running either calcDiff() or calcSamples()
#'
#' @return - a dataframe containing the proportion of probes with power > 0.8 where rows are either mean differences or number of samples tested, and columns are each cell type
#' @export
#'

# Example usage:

# allSamples <- calcSamples(allSDs, dataType = "SDs")
# allProps <-calcProps(allSamples)


calcProps <- function(allCellRes){

  prop_fun <- function(row){
    prop <- sum(row > 0.8)/length(row)
  }

  n <- allCellRes[[1]][,1]

  for(i in 1:length(names(allCellRes))){

    cellRes <- allCellRes[[i]]

    props <- apply(cellRes[,2:ncol(cellRes)], 1, prop_fun)
    n <- cbind(n, props)

  }

  colnames(n)[2:ncol(n)] <- names(allCellRes)
  return(n)

}

