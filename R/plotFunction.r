#' plot power curve
#'
#' @param resFile - the dataframe with proportions of probes with certain level of power output from running calcProps() where each column is a 
#' cell type and each row is a different condition, e.g. different sample size or mean difference
#' @param calcType - character string of the type of power calcualtion run. Either "samples" (if calcSamples() was used) or "difference" (if 
#' calcDiff() was used).
#'
#' @return - ggplot object with power curve plot
#' @export
#'
#' @examples
#' myPlot <- plotPower(allProps, "samples")


plotPower <- function(resFile, calcType){

  resFile <- as.data.frame(resFile)

  if(calcType == "samples"){

    xlabs <- "N Samples (per group)"
    xbreaks <- seq(0, max(allProps[, "n"]), signif(max(allProps[, "n"])/9, 1))

  } else if (calcType == "difference"){

    resFile$n <- resFile$n*100
    xlabs <- "Mean Difference (%)"
    xbreaks <- seq(0, 10, 0.5)

  } else {
    stop("calcType must be one of 'samples' or 'difference'")
  }

  plotdf <- reshape2::melt(resFile, id.vars = "n")

  # create plot
  p <- ggplot(plotdf, aes(x = n, y = value, colour = variable))+
    geom_line()+
    labs(x = xlabs, y = "Proportion CpG's > 80% power", color = "cell types")+
    scale_x_continuous(breaks = xbreaks)+
    scale_y_continuous(breaks = seq(0, 1, by = 0.2))+

    # add in segment lines
    geom_segment(aes(x = 0, y = 0.8, xend = resFile$n[apply(resFile[,2:ncol(resFile)], 1, function(x){all(x>0.8)})][1], yend = 0.8), colour = "black", linetype = "dashed", linewidth = 0.1)


  for(i in 2:ncol(resFile)){
    p <- p + geom_segment(x = resFile$n[which(resFile[,i] > 0.8)[1]], y = 0,
                          xend = resFile$n[which(resFile[,i] > 0.8)[1]], yend = 0.8, colour = "black", linetype = "dashed", linewidth = 0.1)

  }
  print(p)
}

