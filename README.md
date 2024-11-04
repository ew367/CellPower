# CellPower [![DOI](https://zenodo.org/badge/718615958.svg)](https://zenodo.org/doi/10.5281/zenodo.10125075)




Most epigenetic studies generate genome-wide profiles from bulk tissues (e.g. brain) however, this prevents identification of the cell type associated with the observed difference. Technologies such as [fluorescence-activated nuclei sorting (FANS)](https://dx.doi.org/10.17504/protocols.io.bmh2k38e) now allows cell sorting prior to array-based DNA methylation profiling and epigeneome-wide association studies to be performed within populations of specific cells. Not only does this enhance the biological interpretation of these studies, by studing a more homogeneous population of cells it likely increases the power of these studies. 

We have developed a package that performs power analyses for brain cell types and the Illumina EPIC array. We provide tools that can calculates either the necessary sample size to detect a specified mean difference in methylation (```calcSamples()``` function), or the minimum difference that can be detected using a specified sample size (```calcDiff()```). In both settings power calculations can be performed simulateously for multiple cell types to faciliate comparision. This package contains functions to calculate and plot power curves for cell specific methylation array data, building on the functionality of the [pwr package](https://cran.r-project.org/web/packages/pwr/index.html)

Within the package we provide reference data to enable users to perform calculations for brain cell types, the functions will also take user supplied standard deviations and this can be used to perform power calculations for DNA methylation analyses in other tissues or cell types, sublject to there being appropriate reference data.

# Installation

The CellPower package is available via GitHub and can be installed using the devtools package. However, there are some pre-requistite packages that may need to be installed
first.

```
load.lib<-c("doParallel", "dplyr", "foreach", "ggplot2", "pwr", "reshape2", "devtools") # necessary packages
install.lib <- load.lib[!load.lib %in% installed.packages()] # select packages not currently installed
for(lib in install.lib) install.packages(lib,dependencies=TRUE) # install the packages
sapply(load.lib,require,character=TRUE) # load the packages

library(devtools)

install_github("ew367/CellPower")
```
*code adapated from https://stackoverflow.com/questions/29041423/how-to-install-multiple-packages


# Quick Start

Once you have successfully installed the CellPower package you can calculate either the number of samples needed to detect a specified mean difference or the level of difference detectable from the number of samples available, and plot the corresponding power curves. The power calculation functions can perform the analyses using either a matrix of normalised beta values and a corresponding phenotype file, or a matrix of standard deviations derived from the betas matrix.



Within the CellPower package we have provided an example matrix of standard deviations for cell sorted brain data, which can be used as follows:

```
library(CellPower)

data(allSDs) # load example data

# to calculate the number of samples needed to detect a specified difference
allSamples <- calcSamples(allSDs, meanDiff = 5, dataType = "SDs")

# to calculate the proportion of probes with power > 0.8
allProps <-calcProps(allSamples)

# to plot the power curve 
myPlot <- plotPower(allProps, "samples")

# to calculate the level of difference detectable from the number of samples available 
allSamples <- calcDiff(allSDs, nSamples = 100, dataType = "SDs")

# to calculate the proportion of probes with power > 0.8
allProps <-calcProps(allSamples)

# to plot the power curve 

myPlot <- plotPower(allProps, "difference")
```


A longer tutorial that explores the full functionality of the CellPower package will be available shortly.


# Acknowledgements
We are grateful to the developers and contributors of the pwr package. https://CRAN.R-project.org/package=pwr  
