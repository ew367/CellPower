# BrainPowR

The majority of epigenetic epidemiology studies to date have generated genome-wide profiles from bulk tissues (e.g. brain) however these are vulnerable to confounding from variation in cellular composition. FANSs technology now allows cells sorting prior to array-based methylation profiling, potentially allowing smaller differences to be detected. 

We have developed a package that performs power analyses on cell sorted array data, which calculates either the necessary sample size to detect a specified mean difference in methylation (```calcSamples()``` function), or the level of difference that can be detected using a specified sample size (```calcDiff()```).

This package contains functions to calculate and plot power curves for cell specific methylation array data, building on the functionality of the [pwr package](https://cran.r-project.org/web/packages/pwr/index.html)



# Installation

The BrainPower package is available via GitHub and can be installed using the devtools package. However, there are some pre-requistite packages that may need to be installed from Bioconductor first.

```install_github("ew367/BrainPowR")
```


# Quick Start

Once you have successfully installed the BrainPower package you can calculate either the number of samples needed to detect a specified mean difference or the level of difference detectable from the number of samples available, and the corresponding power curves. The power calculation functions can perform the analyses using either a matrix of normalised beta values and a corresponding phenotype file, or a matrix of standard deviations derived from the betas matrix.



Within the BrainPower package we have provided an example matrix of standard deviations for cell sorted brain data, which can be used as follows:

```
library(BrainPower)

data(cellTypeSDswithProbeIDs.rdat)

# to calculate the number of samples needed to detect a specified difference
allSamples <- calcSamples(allSDs, meanDiff = 5, dataType = "SDs")

# to calculate the level of difference detectable from the number of samples available 
allSamples <- calcDiff(allSDs, nSamples = 100, dataType = "SDs")

# to calculate the proportion of probes with power > 0.8
allProps <-calcProps(allSamples)


# to plot the power curve 
myPlot <- plotPower(allProps, "samples")
or
myPlot <- plotPower(allProps, "difference")
```


A longer tutorial that explores the full functionality of the BrainPower package is available here


# Acknowledgements
We are grateful to the developers and contributors of the pwr package. https://CRAN.R-project.org/package=pwr  
