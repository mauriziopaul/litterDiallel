# Litter Diallel

## Summary

This package includes functions for reproducing the analysis in our litter diallel manuscript.

## Installation

You can install `litterDiallel` using the following steps. First, please make sure `devtools` is installed in R.

1. Install MCMCglmm (for compatibility, must be version 2.25) in R:

```
devtools::install_version("MCMCglmm", version = "2.25", repos = "http://cran.us.r-project.org")
```

2. Install `litterDiallel`:

```
devtools::install_github("mauriziopaul/litterDiallel")
```

You should then be able to load the package in R:

```
library(litterDiallel)
``` 

To load the data set, `litters`, use:

```
data("litters")
```

## Notes

A. For the model inclusion probability analysis, BayesDiallel and BayesSpike must be installed. The instructions are [here](http://valdarlab.unc.edu/software/bayesdiallel/BayesDiallel.html) and [here](http://valdarlab.unc.edu/software/bayesdiallel/BayesDiallel.html).

B. For some simple, miscellaneous functions useful for this analysis, install the following package, `PLMcctools`:

```
devtools::install_github("mauriziopaul/PLMcctools")
```

C. To install other suggested packages, use:

```
install.packages(c("tools", "data.table", "xtable"))
```
