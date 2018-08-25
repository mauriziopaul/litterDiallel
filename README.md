# Litter Diallel

## Summary

This package has functions that are useful for reproducing the analysis in our litter diallel manuscript.

## Installation

You can install `litterDiallel` using the following steps. Please make sure `devtools` is installed in R.

1. Install dependency in R:

```
install.packages("MCMCglmm")
```

2. Install `litterDiallel`:

```
devtools::install_github("mauriziopaul/litterDiallel")
```

## Notes

It may also be useful to install BayesDiallel. The instructions are [here](http://valdarlab.unc.edu/software/bayesdiallel/BayesDiallel.html).

You should then be able to load the package using `library(litterDiallel)` in R.

To load the data set, `litters`, use:

```
data("litters")
```
