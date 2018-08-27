# Litter Diallel

## Summary

This package includes functions for reproducing the analysis in our litter diallel manuscript.

## Installation

You can install `litterDiallel` using the following steps. First, please make sure `devtools` is installed in R.

1. Install MCMCglmm (for compatibility, must be version 2.25) in R:

	```R
	devtools::install_version("MCMCglmm", version = "2.25", repos = "http://cran.us.r-project.org")
	```

2. Install `litterDiallel`:

	```R
	devtools::install_github("mauriziopaul/litterDiallel")
	```

## Using the Package

3. You should then be able to load the package in R:

	```R
	library(litterDiallel)
	``` 

4. To load the data set, `litters`, use:

	```R
	data("litters")
	```

5. For an overview of the analysis, see the vignette:

	```R
	vignette("my-vignette")
	```

## Notes

- [ ] For the model inclusion probability analysis, BayesDiallel and BayesSpike must be installed. The instructions are [here](http://valdarlab.unc.edu/software/bayesdiallel/BayesDiallel.html) and [here](http://valdarlab.unc.edu/software/bayesdiallel/BayesDiallel.html).

- [ ] For some simple, miscellaneous functions useful for this analysis, install the following package, `PLMcctools`:

	```R
	devtools::install_github("mauriziopaul/PLMcctools")
	```

- [ ] To install other suggested packages, use:

	```R
	install.packages(c("tools", "data.table", "xtable"))
	```
