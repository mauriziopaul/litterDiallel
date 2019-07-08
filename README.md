# Litter Diallel

## Summary

This package includes functions for reproducing the analysis in our litter diallel manuscript.

## Installation

You can install `litterDiallel` using the following steps. First, please make sure `devtools` is installed in R.

1. Install MCMCglmm (for compatibility, must be version 2.25) in R:

    ```R
    devtools::install_version("MCMCglmm", version = "2.25", repos = "http://cran.us.r-project.org")
    ```

2. Install `litterDiallel` (with the vignette, _recommended_):

    ```R
    devtools::install_github("mauriziopaul/litterDiallel", build_vignettes=TRUE, build = TRUE, build_opts = c("--no-resave-data", "--no-manual"))
    ```
    
    or
    
    Install `litterDiallel` (without the vignette).

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
    browseVignettes("litterDiallel")
    ```

## Analysis Summary

The included data file, `litter-diallel.csv`, has the following column names, where each of 4,448 observed litters is represented by a single row. All counts represent animals that are observed at the time of weaning (~ 3 weeks after birth):

- `Dam_Founder` (factor): The inbred mouse strain name of the female parent, where the following abbreviations are used: A/J (AJ), C57BL/6J (B6), 129S1/SvImJ (129S1), NOD/ShiLtJ (NOD), NZO/HlLtJ (NZO), CAST/EiJ (CAST), PWK/PhJ (PWK), and WSB/EiJ (WSB).
- `Sire_Founder` (factor): The inbred mouse strain name of the male parent.
- `PupGeno` (factor): The strain-cross name of the F1 offspring, given as their dam-by-sire cross (strainDam x strainSire).
- `WeanDate` (factor): The date (DD/MM/YY) that animals were weaned into new cages and separated from their parents. 
- `YearMonth` (factor): The date (YearMonth-YYYY-MM) of weaning.
- `litterorder` (factor): The parity, or litter birth order based on dam, is provided as a factor.
- `litternum` (integer): The litter number (litterorder - 1) is provided as an integer, where the first litter is 0, and each subsequent litter is numbered starting at 1.
- `First_Litter` (integer): A binary variable taking the values 0 or 1, indicating whether it is the first litter born to the given dam.
- `Males` (integer): The number of male pups in the litter.
- `Females` (integer): The number of female pups in the litter.
- `Weaned` (integer): The total number of (male and female) pups in the litter.
- `Male_Prop` (numeric): The proportion of male pups to overall pups in the litter.

The order of the columns in the data set does not matter.

We use the functions `diallelMatrixMaker` and `diallelMatrixMakeAndRotate` to generate design matrices for modeling the different classes of effects. After reading in the data, these functions expect: the name of the data frame object, the dam column name, the sire column name, and two random effect (`batch`, `batch.1`) column names. 

We then use the `MCMCglmm` function from the [`MCMCglmm` package](https://github.com/cran/MCMCglmm) (version 2.25) to analyze our data by fitting (generalized) linear mixed models, including the overdispersed zero-truncated Poisson, binomial, or Gaussian model.

## Notes

- [ ] For the model inclusion probability analysis, BayesDiallel and BayesSpike must be installed. The instructions are [here](https://valdarlab.unc.edu/software/bayesdiallel/).

- [ ] For some simple, miscellaneous functions useful for this analysis, install the following package, `PLMcctools`:

    ```R
    devtools::install_github("mauriziopaul/PLMcctools")
    ```

- [ ] To install other suggested packages, use:

    ```R
    install.packages(c("tools", "data.table", "xtable"))
    ```
