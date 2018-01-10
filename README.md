
<!-- README.md is generated from README.Rmd. Please edit that file -->

# HDTD: Analyzing High-Dimensional Transposable Data

[![Travis-CI Build
Status](https://travis-ci.org/AnestisTouloumis/HDTD.svg?branch=master)](https://travis-ci.org/AnestisTouloumis/HDTD)
[![Project Status: Active The project has reached a stable, usable state
and is being actively
developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)

## Installation

You can install the release version of `HDTD`:

``` r
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("HDTD")
```

The source code for the release version of `HDTD` is available on
Bioconductor at:

  - <http://bioconductor.org/packages/release/bioc/html/HDTD.html>

Or you can install the development version of `HDTD`:

``` r
# install.packages('devtools')
devtools::install_github("AnestisTouloumis/HDTD")
```

To use `HDTD`, you should load the package as follows:

``` r
library("HDTD")
```

## Usage

This package offers functions to estimate and test the matrix parameters
of transposable data in high-dimensional settings. The term
*transposable data* refers to datasets that are structured in a matrix
form such that both the rows and columns correspond to variables of
interest and dependencies are expected to occur *among* rows, *among*
columns and *between* rows and columns. For example, consider microarray
studies in genetics where multiple RNA samples across different tissues
are available per subject. In this case, a data matrix can be created
with row variables the genes, column variables the tissues and
measurements the corresponding expression levels. We expect dependencies
to occur among genes, among tissues and between genes and tissues. For
more examples of trasposable data see references in Touloumis, Marioni
and Tavaré (2017), Touloumis, Tavaré and Marioni (2015) and Touloumis,
Marioni and Tavaré (2016).

There are four core functions:

  - `meanmat.hat` to estimate the mean matrix of the transposable data,
  - `meanmat.ts` to test the overall mean of the row (column) variables
    across groups of column (row) variables,
  - `covmat.hat` to estimate the row and column covariance matrix,
  - `covmat.ts` to test the sphericity, identity and diagonality
    hypothesis test for the row/column covariance matrix.

There are also three utility functions:

  - `transposedata` for interchanching the role of rows and columns,
  - `centerdata` for centering the transposable data around their mean
    matrix,
  - `orderdata` for rearranging the order of the row and/or column
    variables.

## Example

We replicate the analysis that can be found in the vignette based on the
mouse dataset

``` r
data(VEGFmouse)
```

This dataset contains expression levels for \(40\) mice. For each mouse,
the expression levels of \(46\) genes (rows) that belong to the vascular
endothelial growth factor signalling pathway were measured across \(9\)
tissues (adrenal gland, cerebrum, hippocampus, kidney, lung, muscle,
spinal cord, spleen and thymus) that are displayed in the columns.

One can estimate the mean relationship of the gene expression levels
across the \(9\) tissues

``` r
sample_mean <- meanmat.hat(datamat = VEGFmouse,N=40)
sample_mean
#> ESTIMATION OF THE MEAN MATRIX 
#> Sample size      = 40 
#> Row variables    = 46 
#> Column variables = 9 
#> 
#> Estimated mean matrix [1:5, 1:5] =
#>       adrenal.1 cerebrum.1 hippocampus.1 kidney.1  lung.1
#> Akt1     0.8399     1.2157        1.0597   1.1469  1.2673
#> Akt2    -0.2333    -0.6201       -0.3881  -0.5524 -0.5359
#> Akt3    -1.0856    -0.4351       -0.5490  -0.2534 -0.6091
#> Arnt     0.1089     0.1898        0.0968   0.2551 -0.1171
#> Casp9    0.0877     0.2600        0.4812   0.3203  0.7416
```

and test whether the overall gene expression is constant across the
\(9\) tissues:

``` r
tissue_mean_test <- meanmat.ts(datamat = VEGFmouse,N=40,group.sizes=9)
tissue_mean_test
#> MEAN MATRIX TEST 
#> Sample size      = 40 
#> Row variables    = 46 
#> Column variables = 9 
#> 
#> H_0: a constant mean vector across columns 
#> H_1: not H_0 
#> 
#> Test statistic = 373.5277, p-value < 0.0001
```

In this case, the overall gene expression is not conserved.

To analyze the gene-wise and tissue-wise dependence structure, one needs
to estimate the two covariance matrices:

``` r
est_cov_mat <- covmat.hat(datamat=VEGFmouse,N=40)
est_cov_mat
#> ESTIMATION OF THE ROW AND/OR COLUMN COVARIANCE MATRIX 
#> Sample size      = 40 
#> Row variables    = 46 
#> Column variables = 9 
#> Shrinking        = Both sets of variables 
#> Centered data    = FALSE 
#> 
#> ROW VARIABLES
#> Estimated optimal intensity = 0.0115 
#> Estimated covariance matrix [1:5, 1:5] =
#>          Akt1    Akt2    Akt3    Arnt   Casp9
#> Akt1   0.4139 -0.0248  0.0420 -0.0010  0.1084
#> Akt2  -0.0248  0.3341 -0.0240 -0.0029 -0.0151
#> Akt3   0.0420 -0.0240  0.6954  0.1733 -0.0168
#> Arnt  -0.0010 -0.0029  0.1733  0.4746  0.0850
#> Casp9  0.1084 -0.0151 -0.0168  0.0850  0.5337
#> 
#> COLUMN VARIABLES
#> Estimated optimal intensity = 0.3341 
#> Estimated covariance matrix [1:5, 1:5] =
#>               adrenal.1 cerebrum.1 hippocampus.1 kidney.1  lung.1
#> adrenal.1        0.0368    -0.0006        0.0001  -0.0006  0.0010
#> cerebrum.1      -0.0006     0.0432       -0.0002   0.0000 -0.0034
#> hippocampus.1    0.0001    -0.0002        0.0266   0.0019  0.0000
#> kidney.1        -0.0006     0.0000        0.0019   0.0317  0.0012
#> lung.1           0.0010    -0.0034        0.0000   0.0012  0.0809
```

Finally, the package allows users to perform hypothesis tests for the
covariance matrix of the genes

``` r
genes_cov_test <- covmat.ts(VEGFmouse,N=40)
genes_cov_test
#> HYPOTHESES TESTS FOR THE ROW COVARIANCE MATRIX 
#> Sample size      = 40 
#> Row variables    = 46 
#> Column variables = 9 
#> Centered data    = FALSE 
#> 
#> Diagonality hypothesis test:
#> Test Statistic = 8.6324, p-value < 0.0001
#> 
#> Sphericity hypothesis test:
#> Test Statistic = 132.8086, p-value < 0.0001
#> 
#> Identity hypothesis test:
#> Test Statistic = 30.3864, p-value < 0.0001
```

and of the tissues:

``` r
tissues_cov_test <- covmat.ts(VEGFmouse,N=40,voi="columns")
tissues_cov_test
#> HYPOTHESES TESTS FOR THE COLUMN COVARIANCE MATRIX 
#> Sample size      = 40 
#> Row variables    = 46 
#> Column variables = 9 
#> Centered data    = FALSE 
#> 
#> Diagonality hypothesis test:
#> Test Statistic = 1.4866, p-value = 0.0686
#> 
#> Sphericity hypothesis test:
#> Test Statistic = 10.2122, p-value < 0.0001
#> 
#> Identity hypothesis test:
#> Test Statistic = 38.0811, p-value < 0.0001
```

At a \(5\%\) significance level, it appears that the genes are
correlated but we do not have enough evidence to reject the hypothesis
that the tissues are uncorrelated.

## Getting help

The statistical methods implemented in `HDTD` are described in Touloumis
et al. (2017), Touloumis et al. (2015) and Touloumis et al. (2016).
Detailed examples of `HDTD` can be found in Touloumis et al. (2016) or
in the vignette:

``` r
browseVignettes("HDTD")
```

## How to cite

``` 

Please use the following guidelines for citing `HDTD' in
publication:

To cite the mean matrix hypothesis testing methodology, please use

  Touloumis, A., Tavar\'{e}, S. and Marioni, J.C. (2015). Testing
  the Mean Matrix in High-Dimensional Transposable Data.
  Biometrics 71 (1), 157-166

A BibTeX entry for LaTeX users is

  @Article{,
    title = {Testing the Mean Matrix in High-Dimensional Transposable Data},
    author = {Anestis Touloumis and Simon Tavar\'{e} and John C. Marioni},
    journal = {Biometrics},
    year = {2015},
    volume = {71},
    issue = {1},
    pages = {157--166},
    url = {http://onlinelibrary.wiley.com/doi/10.1111/biom.12257/full},
  }

To cite the covariance matrix hypothesis testing methodology,
please use

  Touloumis, A., Marioni, J.C. and Tavar\'{e}, S. (2017).
  Hypothesis Testing for the Covariance Matrix in High-Dimensional
  Transposable Data with Kronecker Product Dependence Structure.

A BibTeX entry for LaTeX users is

  @Article{,
    title = {Hypothesis Testing for the Covariance Matrix in 
         High-Dimensional Transposable Data with Kronecker Product Dependence Structure},
    author = {Anestis Touloumis and John C. Marioni and Simon Tavar\'{e}},
    journal = {Submitted},
    year = {2017},
  }

To cite HDTD or the estimation method for the covariance matrices,
please use

  Touloumis, A., Marioni, J.C. and Tavar\'{e}, S. (2016). HDTD:
  Analyzing multi-tissue gene expression data. Bioinformatics 32
  (14), 2193-2195

A BibTeX entry for LaTeX users is

  @Article{,
    title = {HDTD: Analyzing multi-tissue gene expression data},
    author = {Anestis Touloumis and John C. Marioni and Simon Tavar\'{e}},
    journal = {Bioinformatics},
    year = {2016},
    volume = {32},
    issue = {14},
    pages = {2193--2195},
    url = {https://doi.org/10.1093/bioinformatics/btw224},
  }
```

# References

<div id="refs" class="references">

<div id="ref-Touloumis2016">

Touloumis, A., Marioni, J.C. and Tavaré, S. (2016) HDTD: Analyzing
Multi-tissue Gene Expression Data. *Bioinfomatics*, **32**, 2193–2195.

</div>

<div id="ref-Touloumis2013">

Touloumis, A., Marioni, J.C. and Tavaré, S. (2017) Hypothesis Testing
for the Covariance Matrix in High-Dimensional Transposable Data with
Kronecker Product Dependence Structure. *Submitted*.

</div>

<div id="ref-Touloumis2015">

Touloumis, A., Tavaré, S. and Marioni, J.C. (2015) Testing the Mean
Matrix in High-Dimensional Transposable Data. *Biometrics*, **71**,
157–166.

</div>

</div>
