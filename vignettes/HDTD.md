---
title: "Using `HDTD` to Analyze High-Dimensional Transposable Data: An Application in Genetics"
author: Touloumis, A., Marioni, J.C., and Tavaré, S.
date: "15 February 2018"
output: rmarkdown::html_vignette
vignette: >
    %\VignetteIndexEntry{HDTD to Analyze High-Dimensional Transposable Data}
    %\VignetteEngine{knitr::rmarkdown}
    %\VignetteEncoding{UTF-8}
---


# Introduction
The R/Bioconductor package `HDTD` is designed to analyze high-dimensional transposable data. The term transposable data implies the following structural information in the dataset: 

* the data for each sampling unit (e.g., subject/patient) can be written in a matrix,
* the rows and the columns in each matrix correspond to two distinct sets of variables,
* dependencies might occur among and/or between the row and column variables.

The term high-dimensional implies that the sample size (e.g., the number of subjects/patients) is a lot smaller than the total number of row and column variables.

Since the statistical methods implemented in `HDTD` were primarily motivated by studies in genetics, a microarray dataset is utilized to illustrate the functionality of the package. However, the use of `HDTD` is not limited to gene-expression experiments and we emphasize that it is suitable for analyzing datasets that satisfy the high-dimensional transposable data definition.

# Mouse aging dataset
`HDTD` includes a subset of the tissue study described in Zahn et al. (2007). This dataset contains expression levels for $40$ mice. For each mouse, the expression levels of $46$ genes that belong to the vascular endothelial growth factor signalling pathway were measured across $9$ tissues (adrenal gland, cerebrum, hippocampus, kidney, lung, muscle, spinal cord, spleen and thymus). The experimental design satisfies the definition of transposable data because: 

* the data for each mouse can be written in a matrix form, 
* the rows correspond to genes and the columns to multiple tissues, and 
* we expect dependencies among genes, among tissues and between genes and tissues. 

The dataset is formatted as a single matrix with rows the $46$ genes and columns the $9 \times 40=360$ tissues.


```r
library("HDTD")
data(VEGFmouse)
dim(VEGFmouse)
#> [1]  46 360
rownames(VEGFmouse)
#>  [1] "Akt1"     "Akt2"     "Akt3"     "Arnt"     "Casp9"    "Cdc42"   
#>  [7] "Grb2"     "Hif1a"    "Hras1"    "Hsp90aa1" "Hspb1"    "Map2k1"  
#> [13] "Map2k2"   "Mapk1"    "Mapk13"   "Mapk14"   "Mapk3"    "Mapkapk2"
#> [19] "Nfat5"    "Nfatc3"   "Nfatc4"   "Nos3"     "Nras"     "Nrp1"    
#> [25] "Pdgfc"    "Pik3ca"   "Pik3cb"   "Pik3cd"   "Pik3r1"   "Pik3r3"  
#> [31] "Pla2g12b" "Pla2g4a"  "Pla2g4c"  "Pla2g5"   "Pla2g6"   "Plcg1"   
#> [37] "Plcg2"    "Ppp3ca"   "Ppp3cb"   "Ppp3r1"   "Prkca"    "Ptk2"    
#> [43] "Rac1"     "Rac2"     "Raf1"     "Sphk2"
```

Further, every 9 consecutive columns belong to the same mouse and the tissues are ordered in the same way for each mouse. For example, we can check the column variables for the first two subjects:


```r
colnames(VEGFmouse)[1:18]
#>  [1] "adrenal.1"     "cerebrum.1"    "hippocampus.1" "kidney.1"     
#>  [5] "lung.1"        "muscle.1"      "spinal.1"      "spleen.1"     
#>  [9] "thymus.1"      "adrenal.2"     "cerebrum.2"    "hippocampus.2"
#> [13] "kidney.2"      "lung.2"        "muscle.2"      "spinal.2"     
#> [17] "spleen.2"      "thymus.2"
```

It is crucial to import the data in a particular format when using `HDTD`. In particular, first write the data for each subject (mouse) in a matrix form while preserving the order of the row (genes) and column (tissues) variables within each subject-specific matrix. Then, create a single matrix by stacking column-wise the subject-specific matrices the one after the other. Finally, read the datafile in `HDTD`.

# Mean relationship of genes across tissues
The user can determine the mean relationship of the genes across tissues by testing and estimating the mean matrix. One interesting hypothesis to be tested is conservation of the gene expression levels across tissues, i.e., if the mean gene expression levels vector in the VEGF signaling pathway changes across the $9$ tissues:

```r
meanmat.ts(VEGFmouse, N = 40, group.sizes = 9)
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
Since $p$-value $<0.001$, we have strong evidence against the null hypothesis that there is no tissue effect in the gene expression levels. To explore further the mean gene expression level pattern across tissues, additional tests can be carried out. We illustrate a more complicated hypothesis that also requires data manipulation. Consider testing the hypothesis that the mean gene expression levels vector is constant only across the adrenal glands, the spleen, the kidney and the lung. To do this, we first need to place these $4$ tissues in a successive order in the dataset,   

```r
colnames(VEGFmouse)[1:9]
#> [1] "adrenal.1"     "cerebrum.1"    "hippocampus.1" "kidney.1"     
#> [5] "lung.1"        "muscle.1"      "spinal.1"      "spleen.1"     
#> [9] "thymus.1"
columnorder <- c(1, 4, 5, 8, 2, 3, 6, 7, 9)
VEGForder <- orderdata(VEGFmouse, N = 40, order.cols = columnorder) 
colnames(VEGForder)[1:9]
#> [1] "adrenal.1.1"     "kidney.1.1"      "lung.1.1"        "spleen.1.1"     
#> [5] "cerebrum.1.1"    "hippocampus.1.1" "muscle.1.1"      "spinal.1.1"     
#> [9] "thymus.1.1"
```
and then to perform the test using the ordered dataset

```r
meanmat.ts(VEGForder, N = 40, group.sizes = c(4, 1, 1, 1, 1, 1))
#> MEAN MATRIX TEST 
#> Sample size      = 40 
#> Row variables    = 46 
#> Column variables = 9 
#> 
#> H_0: 6 groups of columns with a constant mean vector within each group 
#> H_1: not H_0 
#> 
#> The number of columns in the 6 successive groups are 4, 1, 1, 1, 1 and 1 respectively.
#> 
#> Test statistic = 218.5071, p-value < 0.0001
```
Note that we have included $5$ additional column groups of size one in the `group.sizes` argument to reflect the fact that the mean gene expression levels vector in each of the remaining $5$ tissues remained unspecified. The null hypothesis is rejected, and hence we may conclude that the mean gene expression levels vector is not constant in the adrenal glands, the spleen, the kidney and the lung.

Apart from hypothesis testing, one can estimate the mean relationship between the genes and the tissues. In this example, the mean matrix seems to be unstructured and thus the mean gene expression levels in the $9$ tissues can be estimated via the sample mean matrix 

```r
sample_mean <- meanmat.hat(VEGFmouse, N = 40)
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
Note that the output preserves the order of the genes and the columns. For example, $0.8399$ is the average $\log_2$ intensity for gene "Akt1" in the adrenal gland based on $40$ mice. The mean matrix for the first 10 genes across the $9$ tissues is

```r
head(round(sample_mean$estmeanmat, 4), n = 10)
#>          adrenal.1 cerebrum.1 hippocampus.1 kidney.1  lung.1 muscle.1
#> Akt1        0.8399     1.2157        1.0597   1.1469  1.2673   0.8459
#> Akt2       -0.2333    -0.6201       -0.3881  -0.5524 -0.5359  -0.2082
#> Akt3       -1.0856    -0.4351       -0.5490  -0.2534 -0.6091  -1.1794
#> Arnt        0.1089     0.1898        0.0968   0.2551 -0.1171  -0.2919
#> Casp9       0.0877     0.2600        0.4812   0.3203  0.7416  -0.2492
#> Cdc42      -0.0538     0.1657       -0.1516  -0.0548  0.0254  -0.3577
#> Grb2       -0.2765    -0.5322        0.0948   0.0162 -0.1499   0.3015
#> Hif1a      -0.5760     1.3233       -3.5652   1.5485  1.0256  -0.1264
#> Hras1      -0.8040    -0.5952       -0.4063  -0.4964 -0.4997  -0.4900
#> Hsp90aa1   -0.3007     0.1292        0.7474   0.5589  0.2163   0.2634
#>          spinal.1 spleen.1 thymus.1
#> Akt1       1.2201   1.2142   1.2025
#> Akt2      -0.3877  -0.5154  -0.6836
#> Akt3      -0.7467  -0.8073  -0.4024
#> Arnt      -0.0861  -0.5188  -0.2219
#> Casp9      0.3691   0.5682   0.0726
#> Cdc42     -0.0816  -0.0888   0.1986
#> Grb2       0.2385  -0.4664  -0.3850
#> Hif1a      1.1011  -2.5046   0.3866
#> Hras1     -0.4475  -0.8511  -0.7226
#> Hsp90aa1   0.7414  -0.2552  -0.2784
```

# Dependence structure of the genes and of the tissues
The user can estimate two covariance matrices, one for the genes (rows) and the other for the multiple tissues (columns). We have developed shrinkage estimators for both covariance matrices but we let the user decide if shrinkage is required to both, one or neither of these matrices. In principle, we recommend shrinking both covariance matrices in order to obtain well-defined and invertible covariance matrix estimators. The `covmat.hat` function provides the corresponding covariance estimators: 

```r
est_cov_mat <- covmat.hat(datamat = VEGFmouse, N = 40)
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
The output summarizes the results but the user can recover the full covariance matrix estimators. For example, the covariance matrix of the tissues is

```r
round(est_cov_mat$cols.covmat, 3)
#>               adrenal.1 cerebrum.1 hippocampus.1 kidney.1 lung.1 muscle.1
#> adrenal.1         0.037     -0.001         0.000   -0.001  0.001    0.001
#> cerebrum.1       -0.001      0.043         0.000    0.000 -0.003    0.001
#> hippocampus.1     0.000      0.000         0.027    0.002  0.000    0.000
#> kidney.1         -0.001      0.000         0.002    0.032  0.001   -0.001
#> lung.1            0.001     -0.003         0.000    0.001  0.081    0.001
#> muscle.1          0.001      0.001         0.000   -0.001  0.001    0.030
#> spinal.1          0.000     -0.001         0.000    0.000  0.016    0.001
#> spleen.1         -0.001     -0.004         0.000   -0.001  0.001    0.000
#> thymus.1         -0.001      0.000         0.001    0.003  0.001    0.000
#>               spinal.1 spleen.1 thymus.1
#> adrenal.1        0.000   -0.001   -0.001
#> cerebrum.1      -0.001   -0.004    0.000
#> hippocampus.1    0.000    0.000    0.001
#> kidney.1         0.000   -0.001    0.003
#> lung.1           0.016    0.001    0.001
#> muscle.1         0.001    0.000    0.000
#> spinal.1         0.039   -0.001    0.001
#> spleen.1        -0.001    0.039   -0.001
#> thymus.1         0.001   -0.001    0.041
```
Moreover, the user can study the gene-wise or tissue-wise correlation by using the `covmat.ts` function. For example, the results from the identity, sphericity and diagonality hypothesis tests applied on the column variables

```r
covmat.ts(datamat = VEGFmouse, N = 40, voi = "columns")
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
suggest that the tissues might be uncorrelated at a $5\%$ significance level. 

# How to cite

```r
print(citation("HDTD"), bibtex = TRUE)
#> 
#> Please use the following guidelines for citing `HDTD' in
#> publication:
#> 
#> To cite the mean matrix hypothesis testing methodology, please use
#> 
#>   Touloumis, A., Tavar\'{e}, S. and Marioni, J.C. (2015). Testing
#>   the Mean Matrix in High-Dimensional Transposable Data.
#>   Biometrics 71 (1), 157-166
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Article{,
#>     title = {Testing the Mean Matrix in High-Dimensional Transposable Data},
#>     author = {Anestis Touloumis and Simon Tavar\'{e} and John C. Marioni},
#>     journal = {Biometrics},
#>     year = {2015},
#>     volume = {71},
#>     issue = {1},
#>     pages = {157--166},
#>     url = {http://onlinelibrary.wiley.com/doi/10.1111/biom.12257/full},
#>   }
#> 
#> To cite the covariance matrix hypothesis testing methodology,
#> please use
#> 
#>   Touloumis, A., Marioni, J.C. and Tavar\'{e}, S. (2017).
#>   Hypothesis Testing for the Covariance Matrix in High-Dimensional
#>   Transposable Data with Kronecker Product Dependence Structure.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Article{,
#>     title = {Hypothesis Testing for the Covariance Matrix in 
#>          High-Dimensional Transposable Data with Kronecker Product Dependence Structure},
#>     author = {Anestis Touloumis and John C. Marioni and Simon Tavar\'{e}},
#>     journal = {Submitted},
#>     year = {2017},
#>   }
#> 
#> To cite HDTD or the estimation method for the covariance matrices,
#> please use
#> 
#>   Touloumis, A., Marioni, J.C. and Tavar\'{e}, S. (2016). HDTD:
#>   Analyzing multi-tissue gene expression data. Bioinformatics 32
#>   (14), 2193-2195
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Article{,
#>     title = {HDTD: Analyzing multi-tissue gene expression data},
#>     author = {Anestis Touloumis and John C. Marioni and Simon Tavar\'{e}},
#>     journal = {Bioinformatics},
#>     year = {2016},
#>     volume = {32},
#>     issue = {14},
#>     pages = {2193--2195},
#>     url = {https://doi.org/10.1093/bioinformatics/btw224},
#>   }
```



