#' Estimation and Hypothesis Testing in High-Dimensional Transposable Data
#' 
#' The package HDTD offers functions to estimate and test the matrix parameters
#' of transposable data in high-dimensional settings.
#' 
#' The term transposable data refers to datasets that are structured in a
#' matrix form such that both the rows and columns correspond to variables of
#' interest. For example, consider microarray studies in genetics where
#' multiple RNA samples across different tissues are available per subject. In
#' this case, a data matrix can be created with row variables the genes, column
#' variables the tissues and measurements the corresponding expression levels.
#' 
#' The function \code{\link{meanmat.hat}} estimates the mean matrix of the
#' transposable data.
#' 
#' The mean relationship of the row and column variables can be tested using
#' the function \code{\link{meanmat.ts}}. The implemented test is nonparametric
#' and not seriously restricted by the dependence structure among and/or
#' between the row and column variables.
#' 
#' The function \code{\link{covmat.hat}} provides Stein-type shrinkage
#' estimators for the row covariance matrix and/or for the column covariance
#' matrix under a matrix-variate normal model.
#' 
#' The sphericity and identity hypothesis for the row or column covariance
#' matrix can be tested using the function \code{\link{covmat.ts}}. Both tests
#' are nonparametric, i.e., they do not rely on a normality assumption.
#' 
#' There are three utility functions that allow the user to change to
#' interchange the role of row and column variables
#' (\code{\link{transposedata}}), to center the transposable data
#' (\code{\link{centerdata}}) or to rearrange the order of the row and/or
#' column variables (\code{\link{orderdata}}).
#' 
#' @name HDTD-package
#' @aliases HDTD-package HDTD
#' @docType package
#' @author Anestis Touloumis, John Marioni, Simon Tavare.
#' 
#' Maintainer: Anestis.Touloumis <A.Touloumis@@brighton.ac.uk>
#' @references Touloumis, A., Tavare, S. and Marioni, J. C. (2015) Testing the
#' Mean Matrix in High-Dimensional Transposable Data. \emph{Biometrics} \bold{71},
#' 157--166
#' 
#' Touloumis, A., Marioni, J. C. and Tavare, S. (2016) HDTD: Analyzing 
#' multi-tissue gene expression data. \emph{Bioinformatics} \bold{32}, 2193--2195.
#' 
#' Touloumis, A., Marioni, J. C. and Tavare, S. (2013) Hypothesis Testing for
#' the Covariance Matrix in High-Dimensional Transposable Data with Kronecker
#' Product Dependence Structure. \url{http://arxiv.org/abs/1404.7684}.
#' @keywords package
#' @examples
#' data(VEGFmouse)
#' ## The sample mean matrix.
#' sample.mean <- meanmat.hat(VEGFmouse,40)
#' sample.mean
#' ## Testing if there is no tissue effect on the mean expression level for each gene.
#' tistest <- meanmat.ts(VEGFmouse,40,group.sizes=9,voi='columns')
#' tistest
#' # Estimating the covariance matrices of the genes (rows) and of the tissues (columns).
#' estcovmat <- covmat.hat(VEGFmouse,40,shrink='both',centered=FALSE)
#' estcovmat
#' ## Hypothesis tests for the covariance matrix of the genes (rows).
#' genestest <- covmat.ts(VEGFmouse,40,'rows',FALSE)
#' genestest
#' ## Hypothesis tests for the covariance matrix of the tissues (columns).
#' tissuestest <- covmat.ts(VEGFmouse,40,'columns',FALSE)
#' tissuestest
#' @useDynLib HDTD, .registration = TRUE  
#' @import Rcpp
#' @importFrom stats na.omit pnorm
NULL
