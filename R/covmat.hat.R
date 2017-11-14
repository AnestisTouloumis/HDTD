#' Estimation of the Row and of the Column Covariance Matrices.
#' 
#' This function provides the row and/or column covariance matrix estimators.
#' 
#' It is assumed that there are \code{nrow(datamat)} row variables and
#' \code{ncol(datamat)}/\code{N} column variables in \code{datamat}. Further,
#' \code{datamat} should be written in such a way that every
#' \code{ncol(datamat)}/\code{N} consecutive columns belong to the same subject
#' and the order of the column variables in each block is preserved across
#' subjects.
#' 
#' For identifiability reasons, the trace of the row covariance matrix is set
#' equal to its dimension. If you want to place the equivalent restriction on
#' the column covariance matrix, interchange the role of row and column
#' variables by utilizing the function \code{\link{transposedata}}.
#' 
#' @aliases covmat.hat print.covmat.hat
#' @param datamat numeric matrix containing the transposable data.
#' @param N positive integer number indicating the sample size, i.e., the
#' number of subjects.
#' @param shrink character indicating if shrinkage estimation should be
#' performed. Options include '\code{rows}', '\code{columns}', '\code{both}'
#' and '\code{none}'.
#' @param centered logical indicating if the transposable data are centered.
#' Options include \code{TRUE} or \code{FALSE}.
#' @param voi character indicating if the row, column or both covariance
#' matrices should be printed. Options include '\code{rows}', '\code{columns}'
#' and '\code{both}'.
#' @return Returns a list with components: \item{rows.covmat}{the estimated row
#' covariance matrix.} \item{rows.intensity}{the estimated row intensity.}
#' \item{cols.covmat}{the estimated column covariance matrix.}
#' \item{cols.intensity}{the estimated column intensity.} \item{N}{the sample
#' size.} \item{n.rows}{the number of row variables.} \item{n.cols}{the number
#' of column variables.} \item{shrink}{character indicating if shrinkage
#' estimation was performed.} \item{centered}{logical indicating if the
#' transposable data were centered.}
#' 
#' @references  Touloumis, A., Marioni, J. C. and Tavare, S. (2016) HDTD: Analyzing 
#' multi-tissue gene expression data, \emph{Bioinformatics} \bold{32}, 2193--2195.
#' @author Anestis Touloumis
#' @examples
#' data(VEGFmouse)
#' # Estimating the covariance matrices of the genes (rows) and of the tissues (columns).
#' estcovmat <- covmat.hat(VEGFmouse,40,shrink='both',centered=FALSE)
#' estcovmat
#' @export covmat.hat
covmat.hat <- function(datamat, N, shrink = "both", centered = FALSE, voi = "both") {
    if (!is.matrix(datamat)) 
        datamat <- as.matrix(datamat, dimnames = list(rownames(datamat), colnames(datamat)))
    datamat <- na.omit(datamat)
    N <- as.numeric(N)
    if (length(N) != 1 | ((N - round(N)) != 0) | (N <= 0)) 
        stop("'N' must be a positive integer number")
    shrink <- as.character(shrink)
    shrinks <- c("none", "rows", "columns", "both")
    icheck <- as.integer(match(shrink, shrinks, -1))
    if (icheck < 1) 
        stop("'shrink' must be one of 'rows', 'columns', 'none', 'both'")
    centered <- as.logical(centered)
    if (centered != TRUE & centered != FALSE) 
        stop("'centered' must be either 'TRUE' or 'FALSE'")
    if (!centered & N <= 3) 
        stop("'N' must be greater than or equal to 4")
    voi <- as.character(voi)
    vois <- c("rows", "columns", "both")
    icheck <- as.integer(match(voi, vois, -1))
    if (icheck < 1) 
        stop(" 'voi' must be one of 'rows', 'columns', 'both'")
    p1 <- nrow(datamat)
    p2 <- ncol(datamat)/N
    if ((p2 - round(p2)) != 0) 
        stop(" The number of column variables is not a positive integer number")
    pars <- covmat.hat.generic(datamat, N, shrink, centered, p1, p2, voi)
    shrink <- switch(shrink, none = "None", both = "Both sets of variables", rows = "Row variables", 
        columns = "Column variables")
    ans <- list(rows.covmat = pars$rowcovmat, rows.intensity = pars$lambdaS, cols.covmat = pars$colcovmat, 
        cols.intensity = pars$lambdaD, N = N, n.rows = p1, n.cols = p2, shrink = shrink, 
        centered = centered)
    class(ans) <- "covmat.hat"
    ans
}
