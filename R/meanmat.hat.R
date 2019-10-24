#' Estimation the Mean Matrix
#'
#' This function estimates the mean matrix.
#'
#' It is assumed that there are \code{nrow(datamat)} row variables and
#' \code{ncol(datamat)}/\code{N} column variables in \code{datamat}. Further,
#' \code{datamat} should be written in such a way that every
#' \code{ncol(datamat)}/\code{N} consecutive columns belong to the same subject
#' and the order of the column variables in each block is preserved across
#' subjects.
#'
#' @aliases meanmat.hat print.meanmat.hat
#' @param datamat numeric matrix containing the transposable data.
#' @param N positive integer number indicating the sample size, i.e., the
#' number of subjects.
#' @param group.sizes numeric vector indicating the size of the row or column
#' groups that share the same mean vector. It should be used only when
#' \code{group.vars=}'\code{rows}' or '\code{columns}'.
#' @param group.vars character indicating that the mean matrix can be
#' simplified over the row or column variables. Options include '\code{rows}'
#' or '\code{columns}'.
#' @return Returns a list with components: \item{estmeanmat}{the estimated mean
#' matrix.} \item{N}{the sample size.} \item{n.rows}{the number of row
#' variables.} \item{n.cols}{the number of column variables.}
#' @author Anestis Touloumis
#' @references  Touloumis, A., Marioni, J. C. and Tavare, S. (2016)
#' HDTD: Analyzing multi-tissue gene expression data.
#' \emph{Bioinformatics} \bold{32}, 2193--2195.
#' @examples
#' data(VEGFmouse)
#' ## The sample mean matrix of the VEGF mouse data.
#' sample_mean <- meanmat.hat(datamat = VEGFmouse, N = 40)
#' sample_mean
#' sample_mean$estmeanmat
#' @export
meanmat.hat <- function(datamat, N, group.sizes = NULL, group.vars = NULL) {
    if (!is.matrix(datamat))
        datamat <- as.matrix(datamat,
                             dimnames =
                                 list(rownames(datamat), colnames(datamat)))
    datamat <- na.omit(datamat)
    rowsnam <- rownames(datamat)
    colsnam <- colnames(datamat)
    N <- as.numeric(N)
    if (length(N) != 1 | ((N - round(N)) != 0) | (N <= 1))
        stop("'N' must be a positive integer number greater than 1")
    p1 <- nrow(datamat)
    p2 <- ncol(datamat) / N
    colsnam <- colsnam[seq_len(p2)]
    if ((p2 - round(p2)) != 0)
        stop("The number of column variables is not a positive integer number")
    ans <- sumdatamatrix(datamat, N)
    if (!is.null(group.sizes)) {
        group.sizes <- as.numeric(group.sizes)
        if (any(group.sizes <= 0 | (group.sizes - round(group.sizes)) != 0))
            stop("'group.sizes' must be a vector of positive integer numbers")
        if (group.vars != "columns" & group.vars != "rows")
            stop("'group.vars' must be either 'rows' or 'columns'")
        if (group.vars == "columns" & sum(group.sizes) != p2)
            stop("Total group sizes does not match number of column variables")
        if (group.vars == "rows" & sum(group.sizes) != p1)
            stop("Total group sizes does not match number of row variables")
        projmat <- projectionmatrix(cumsum(group.sizes))
        ans <- if (group.vars == "columns")
            ans %*% (diag(p2) - projmat) else (diag(p1) - projmat) %*% ans
    }
    ans <- ans / N
    dimnames(ans) <- list(rowsnam, colsnam)
    ans <- list(estmeanmat = ans, N = N, n.rows = p1, n.cols = p2)
    class(ans) <- "meanmat.hat"
    ans
}
