#' Centering Transposable Data
#' 
#' This function centers the transposable data around their sample mean matrix.
#' 
#' It is assumed that there are \code{nrow(datamat)} row variables and
#' \code{ncol(datamat)}/\code{N} column variables in \code{datamat}. Further,
#' \code{datamat} should be written in such a way that every
#' \code{ncol(datamat)}/\code{N} consecutive columns belong to the same subject
#' and the order of the column variables in each block is preserved across
#' subjects.
#' 
#' @param datamat numeric matrix containing the transposable data.
#' @param N positive integer number indicating the sample size, i.e., the
#' number of subjects.
#' @return Returns a matrix of the same size as \code{datamat}.
#' @author Anestis Touloumis
#' @seealso \code{\link{covmat.hat}} and \code{\link{covmat.ts}}.
#' @examples
#' data(VEGFmouse)
#' ## Centering the VEGF dataset around the sample mean matrix.
#' VEGFcen <- centerdata(VEGFmouse,40)
#' @export centerdata
centerdata <- function(datamat, N) {
    if (!is.matrix(datamat)) 
        datamat <- as.matrix(datamat, dimnames = list(rownames(datamat), colnames(datamat)))
    datamat <- na.omit(datamat)
    N <- as.numeric(N)
    if (length(N) != 1 | ((N - round(N)) != 0) | (N <= 1)) 
        stop("'N' must be an integer number greater than 1")
    p1 <- nrow(datamat)
    p2 <- ncol(datamat)/N
    if ((p2 - round(p2)) != 0) 
        stop("The number of column variables is not a positive integer number")
    meanmat <- sumdatamatrix(datamat, N)/N
    datamat - matrix(meanmat, p1, N * p2)
}
