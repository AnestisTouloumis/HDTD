#' Interchanging the Row and Column Variables in Transposable Data
#' 
#' This function interchanges the row and column variables in transposable data
#' so that the original row variables will be treated as column variables and
#' the original column variables as row variables.
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
#' @return Returns a matrix with \code{ncol(datamat)} rows and
#' \code{nrow(datamat)}/\code{N} columns.
#' @author Anestis Touloumis
#' @seealso \code{\link{centerdata}} and \code{\link{orderdata}}.
#' @examples
#' data(VEGFmouse)
#' ## Transposing the VEGF dataset.
#' VEGFtr <- transposedata(VEGFmouse, N = 40)
#' @export transposedata
transposedata <- function(datamat, N) {
    if (!is.matrix(datamat)) 
        datamat <- as.matrix(datamat, dimnames = list(rownames(datamat), colnames(datamat)))
    datamat <- na.omit(datamat)
    row.vars.names <- rownames(datamat)
    col.vars.names <- colnames(datamat)
    N <- as.numeric(N)
    if (length(N) != 1 | ((N - round(N)) != 0) | (N <= 0)) 
        stop("'N' must be a positive integer number")
    p1 <- nrow(datamat)
    p2 <- ncol(datamat)/N
    if ((p2 - round(p2)) != 0) 
        stop("The number of column variables is not a positive integer number")
    ans <- transposedatamatrix(datamat, N)
    rownames(ans) <- col.vars.names[1:p2]
    colnames(ans) <- paste(row.vars.names, rep(seq_len(N), each = ncol(ans)/N), sep = ".")
    ans
}
