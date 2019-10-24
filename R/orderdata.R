#' Reordering Row and Column Variables
#'
#' This utility function rearranges the row and/or the column variables in a
#' desired order.
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
#' @param order.rows numeric vector displaying the desired order of the row
#' variables.
#' @param order.cols numeric vector displaying the desired order of the column
#' variables.
#' @return Returns a matrix of the same size as \code{datamat}.
#' @author Anestis Touloumis
#' @seealso \code{\link{meanmat.ts}} and \code{\link{meanmat.hat}}.
#' @examples
#' data(VEGFmouse)
#' set.seed(1)
#' tissuesold <-  colnames(VEGFmouse[ ,1:9])
#' ## Suppose that you want to order the tissues in the folowing order.
#' tissuesnew <- colnames(VEGFmouse[ ,1:9])[sample(9)]
#' tissuesnew
#' ## To do this, create a numeric vector with the desired order.
#' ordtis <- pmatch(tissuesnew, tissuesold)
#' VEGFmousenew <- orderdata(datamat = VEGFmouse, N = 40, order.cols = ordtis)
#' colnames(VEGFmousenew)[1:9]
#' @export orderdata
orderdata <- function(datamat, N, order.rows = NULL, order.cols = NULL) {
    if (!is.matrix(datamat))
        datamat <- as.matrix(datamat,
                             dimnames =
                                 list(rownames(datamat), colnames(datamat)))
    datamat <- na.omit(datamat)
    row.vars.names <- rownames(datamat)
    col.vars.names <- colnames(datamat)
    N <- as.numeric(N)
    if (length(N) != 1 | ((N - round(N)) != 0) | (N <= 0))
        stop("'N' must be a positive integer number")
    p2 <- ncol(datamat) / N
    if ((p2 - round(p2)) != 0)
        stop("The number of column variables is not a positive integer number")
    if ((is.null(order.rows) & is.null(order.cols)))
        stop("At least one of 'order.rows' or 'order.cols' must be provided")
    if (!is.null(order.rows)) {
        if (!is.numeric(order.rows))
            order.rows <- as.numeric(order.rows)
        if (length(unique(order.rows)) < length(order.rows))
            stop("'order.rows' contains duplicated row variables")
        if (any((order.rows - round(order.rows)) != 0) || any(order.rows <= 0))
            stop("'order.rows' must be a vector of positive integer numbers")
        ans <- datamat[order.rows, ]
        if (!is.null(row.vars.names))
            rownames(ans) <- row.vars.names[order.rows]
    }
    if (!is.null(order.cols)) {
        if (!is.numeric(order.cols))
            order.cols <- as.numeric(order.cols)
        if (length(unique(order.cols)) < length(order.cols))
            stop("'order.cols' contains duplicated column variables")
        if (any((order.cols - round(order.cols)) != 0) || any(order.cols <= 0))
            stop("'order.cols' must be a vector of positive integer numbers")
        if (!is.null(order.rows)) {
            ans <- transposedata(ans, N)[order.cols, ]
            } else {
                ans <- transposedata(datamat, N)[order.cols, ]
                }
        ans <- transposedata(ans, N)
        if (!is.null(col.vars.names))
            colnames(ans) <- paste(col.vars.names[order.cols],
                                    rep(seq_len(N), each = ncol(ans) / N),
                sep = ".")
    }
    ans
}
