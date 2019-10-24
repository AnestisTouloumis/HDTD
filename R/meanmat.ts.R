#' Nonparametric Tests for the Mean Matrix
#'
#' This function performs hypothesis testing for the mean matrix.
#'
#' It is assumed that there are \code{nrow(datamat)} row variables and
#' \code{ncol(datamat)}/\code{N} column variables in \code{datamat}. Further,
#' \code{datamat} should be written in such a way that every
#' \code{ncol(datamat)}/\code{N} consecutive columns belong to the same subject
#' and the order of the column variables in each block is preserved across
#' subjects.
#'
#' @aliases meanmat.ts print.meanmat.ts
#' @param datamat numeric matrix containing the transposable data.
#' @param N positive integer number indicating the sample size, i.e., the
#' number of subjects.
#' @param group.sizes numeric vector indicating the group sizes under the null
#' hypothesis.
#' @param voi character indicating if the test will be applied to the row or
#' column variables. Options include '\code{rows}' or '\code{columns}'.
#' @return Returns a list with components: \item{statistic}{the value of the
#' test statistic.} \item{p.value}{the corresponding p-value.} \item{voi}{the
#' set of variables that the test was applied to.} \item{n.groups}{the number
#' of groups under the null hypothesis.} \item{group.sizes}{the size of each
#' group under the null hypothesis.} \item{N}{the sample size.}
#' \item{n.rows}{the number of row variables.} \item{n.cols}{the number of
#' column variables.}
#' @author Anestis Touloumis
#' @references Touloumis, A., Tavare, S. and Marioni, J. C. (2015) Testing the
#' Mean Matrix in High-Dimensional Transposable Data.
#' \emph{Biometrics} \bold{71}, 157--166.
#'
#' @examples
#' data(VEGFmouse)
#' ## Testing conservation of the overall gene expression across tissues.
#' tissues_mean_test <- meanmat.ts(datamat = VEGFmouse, N = 40, group.sizes = 9)
#' tissues_mean_test
#' ## Testing if the adrenal and the cerebrum tissues have the same mean vector.
#' test2 <- meanmat.ts(VEGFmouse, N = 40, group.sizes = c(2, rep(1,7)))
#' test2
#' @export
meanmat.ts <- function(datamat, N, group.sizes, voi = "columns") {
    if (!is.matrix(datamat))
        datamat <- as.matrix(datamat)
    datamat <- na.omit(datamat)
    N <- as.numeric(N)
    if (length(N) != 1 | ((N - round(N)) != 0) | (N <= 3))
        stop("'N' must be a positive integer number greater than 3")
    group.sizes <- as.numeric(group.sizes)
    if (any(group.sizes <= 0) || any((group.sizes - round(group.sizes)) != 0))
        stop("'group.sizes' must be a vector of positive integer numbers")
    voi <- as.character(voi)
    if (voi != "columns" & voi != "rows")
        stop("'voi' must be either 'rows' or 'columns'")
    p1 <- nrow(datamat)
    p2 <- ncol(datamat) / N
    if ((p2 - round(p2)) != 0)
        stop("The number of column variables is not a positive integer number")
    if (voi == "columns" & sum(group.sizes) != p2)
        stop("Total group sizes is less than the number of column variables")
    if (voi == "rows" & sum(group.sizes) != p1)
        stop("Total group sizes is less than the number of row variables")
    if (all(group.sizes == 1))
        stop("The size of at least one group must be greater than one")
    if (any(group.sizes == 1)) {
        rmvars <- cumsum(group.sizes)[group.sizes == 1]
        order.rows <- order.cols <- NULL
        if (voi == "rows")
            order.rows <- (1:p1)[-rmvars] else order.cols <- (1:p2)[-rmvars]
        datamat <- orderdata(datamat, N, order.rows, order.cols)
        group.sizes1 <- group.sizes[group.sizes != 1]
        projmat <- projectionmatrix(cumsum(group.sizes1))
    } else projmat <- projectionmatrix(cumsum(group.sizes))
    ts <- meanmat.ts.generic(datamat, N, projmat, voi)
    ans <- list(statistic = ts, p.value = 1 - pnorm(ts), voi = voi,
                n.groups = length(group.sizes), group.sizes = group.sizes,
                N = N, n.rows = p1, n.cols = p2)
    class(ans) <- "meanmat.ts"
    ans
}
