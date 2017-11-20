#' Nonparametric Tests for the Row or Column Covariance Matrix
#' 
#' Testing the sphericity, identity and diagonality hypotheses for 
#' the row or column covariance matrix.
#' 
#' It is assumed that there are \code{nrow(datamat)} row variables and
#' \code{ncol(datamat)}/\code{N} column variables in \code{datamat}. Further,
#' \code{datamat} should be written in such a way that every
#' \code{ncol(datamat)}/\code{N} consecutive columns belong to the same subject
#' and the order of the column variables in each block is preserved across
#' subjects.
#' 
#' The tests are nonparametric and thus robust to departures from the
#' matrix-variate normal model.
#' 
#' @aliases covmat.ts print.covmat.ts
#' @param datamat numeric matrix containing the transposable data.
#' @param N positive integer number indicating the sample size, i.e., the
#' number of subjects.
#' @param voi character indicating if the test should be applied on the row or
#' column covariance matrix. Options include '\code{rows}' or '\code{columns}'.
#' @param centered logical indicating if the transposable data are centered.
#' Options include \code{TRUE} or \code{FALSE}.
#' @return It returns a list with components: \item{sphericity.ts}{a list
#' containing the test statistic and p-value of the sphericity test.}
#' \item{identity.ts}{a list containing the test statistic and p-value of the
#' identity test.} \item{N}{the sample size.} \item{n.rows}{the number of row
#' variables.} \item{n.cols}{the number of column variables.}
#' \item{variables}{character indicating if the tests were applied to the row
#' or column covariance matrix.} \item{centered}{logical indicating if the
#' transposable data were centered.}
#' @author Anestis Touloumis
#' @references Touloumis, A., Marioni, J.C. and Tavare, S. (2017). Hypothesis
#' Testing for the Covariance Matrix in High-Dimensional Transposable Data with
#' Kronecker Product Dependence Structure.
#' @examples
#' data(VEGFmouse)
#' ## Hypothesis tests for the covariance matrix of the genes (rows).
#' genes_test <- covmat.ts(datamat = VEGFmouse, N = 40)
#' genes_test
#' ## Hypothesis tests for the covariance matrix of the tissues (columns).
#' tis_test <- covmat.ts(datamat = VEGFmouse, N = 40, voi = 'columns')
#' tis_test
#' @export covmat.ts
covmat.ts <- function(datamat = datamat, N = N, voi = "rows", 
                        centered = FALSE) {
    if (!is.matrix(datamat)) 
        datamat <- as.matrix(datamat)
    datamat <- na.omit(datamat)
    N <- as.numeric(N)
    if (length(N) != 1 | ((N - round(N)) != 0) | (N <= 0)) 
        stop("'N' must be a positive integer number")
    voi <- as.character(voi)
    if (voi != "columns" & voi != "rows") 
        stop("'voi' must be either 'rows' or 'columns'")
    centered <- as.logical(centered)
    if (centered != TRUE & centered != FALSE) 
        stop("'centered' must be either 'TRUE' or 'FALSE'")
    if (!centered & N <= 3) 
        stop("'N' must be greater than or equal to 4")
    p1 <- nrow(datamat)
    p2 <- ncol(datamat)/N
    if ((p2 - round(p2)) != 0) 
        stop("The number of column variables is not a positive integer number")
    if (voi == "columns") 
        datamat <- transposedata(datamat, N)
    ts <- covmat.ts.generic(datamat, N, centered)
    ans <- list(diagonality.ts = list(statistic = ts$Dtest,
                                        p.value = 1 - pnorm(ts$Dtest)), 
        sphericity.ts = list(statistic = ts$Utest, 
                                p.value = 1 - pnorm(ts$Utest)), 
        identity.ts = list(statistic = ts$Vtest, 
                                p.value = 1 - pnorm(ts$Vtest)), 
        N = N, n.rows = p1, n.cols = p2, 
        variables = if (voi == "rows") "Rows" else "Columns", 
        centered = centered)
    class(ans) <- "covmat.ts"
    ans
}
