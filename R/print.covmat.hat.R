#' @export
print.covmat.hat <- function(x, ...) {
    cat("ESTIMATION OF THE ROW AND/OR COLUMN COVARIANCE MATRIX", "\n")
    cat("Sample size      =", x$N, "\n")
    cat("Row variables    =", x$n.rows, "\n")
    cat("Column variables =", x$n.cols, "\n")
    cat("Shrinking        =", x$shrink, "\n")
    cat("Centered data    =", x$centered, "\n")
    if (!is.null(x$rows.covmat)) {
        cat("\nROW VARIABLES\n")
        if (!is.null(x$rows.intensity)) 
            cat("Estimated optimal intensity =",
                round(x$rows.intensity, 4), "\n")
        cat("Estimated covariance matrix [1:5, 1:5] =\n")
        print(round(x$rows.covmat[1:min(5, x$n.rows), 1:min(5, x$n.rows)], 4))
    }
    if (!is.null(x$cols.covmat)) {
        cat("\nCOLUMN VARIABLES\n")
        if (!is.null(x$cols.intensity)) 
            cat("Estimated optimal intensity =", 
                round(x$cols.intensity, 4), "\n")
        cat("Estimated covariance matrix [1:5, 1:5] =\n")
        print(round(x$cols.covmat[1:min(5, x$n.cols), 1:min(5, x$n.cols)], 4))
    }
}
