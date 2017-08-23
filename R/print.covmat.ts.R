#' @export
print.covmat.ts <- function(x, ...) {
    cat("SPHERICITY AND IDENTITY TESTS FOR THE ROW OR COLUMN VARIABLES", 
        "\n")
    cat("Sample size           = ", x$N, "\n")
    cat("Row variables         = ", x$n.rows, "\n")
    cat("Column variables      = ", x$n.cols, "\n")
    cat("Variables tested      = ", x$variables, "\n")
    cat("Centered data         = ", x$centered, "\n")
    cat("\nSphericity test for the covariance matrix of the", 
        x$variables, "\n")
    if (x$sphericity.ts$p.value < 1e-04) 
        cat("Test Statistic =", round(x$sphericity.ts$statistic, 
            4), ", p-value < 0.0001\n")
    if (x$sphericity.ts$p.value > 1e-04) 
        cat("Test Statistic =", round(x$sphericity.ts$statistic, 
            4), ", p-value =", round(x$sphericity.ts$p.value, 
            4), "\n")
    cat("\nIdentity test for the covariance matrix of the", x$variables, 
        "\n")
    if (x$identity.ts$p.value < 1e-04) 
        cat("Test Statistic =", round(x$identity.ts$statistic, 
            4), ", p-value < 0.0001\n")
    if (x$identity.ts$p.value > 1e-04) 
        cat("Test Statistic =", round(x$identity.ts$statistic, 
            4), ", p-value =", round(x$identity.ts$p.value, 4), 
            "\n")
}
