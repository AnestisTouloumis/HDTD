#' @export
print.covmat.ts <- function(x, ...) {
    print_vars <- ifelse(x$variables == "Rows", "ROW", "COLUMN")
    cat("HYPOTHESES TESTS FOR THE",print_vars, "COVARIANCE MATRIX", "\n")
    cat("Sample size      =", x$N, "\n")
    cat("Row variables    =", x$n.rows, "\n")
    cat("Column variables =", x$n.cols, "\n")
    cat("Centered data    =", x$centered, "\n")
    cat("\nDiagonality hypothesis test:\n")
    if (x$diagonality.ts$p.value < 1e-04) 
        cat("Test Statistic = ", round(x$diagonality.ts$statistic, 4),
            ", p-value < 0.0001\n", sep = "")
    if (x$diagonality.ts$p.value > 1e-04) 
        cat("Test Statistic = ", round(x$diagonality.ts$statistic, 4),
            ", p-value = ", round(x$diagonality.ts$p.value, 4), "\n", sep = "")
    cat("\nSphericity hypothesis test:\n")
    if (x$sphericity.ts$p.value < 1e-04) 
        cat("Test Statistic = ", round(x$sphericity.ts$statistic, 4),
            ", p-value < 0.0001\n", sep = "")
    if (x$sphericity.ts$p.value > 1e-04) 
        cat("Test Statistic = ", round(x$sphericity.ts$statistic, 4),
            ", p-value = ", round(x$sphericity.ts$p.value, 4), "\n", sep = "")
    cat("\nIdentity hypothesis test:\n")
    if (x$identity.ts$p.value < 1e-04) 
        cat("Test Statistic = ", round(x$identity.ts$statistic, 4),
            ", p-value < 0.0001\n", sep = "")
    if (x$identity.ts$p.value > 1e-04) 
        cat("Test Statistic = ", round(x$identity.ts$statistic, 4),
            ", p-value = ", 
            round(x$identity.ts$p.value, 4), "\n", sep =" ")
}