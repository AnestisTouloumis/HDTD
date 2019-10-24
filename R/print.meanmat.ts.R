#' @export
print.meanmat.ts <- function(x, ...) {
    cat("MEAN MATRIX TEST", "\n")
    cat("Sample size      =", x$N, "\n")
    cat("Row variables    =", x$n.rows, "\n")
    cat("Column variables =", x$n.cols, "\n")
    if (x$n.groups > 1) {
        cat("\nH_0:", x$n.groups, "groups of", x$voi,
            "with a constant mean vector within each group \n")
    } else {
        cat("\nH_0: a constant mean vector across", x$voi,
            "\n")
    }
    cat("H_1: not H_0 \n")
    if (x$n.groups == 2)
        cat("\nThe number of ", x$voi,
            " in the two successive groups are ",
            paste(x$group.sizes, collapse = " and "),
            " respectively.\n", sep = "")
    if (x$n.groups > 2)
        cat("\nThe number of ", x$voi,
            " in the ", x$n.groups, " successive groups are ",
            paste(c(paste(x$group.sizes[1:(x$n.groups - 1)],
                    collapse = ", "), x$group.sizes[x$n.groups]),
                    collapse = " and "), " respectively.\n",
            sep = "")
    if (x$p.value < 1e-04)
        cat("\nTest statistic = ", round(x$statistic, 4),
            ", p-value < 0.0001\n", sep = "")
    if (x$p.value >= 1e-04)
        cat("\nTest statistic = ", round(x$statistic, 4),
            ", p-value = ", round(x$p.value, 4), "\n",
            sep = "")
}
