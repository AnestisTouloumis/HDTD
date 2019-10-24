covmat.hat.generic <- function(datamat, N, shrink, centered, p1, p2, voi) {
    lambdaS <- lambdaD <- 0
    datamat1 <- transposedatamatrix(datamat, N)
    if (!centered) {
        if (p1 < p2) {
            shrink_stats <- covmathat_statistics(datamat1, N)
            } else {
                shrink_stats <- covmathat_statistics_trans(datamat1, N)
                }
        } else {
            if (p1 < p2) {
                shrink_stats <- covmathat_statistics_centered(datamat1, N)
                } else {
                shrink_stats <- covmathat_statistics_trans_centered(datamat1, N)
                }
            }
    trDeltahat <- shrink_stats[1]
    trDelta2hat <- shrink_stats[2]
    trOmega2hat <- shrink_stats[3]
    trSigma2hat <- trOmega2hat / trDelta2hat
    if (centered) n <- N else n <- N - 1
    if (shrink == "both" | shrink == "columns") {
        lambdaD_a <- trSigma2hat / (p1 ^ 2) / n * (trDeltahat^2 + trDelta2hat)
        lambdaD_b <- trDelta2hat - trDeltahat ^ 2 / p2
        lambdaD <- lambdaD_a / (lambdaD_a + lambdaD_b)
        lambdaD <- max(0, min(1, lambdaD))
        }
    if (shrink == "both" | shrink == "rows") {
        lambdaS_a <- trDelta2hat / trDeltahat ^ 2 / n * (p1 ^ 2 + trSigma2hat)
        lambdaS_b <- trSigma2hat - p1
        lambdaS <- lambdaS_a / (lambdaS_a + lambdaS_b)
        lambdaS <- max(min(1, lambdaS), 0)
        }
    if (centered) {
        if (voi != "columns")
            Sigma <- tcrossprod(datamat) / N
        Delta <- tcrossprod(transposedatamatrix(datamat, N)) / N / p1
    } else {
        datacen <- centerdatamatrix(datamat, N)
        if (voi != "columns")
            Sigma <- matrix(0, p1, p1)
        Delta <- matrix(0, p2, p2)
        id <- rep(1:N, each = p2)
        for (i in seq_len(N)) {
            if (voi != "columns")
                Sigma <- Sigma + tcrossprod(datacen[, id == i])
            Delta <- Delta + crossprod(datacen[, id == i])
        }
        if (voi != "columns")
            Sigma <- Sigma / (N - 1)
        Delta <- Delta / (N - 1) / p1
    }
    rowcovmat <- colcovmat <- NULL
    if (voi != "columns") {
        Sigma <- Sigma / trDeltahat
        rowcovmat <- Sigma * (1 - lambdaS) + diag(lambdaS, p1)
    }
    if (voi != "rows") {
        colcovmat <- Delta * (1 - lambdaD) + diag(lambdaD * trDeltahat / p2, p2)
    }
    if ((shrink == "none" | shrink == "rows") | voi == "rows")
        lambdaD <- NULL
    if ((shrink == "none" | shrink == "columns") | voi == "columns")
        lambdaS <- NULL
    ans <- list(rowcovmat = rowcovmat, lambdaS = lambdaS,
                colcovmat = colcovmat, lambdaD = lambdaD)
    ans
}
