covmat.ts.generic <- function(datamat, N, centered) {
    p1 <- nrow(datamat)
    p2 <- ncol(datamat)/N
    if (!centered) {
        test_statistics <- if (p1 <= p2) 
            statistics(datamat, N) else statistics_trans(datamat, N)
        } else test_statistics <- 
        if (p1 <= p2){
            statistics_centered(datamat, N)
            } else {
            statistics_trans_centered(datamat, N)
            }
    Ustat <- p1 * test_statistics[2]/(test_statistics[1]^2) - 1
    Vstat <- test_statistics[2] - 2 * test_statistics[1] + p1
    Dstat <- test_statistics[2] - test_statistics[3]
    trSigmaC2 <- test_statistics[4]/test_statistics[2]
    sigmaU <- (2/(N - 1)) * (trSigmaC2/(p2^2))
    sigmaV <- sigmaU * p1
    sigmaD <- sigmaU * test_statistics[3]
    Utest <- Ustat/sigmaU
    Vtest <- Vstat/sigmaV
    Dtest <- Dstat/sigmaD
    ans <- list(Utest = Utest, Vtest = Vtest, Dtest = Dtest)
    ans
}
