meanmat.ts.generic <- function(datamat, N, projmat, voi) {
    dims <- prod(dim(datamat))/N
    if (voi == "columns") 
        datamat <- transposedatamatrix(datamat, N)
    datamat <- projmat %*% datamat
    test_stats <- meanmatts_statistics(datamat, N)
    ans <- test_stats[1]/sqrt(2 * test_stats[2]/N/(N - 1))
    ans
}
