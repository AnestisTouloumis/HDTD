projmatrix <- function(x) {
    lenx <- length(x)
    if (lenx != 1) {
        ans <- diag(x[lenx])
        index <- c(0, x)
        for (i in (seq_len(lenx) + 1)) {
            index1 <- seq(index[i - 1] + 1, index[i])
            lenind1 <- length(index1)
            if (lenind1 != 1) {
                ans[index1, index1] <-
                    ans[index1, index1] - rowMeans(ans[index1, index1])
            } else ans[index1, index1] <- 0
        }
    } else {
        ans <- diag(x) - 1 / x
    }
    ans
}
