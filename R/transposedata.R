transposedata <-
function(datamat,N){
p1 <- nrow(datamat)
p2 <- ncol(datamat)/N
idold <- rep(1:N,each=p2)
idnew <- rep(1:N,each=p1)
ans <- matrix(0,p2,p1*N)
for(i in 1:N) ans[,idnew==i] <- t(datamat[,idold==i])
ans
}
