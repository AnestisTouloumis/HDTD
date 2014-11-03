transposedata <-
function(datamat,N){
if(!is.matrix(datamat)) datamat <- as.matrix(datamat,dimnames=list(rownames(datamat),colnames(datamat)))
datamat <- na.omit(datamat)
row.vars.names <- rownames(datamat)
col.vars.names <- colnames(datamat)
N <- as.numeric(N)
if(length(N)!=1 | ((N-round(N))!=0) | (N<=0)) 
  stop("'N' must be a positive integer number")
p1 <- nrow(datamat)
p2 <- ncol(datamat)/N
if((p2-round(p2))!=0) 
  stop("The number of column variables in `datamat' is not a positive integer number")
idold <- rep(1:N,each=p2)
idnew <- rep(1:N,each=p1)
ans <- matrix(0,p2,p1*N)
for(i in 1:N) ans[,idnew==i] <- t(datamat[,idold==i])
rownames(ans) <- col.vars.names[1:p2]
colnames(ans) <- paste(row.vars.names,rep(1:N,each=ncol(ans)/N),sep=".")
ans
}
