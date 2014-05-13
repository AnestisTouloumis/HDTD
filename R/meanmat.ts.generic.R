meanmat.ts.generic <-
function(datamat,N,projmat,voi)
{
dims <- prod(dim(datamat))/N
if(voi=="columns") datamat <- transposedata(datamat,N)
datamat <- projmat%*%datamat
datamat <- matrix(datamat,dims,N)
dataveccen <- datamat - rowMeans(datamat)
Gn <- helpcon <- 0
for(i in 1:(N-1)){
 Gn <- sum(crossprod(datamat[,i],datamat[,(i+1):N])) + Gn 
 helpcon <- sum(crossprod(dataveccen[,i],dataveccen[,(i+1):N])^2) + helpcon 
                 }
Gn <- 2*Gn/N/(N-1)
colssums <- colSums(dataveccen^2)
Q <- sum(colssums^2)
tr2S <- (sum(colssums)/(N-1))^2
trS2 <- (helpcon*2+Q)/(N-1)^2
Q <- Q/(N-1)
trOmega <- (N-1)/(N*(N-2)*(N-3))*((N-1)*(N-2)*trS2+tr2S-N*Q)
sdGn <- sqrt(2*trOmega/(N-1)/N)
ans <- Gn/sdGn
ans
}
