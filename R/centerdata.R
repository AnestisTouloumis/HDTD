centerdata <-
function(datamat,N)
{
if(!is.matrix(datamat)) datamat <- as.matrix(datamat)
datamat <- na.omit(datamat)
N <- as.numeric(N)
if(length(N)!=1 | ((N-round(N))!=0) | (N<=1)) 
   stop("'N' must be an integer number greater than 1")
p1 <- nrow(datamat)
p2 <- ncol(datamat)/N
if((p2-round(p2))!=0) 
   stop("The number of column variables is not a positive integer number")
meanmat <- rowMeans(matrix(datamat,p1*p2,N))
datamat <- datamat-meanmat
datamat
}
