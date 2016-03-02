meanmat.hat <-
function(datamat,N,group.sizes=NULL,group.vars=NULL)
{
if(!is.matrix(datamat)) datamat <- as.matrix(datamat,dimnames=list(rownames(datamat),colnames(datamat)))
datamat <- na.omit(datamat)
rowsnam <- rownames(datamat)
colsnam <- colnames(datamat)
N <- as.numeric(N)
if(length(N)!=1 | ((N-round(N))!=0) | (N<=1)) 
   stop("'N' must be a positive integer number greater than 1")
p1 <- nrow(datamat)
p2 <- ncol(datamat)/N
colsnam <- colsnam[seq_len(p2)]
if((p2-round(p2))!=0) 
   stop("The number of column variables is not a positive integer number")
ans <- matrix(rowSums(matrix(datamat,p1*p2,N)),p1,p2)
if(!is.null(group.sizes)){
   group.sizes <- as.numeric(group.sizes)
   if(any(group.sizes<=0) | any((group.sizes-round(group.sizes))!=0))
      stop("'group.sizes' must be a vector of positive integer numbers")
   if(group.vars!="columns" & group.vars!="rows")
      stop("'group.vars' must be either 'rows' or 'columns'")
   if(group.vars=="columns" & sum(group.sizes)!=p2)
      stop("The number of variables in the column groups is less than the number of column variables") 
   if(group.vars=="rows" & sum(group.sizes)!=p1)
      stop("The number of variables in the row groups is less than the number of row variables") 
   projmat <- projmatrix(cumsum(group.sizes))
   ans <- if(group.vars=="columns") ans%*%(diag(p2)-projmat) else (diag(p1)-projmat)%*%ans
                        }
ans <- ans/N
dimnames(ans) <- list(rowsnam,colsnam)
ans <- list(estmeanmat=ans,N=N,n.rows=p1,n.cols=p2)
class(ans) <- "meanmat.hat"
ans
}
