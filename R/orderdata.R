orderdata <-
function(datamat,N,order.rows=NULL,order.cols=NULL)
{
if(!is.matrix(datamat)) datamat <- as.matrix(datamat,dimnames=list(rownames(datamat),colnames(datamat)))
datamat <- na.omit(datamat)
row.vars.names <- rownames(datamat)
col.vars.names <- colnames(datamat)
N <- as.numeric(N)
if(length(N)!=1 | ((N-round(N))!=0) | (N<=0)) 
   stop("'N' must be a positive integer number")
p2 <- ncol(datamat)/N
if((p2-round(p2))!=0) 
   stop("The number of column variables is not a positive integer number")
if((is.null(order.rows) & is.null(order.cols)))
   stop("At least one of 'order.rows' and/or 'order.cols' must be provided")  
if(!is.null(order.rows)){
   if(!is.numeric(order.rows)) order.rows <- as.numeric(order.rows)
   if(length(unique(order.rows))>length(order.rows))
      stop("'order.rows' contains duplicated row variables")
   if(((order.rows-round(order.rows))!=0)|| any(order.rows<=0)) 
      stop("'order.rows' must be a vector of positive integer numbers")
   ans <- datamat[order.rows,]
   if(!is.null(row.vars.names)) rownames(ans) <- row.vars.names[order.rows]
                        } 
if(!is.null(order.cols)){
   if(!is.numeric(order.cols)) order.cols <- as.numeric(order.cols)
   if(length(unique(order.cols))>length(order.cols))
      stop("'order.cols' contains duplicated column variables")
   if(((order.cols-round(order.cols))!=0)|| any(order.cols<=0)) 
      stop("'order.cols' must be a vector of positive integer numbers")
   if(!is.null(order.rows)) ans <- transposedata(ans,N)[order.cols,] else ans <- transposedata(datamat,N)[order.cols,]
   ans <- transposedata(ans,N)
if(!is.null(col.vars.names)) colnames(ans) <- paste(col.vars.names[order.cols],rep(seq_len(N),each=ncol(ans)/N),sep=".")
                        }
  ans 
}
