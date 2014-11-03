covmat.hat <-
function(datamat,N,shrink="both",centered=FALSE)
{
if(!is.matrix(datamat)) datamat <- as.matrix(datamat,dimnames=list(rownames(datamat),colnames(datamat)))
datamat <- na.omit(datamat)
N <- as.numeric(N)
if(length(N)!=1 | ((N-round(N))!=0) | (N<=0)) 
   stop("'N' must be a positive integer number")
shrink <- as.character(shrink)
shrinks <- c("none", "rows", "columns", "both")
icheck <- as.integer(match(shrink,shrinks, -1))
if(icheck<1)
   stop("'shrink' must be one of 'rows', 'columns', 'none', 'both'")
centered <- as.logical(centered)
if(centered!=TRUE & centered!=FALSE)
   stop("'centered' must be either 'TRUE' or 'FALSE'")
if(!centered & N<=3)
  stop("'N' must be greater than or equal to 4")
p1 <- nrow(datamat)
p2 <- ncol(datamat)/N
if((p2-round(p2))!=0) 
   stop("The number of column variables is not a positive integer number")
pars <- covmat.hat.generic(datamat,N,shrink,centered,p1,p2)
shrink <- switch(shrink, "none"="None", "both"="Both Sets of Variables", "rows"="Row Variables", "columns"="Column Variables")
ans <- list(rows.covmat=pars$rowcovmat,rows.intensity=pars$lambdaS,
            cols.covmat=pars$colcovmat,cols.intensity=pars$lambdaD,
            N=N,n.rows=p1,n.cols=p2,shrink=shrink,centered=centered)
class(ans) <- "covmat.hat"
ans
}
