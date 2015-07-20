covmat.hat <-
function(datamat,N,shrink="both",centered=FALSE,voi="both")
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
voi <- as.character(voi)
vois <- c("rows","columns","both")
icheck <- as.integer(match(voi,vois, -1))
if(icheck<1)
  stop(" 'voi' must be one of 'rows', 'columns', 'both'")
p1 <- nrow(datamat)
p2 <- ncol(datamat)/N
if((p2-round(p2))!=0) 
   stop(" The number of column variables is not a positive integer number")
pars <- covmat.hat.generic(datamat,N,shrink,centered,p1,p2,voi)
shrink <- switch(shrink, "none" = "None", "both" = "Both sets of variables", "rows" = "Row variables", "columns" = "Column variables")
ans <- list(rows.covmat=pars$rowcovmat,rows.intensity=pars$lambdaS,
            cols.covmat=pars$colcovmat,cols.intensity=pars$lambdaD,
            N=N,n.rows=p1,n.cols=p2,shrink=shrink,centered=centered)
class(ans) <- "covmat.hat"
ans
}
