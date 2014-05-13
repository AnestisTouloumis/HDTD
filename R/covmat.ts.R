covmat.ts <-
function(datamat=datamat,N=N,voi="rows",centered=FALSE)
{
if(!is.matrix(datamat)) datamat <- as.matrix(datamat)
datamat <- na.omit(datamat)
N <- as.numeric(N)
if(length(N)!=1 | ((N-round(N))!=0) | (N<=0)) 
   stop("'N' must be a positive integer number")
voi <- as.character(voi)
if(voi!="columns" & voi!="rows")
   stop("'voi' must be either 'rows' or 'columns'")
centered <- as.logical(centered)
if(centered!=TRUE & centered!=FALSE)
   stop("'centered' must be either 'TRUE' or 'FALSE'")
if(!centered & N<=3)
  stop("'N' must be greater than or equal to 4")
p1 <- nrow(datamat)
p2 <- ncol(datamat)/N
if((p2-round(p2))!=0) 
   stop("The number of column variables is not a positive integer number")
if(voi=="rows")  datamat <- transposedata(datamat,N) 
ts <-  covmat.ts.generic(datamat,N,centered)
ans <- list(sphericity.ts=list(statistic=ts$Utest,p.value=1-pnorm(ts$Utest)),
            identity.ts=list(statistic=ts$Vtest,p.value=1-pnorm(ts$Vtest)),
            N=N,n.rows=p1,n.cols=p2,variables=voi,centered=centered) 
class(ans) <- "covmat.ts" 
ans
}
