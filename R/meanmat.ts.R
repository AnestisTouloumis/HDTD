meanmat.ts <-
  function(datamat,N,group.sizes,voi="columns")
  {
    if(!is.matrix(datamat)) datamat <- as.matrix(datamat)
    datamat <- na.omit(datamat)
    N <- as.numeric(N)
    if(length(N)!=1 | ((N-round(N))!=0) | (N<=3)) 
      stop("'N' must be a positive integer number greater than 3")
    group.sizes <- as.numeric(group.sizes)
    if(any(group.sizes<=0) || any((group.sizes-round(group.sizes))!=0))
      stop("'group.sizes' must be a vector of positive integer numbers")
    voi <- as.character(voi)
    if(voi!="columns" & voi!="rows")
      stop("'voi' must be either 'rows' or 'columns'")
    p1 <- nrow(datamat)
    p2 <- ncol(datamat)/N
    if((p2-round(p2))!=0) 
      stop("The number of column variables is not a positive integer number")
    if(voi=="columns" & sum(group.sizes)!=p2)
      stop("The total number of variables in the prespecified groups is less than the number of column variables") 
    if(voi=="rows" & sum(group.sizes)!=p1)
      stop("The total number of variables in the prespecified groups is less than the number of row variables") 
    if(all(group.sizes==1))
      stop("The size of at least one of the prespecified groups needs to be greater than one") 
    if(any(group.sizes==1)) {
      rmvars <- cumsum(group.sizes)[group.sizes==1]
      order.rows <- order.cols <- NULL
      if(voi=="rows") order.rows <- (1:p1)[-rmvars] else order.cols <- (1:p2)[-rmvars]
      datamat <- orderdata(datamat,N,order.rows,order.cols) 
      group.sizes1 <- group.sizes[group.sizes!=1]
      projmat <- projmatrix(cumsum(group.sizes1))
    } else  projmat <- projmatrix(cumsum(group.sizes))
    ts <-  meanmat.ts.generic(datamat,N,projmat,voi)
    ans <- list(statistic=ts,p.value=1-pnorm(ts),voi=voi,n.groups=length(group.sizes),group.sizes=group.sizes,
                N=N,n.rows=p1,n.cols=p2) 
    class(ans) <- "meanmat.ts"
    ans
  }
