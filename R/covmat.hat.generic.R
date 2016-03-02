covmat.hat.generic <-
  function(datamat,N,shrink,centered,p1,p2,voi)
  {
    id <- rep(seq_len(N),each=p2)
    if(centered) {
      if(voi!="columns")  Sigma <- tcrossprod(datamat)/N
      Delta <- tcrossprod(transposedata(datamat,N))/N/p1
    } else {
      datacen <- centerdata(datamat,N)
      if(voi!="columns") Sigma <- matrix(0,p1,p1)
      Delta <- matrix(0,p2,p2)
      for(i in seq_len(N)) {
        if(voi!="columns") Sigma <- Sigma + tcrossprod(datacen[,id==i])
        Delta <- Delta + crossprod(datacen[,id==i])    
      }
      if(voi!="columns") Sigma <- Sigma/(N-1)
      Delta <- Delta/(N-1)/p1
    }
    trDeltahat <- sum(diag(Delta)) 
    if(voi!="columns") Sigma <- Sigma/trDeltahat   
    if(centered) {
      if(shrink!="none"){
        trDelta2hat <- 0
        if(p1>=p2) { 
          y <- crossprod(datamat)
          for(i in seq_len(N-1))
          {
            helpmat <- y[id==i,id==i]
            sum1 <- 0 
            for(j in seq(i+1,N))
            {
              sum1 <- sum1 + y[id==j,id==j]
            }
            trDelta2hat <- sum(helpmat*sum1)+ trDelta2hat 
          }
        } else {
          id1 <- rep(seq_len(N),each=p1)
          datamat1 <- t(transposedata(datamat,N))
          y1 <- tcrossprod(datamat1)
          for(i in seq_len(N-1))  trDelta2hat <- sum(y1[id1==i,id1>i]^2)+ trDelta2hat 
        }
        trDelta2hat  <- 2*trDelta2hat/N/(N-1)/p1^2
        datavec <- matrix(datamat,p1*p2,N)
        helpcon <- 0
        for(i in seq_len(N-1)) helpcon <- sum(crossprod(datavec[,i],datavec[,seq(i+1,N)])^2) + helpcon 
        trOmega2hat <- 2*helpcon/N/(N-1)
        trSigma2hat <- trOmega2hat/trDelta2hat 
        if(shrink=="both" | shrink=="columns"){
          lambdaD_a <- trSigma2hat/(p1^2)/N*(trDeltahat^2+trDelta2hat)
          lambdaD_b <- trDelta2hat-trDeltahat^2/p2
          lambdaD <- lambdaD_a/(lambdaD_a+lambdaD_b)
          lambdaD <- max(min(1,lambdaD),0)} 
        if(shrink=="both" | shrink=="rows"){
          lambdaS_a <- trDelta2hat/trDeltahat^2/N*(p1^2+trSigma2hat)
          lambdaS_b <- trSigma2hat-p1
          lambdaS <- lambdaS_a/(lambdaS_a+lambdaS_b)
          lambdaS <- max(min(1,lambdaS),0)
        }}}    
    else {
      dataveccen <- matrix(datacen,p1*p2,N)
      if(shrink!="none"){
        Y2 <- Y4 <- Y5 <- 0  
        if(p1>=p2) { 
          y <- crossprod(datamat)
          for(i in seq_len(N-1))
          {
            helpmat <- y[id==i,id==i]
            sum1 <- 0 
            for(j in seq(i+1,N))
            {
              sum1 <- sum1 + y[id==j,id==j]
            }
            Y2 <- sum(helpmat*sum1)+ Y2
          }
          Y2 <- 2*Y2
          sum41 <- sum42 <- sum43 <- 0
          sum51 <- sum52 <- sum53 <- sum54 <- sum55 <- 0
          sum541 <- sum551 <- 0
          for(i in seq_len(N)){
            helpmat <- datamat[,id==i]
            helpmat2 <- y[id==i,id==i]
            helpmat3 <- datacen[,id==i]
            helpmat4 <- crossprod(helpmat,helpmat3)
            helpmat5 <- crossprod(helpmat3)
            sum41 <- sum41 + sum(c(helpmat5*c(helpmat2)))
            sum42 <- sum42 + sum(helpmat2^2)
            sum43 <- sum43 + sum(c(helpmat2)*y[id==i,id!=i])
            sum51 <- sum51 + sum(helpmat5^2)
            sum541 <- sum541 + sum(helpmat4^2)
            sum551 <- sum551 + sum(helpmat4*t(helpmat4))
            for(j in seq(i+1,N)){
              if(i!=N){
                sum52 <- sum52 + sum(y[id==i,id==j]^2)
                sum53 <- sum53 + sum(y[id==i,id==j]*y[id==j,id==i])        
              }     
            }   
          }
          Y4 <- (N^2)*sum41-((N-1)^2)*sum42-Y2+2*(N-1)*sum43
          sum52 <- 2*sum52
          sum53 <- 2*sum53
          sum54 <- (N^2)*sum541+2*(N-1)*sum43-((N-1)^2)*sum42-sum52
          sum55 <- (N^2)*sum551-((N-1)^2)*sum42+2*(N-1)*sum43-sum53
          Y5 <- ((N-1)*(N^2-3*N+3)*sum42+(2*N-3)*(Y2+sum52+sum53)+2*(N-3)*(Y4+sum54+sum55)-4*(N^2-3*N+3)*sum43-sum51*(N^3))/3     
        } else {
          id1 <- rep(seq_len(N),each=p1)
          datamat1 <- t(transposedata(datamat,N))
          y1 <- tcrossprod(datamat1)
          for(i in seq_len(N-1))  Y2 <- sum(y1[id1==i,id1>i]^2)+Y2
          Y2 <- 2*Y2
          sum41 <- sum42 <- sum43 <- 0
          sum51 <- sum52 <- sum53 <- sum54 <- sum55 <- 0
          sum541 <- sum551 <- 0
          for(i in seq_len(N)){
            helpmat <- datamat1[id1==i,]
            helpmat2 <- y1[id1==i,id1==i]
            helpmat3 <- datacen[,id==i]
            helpmat4 <- tcrossprod(helpmat3,helpmat)
            helpmat5 <- tcrossprod(helpmat3) 
            helpmat6 <- tcrossprod(helpmat)
            sum41 <- sum41 + sum(helpmat4^2)
            sum42 <- sum42 + sum(helpmat2^2)
            sum43 <- sum43 + sum(c(helpmat2)*y1[id1==i,id1!=i])
            sum51 <- sum51 + sum(helpmat5^2)
            sum541 <- sum541 + sum(helpmat6*helpmat5)
            sum551 <- sum551 + sum(helpmat4*t(helpmat4))
            for(j in seq(i+1,N)){
              if(i!=N){
                sum52 <- sum52 + sum(y1[id1==i,id1==i]*y1[id1==j,id1==j])
                sum53 <- sum53 + sum(y1[id1==i,id1==j]*y1[id1==j,id1==i])
              }     
            }   
          }
          Y4 <- (N^2)*sum41-((N-1)^2)*sum42-Y2+2*(N-1)*sum43
          sum52 <- 2*sum52
          sum53 <- 2*sum53
          sum54 <- (N^2)*sum541+2*(N-1)*sum43-((N-1)^2)*sum42-sum52
          sum55 <- (N^2)*sum551+2*(N-1)*sum43-((N-1)^2)*sum42-sum53
          Y5 <- ((N-1)*(N^2-3*N+3)*sum42+(2*N-3)*(Y2+sum52+sum53)+2*(N-3)*(Y4+sum54+sum55)-4*(N^2-3*N+3)*sum43-sum51*(N^3))/3               
        }
        Y2 <- Y2/N/(N-1)
        Y4 <- Y4/N/(N-1)/(N-2)
        Y5 <- Y5/N/(N-1)/(N-2)/(N-3)
        trDelta2hat <- (Y2-2*Y4+Y5)/(p1^2)
        helpcon <- 0
        for(i in seq_len(N-1)){
          helpcon <- sum(crossprod(dataveccen[,i],dataveccen[,seq(i+1,N)])^2) + helpcon 
        }
        colssums <- colSums(dataveccen^2)
        Q <- sum(colssums^2)
        tr2S <- (sum(colssums)/(N-1))^2
        trS2 <- (helpcon*2+Q)/(N-1)^2
        Q <- Q/(N-1)
        trOmega <- (N-1)/(N*(N-2)*(N-3))*((N-1)*(N-2)*trS2+tr2S-N*Q) 
        trSigma2hat  <- trOmega/trDelta2hat
        if(shrink=="both" | shrink=="columns"){
          lambdaD_a <- trSigma2hat/(p1^2)/(N-1)*(trDeltahat^2+trDelta2hat)
          lambdaD_b <- trDelta2hat-trDeltahat^2/p2
          lambdaD <- lambdaD_a/(lambdaD_a+lambdaD_b)
          lambdaD <- max(0,min(1,lambdaD))
        }  
        if(shrink=="both" | shrink=="rows"){
          lambdaS_a <- trDelta2hat/trDeltahat^2/(N-1)*(p1^2+trSigma2hat)
          lambdaS_b <- trSigma2hat-p1
          lambdaS <- lambdaS_a/(lambdaS_a+lambdaS_b)
          lambdaS <- max(min(1,lambdaS),0)
        }
      }
    }
    if(shrink=="both" | shrink=="rows"){ 
      if(voi!="columns") Sigma <- Sigma*(1-lambdaS) + diag(lambdaS,p1) 
    }
    if(shrink=="both" | shrink=="columns"){
      if(voi!="rows")  Delta <- Delta*(1-lambdaD) + diag(lambdaD*trDeltahat/p2,p2)
    }
    if(voi=="rows") ans <- list(rowcovmat=Sigma,lambdaS=if(shrink=="both" | shrink=="rows") lambdaS else NULL)
    if(voi=="columns") ans <- list(colcovmat=Delta,lambdaD=if(shrink=="both" | shrink=="columns") lambdaD else NULL)
    if(voi=="both") ans <- list(rowcovmat=Sigma,lambdaS=if(shrink=="both" | shrink=="rows") lambdaS else NULL,colcovmat=Delta,lambdaD=if(shrink=="both" | shrink=="columns") lambdaD else NULL)
    ans
  }