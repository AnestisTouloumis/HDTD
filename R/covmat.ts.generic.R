covmat.ts.generic <-
  function(datamat,N,centered){ 
    p1 <- nrow(datamat)
    p2 <- ncol(datamat)/N
    id <- rep(seq_len(N),each=p2) 
    Y1 <- sum(datamat^2)/N
    Y2 <- Y3 <- Y4 <- Y5 <- 0  
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
      if(centered==FALSE){
        xmean <- matrix(rowMeans(matrix(datamat,p1*p2,N)),p1,p2)
        sum41 <- sum42 <- sum43 <- 0
        sum51 <- sum52 <- sum53 <- sum54 <- sum55 <- 0
        sum541 <- sum551 <- 0
        for(i in seq_len(N)){
          helpmat <- datamat[,id==i]
          helpmat2 <- y[id==i,id==i]
          helpmat3 <- helpmat-xmean
          helpmat4 <- crossprod(helpmat,helpmat3)
          helpmat5 <- crossprod(helpmat3)
          Y3 <- sum(c(helpmat)*datamat[,id>i])+Y3
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
      }  
    }
    else {
      id1 <- rep(seq_len(N),each=p1)
      datamat1 <- t(transposedata(datamat,N))
      y1 <- tcrossprod(datamat1)
      for(i in seq_len(N-1))  Y2 <- sum(y1[id1==i,id1>i]^2)+Y2
      Y2 <- 2*Y2
      if(centered==FALSE){
        xmean <- matrix(rowMeans(matrix(datamat,p1*p2,N)),p1,p2)
        sum41 <- sum42 <- sum43 <- 0
        sum51 <- sum52 <- sum53 <- sum54 <- sum55 <- 0
        sum541 <- sum551 <- 0
        for(i in seq_len(N)){
          Y3 <- sum(c(datamat[,id==i])*datamat[,id>i])+Y3
          helpmat <- datamat1[id1==i,]
          helpmat2 <- y1[id1==i,id1==i]
          helpmat3 <- helpmat-xmean
          helpmat4 <- tcrossprod(helpmat3,helpmat)
          helpmat5 <- tcrossprod(helpmat3) 
          helpmat6 <- tcrossprod(helpmat)
          sum41 <- sum41 + sum(helpmat4^2)
          sum42 <- sum42 + sum(helpmat2^2)
          sum43 <- sum43 + sum(c(helpmat2)*y1[id1==i,id1!=i])
          sum51 <- sum51 + sum(helpmat5^2)
          sum541 <- sum541 + sum(helpmat6*helpmat5)
          sum551 <- sum551 + sum(helpmat4*t(helpmat4))
          if(i!=N){
          for(j in seq(i+1,N)){
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
    }
    Y2 <- Y2/N/(N-1)    
    Y3 <- 2*Y3/(N-1)/N
    Y4 <- Y4/N/(N-1)/(N-2)
    Y5 <- Y5/N/(N-1)/(N-2)/(N-3)
    T1 <- (Y1-Y3)/p1
    T2 <- (Y2-2*Y4+Y5)/(p1^2)
    dataveccen <- matrix(datamat,p1*p2,N)
    if(!centered){
      dataveccen <- dataveccen - rowMeans(dataveccen)
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
    } else {
      helpcon <- 0
      for(i in seq_len(N-1)){
        helpcon <- sum(crossprod(dataveccen[,i],dataveccen[,seq(i+1,N)])^2) + helpcon 
      }
      trOmega <- 2*helpcon/N/(N-1)
    }
    con <- trOmega/T2/(p1^2)
    Utest <- N/2/con*(p2*T2/(T1^2)-1)
    Vtest <- N/2/con*(T2/p2-2*T1/p2+1)
    ans <- list(Utest=Utest,Vtest=Vtest)
    ans
  }
