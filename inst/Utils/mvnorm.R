#Copyright (C) 2004-2019  Bernard Boulerice 

#require(moc)

# Multivariate normal for a sample matrix of observations x
# with parameters means=mu and Cholesky decomposition of the
# inverse of the Variances-Covariances matrix

mvnorm.invChol <- function(x,mu,shape=1,invChol)
{
  xdim <- dim(x)
  y <- (x-mu)/shape
  ss <- diag(xdim[2])
  lind <- upper.tri(ss,diag=TRUE)
  de <- cumsum(1:xdim[2])
  for ( i in 1:xdim[1]) {
    ss[lind] <- invChol[i,]
    nona <- !is.na(y[i,])
    y[i,nona] <- ss[nona,nona]%*%y[i,nona]
  }
  apply(dnorm(y)/shape*abs(invChol[,de]),1,prod,na.rm=TRUE)
}


# utility functions to convert Covariance matrix to Cholesky
# decomposition of its inverse and the converse function

Cov2invChol <- function(Sigma)
  {
    U <- chol(solve(Sigma))
    rbind(c(U[upper.tri(U,diag=TRUE)]))
  }

invChol2Cov <- function(invChol)
  {
    p <- (sqrt(length(invChol)*8+1)-1)/2
    tmp <- diag(p)
    tmp[upper.tri(tmp,diag=TRUE)] <- invChol
    solve(t(tmp)%*%tmp)[upper.tri(diag(p),diag=TRUE)]
  }


# Standard parametrization of Multivariate normal using mvnorm.invChol
# Sigma should be the vector form of the variance-covariance matrix 

mvnorm <- function(x,mu,shape=1,Sigma=NULL)
  {
    invChol <- t(apply(Sigma,1,
               function(s) Cov2invChol(matrix(s,p <- sqrt(length(s)),p))))
    mvnorm.invChol(x,mu,shape,invChol)
  }


# Multivariate normal likelihood computed from the sample means and covariances
# with parameters means and Cholesky decomposition of the
# inverse of the Variance-Covariance matrix.
# The maximum occurs when  muCov and muCholSigma corresponds

mvnorm.mucov <- function(muCov,muCholSigma,shape=1,extra=1)
  {ldim <- dim(muCov)
   p <- (sqrt(9+8*ldim[2])-3)/2
   if(round(p)-p != 0) stop("Ouch! Impossible muCov length!")
   tmp <- cbind(rep(NA,ldim[1]))
   for(i in 1:ldim[1]) {
     nona <- which(!is.na(c(muCov[i,(1:p)])))
     xbarmu <- rbind(muCov[i,(1:p)])-rbind(muCholSigma[i,(1:p)])
     xbarmu <- t(xbarmu)%*%xbarmu
     S <- diag(p)
     S[upper.tri(S,diag=TRUE)] <-  muCov[i,-(1:p)]
     S[lower.tri(S)] <- S[upper.tri(S)]
     invSigma <- diag(p)
     invSigma[upper.tri(invSigma,diag=TRUE)] <- muCholSigma[i,-(1:p)]
     invSigma <- t(invSigma)%*%invSigma
     tmp[i,] <- (2*pi)^(-p/2)*det(cbind(invSigma[nona,nona]))^(1/2)*
       exp(-sum(diag(invSigma[nona,nona]%*%(S[nona,nona]+xbarmu[nona,nona])))/2)}
   tmp
 }
