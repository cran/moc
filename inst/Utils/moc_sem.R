# Copyright (C) 2006-2019  Bernard Boulerice 

require(moc)
source(system.file("Utils","mvnorm.R",package="moc"),local=TRUE)

# function for direct computation of the log-likelihood of the saturated
# model with mvnorm.mucov or mvnorm.invchol

Sat2logLik <- function(Cov,n,full=FALSE)
  { p <- dim(Cov)[1]
    cbind("-2*lokLik"=n*(p+p*log(2*pi)+log(det(Cov))),Df=ifelse(full,n-p*(p+3)/2,0))}


# Some utilitary functions

coef.moc <- function(object,split=FALSE,...)
  {
    if (!inherits(object,"moc")) stop("  Not a moc object!")
    if(split) return(split(object$coef,rep(c("mu","shape","extra","mixture"),object$npar)))
    else return(object$coef)
  }

Cov2stdChol <- function(Sigma)
  {
    U <- chol(Sigma)
    sdim <- dim(Sigma)
    if(sdim[1]!=sdim[2]) stop("Sigma should be a square matrix")
    omega <- diag(U)
    U <-solve(U%*%diag(1/omega))
    list(omega=omega,UChol=U[upper.tri(diag(sdim[1]))])
  }

stdChol2Cov <- function(omega,UChol)
  {
    p <- length(omega)
    tmp <- diag(rep(1,p))
    tmp[upper.tri(diag(p))] <- UChol
    tmp2 <- solve(tmp)
    diag(omega)%*%t(tmp2)%*%tmp2%*%diag(omega)
  }

ParCor <- function(cor) {log((1+cor)/(1-cor))}

             # dont return cor of 1 or -1 to avoid numerical instabilities
invParCor <- function(z) {0.999*(exp(z)-1)/(exp(z)+1)}


# Functions for constructing SEM in conjonction with moc

# The RAM path diagram approach to SEM specifies three matrices:

# 1 B  for factor loadings of theorized latent variables (factors) and observed variables
#      on each other. These correspond to singles arrows in path diagram.
#      This is thus a square matrix with nl+no rows and columns.
#      (no: number of observed variables, nl: number of latent variables)

# 2 S  is for symmetric paths or two headed arrows and specify correlations
#      and variances among the variables.
#      Thus it is a square symmetric matrix of dimension (nl+n0)x(nl+n0)

# 3 F  is for filtering the observed variables from the latent ones and is typically
#      fixed. It has dimension no x (nl+no).

# The resulting Variance-Covariance matrix for the observed variables is
#      F * (I-B)^{-1} * S * (I-B)^{-1}' F'    (' denotes matrix transpose)

RAM.muInvChol <- function(F=NULL,S=NULL,B=NULL,mu=NULL,db=dim(B))
  {
    IB <- solve(diag(dB[1])-B)
    cov <- F %*% IB %*% S %*% t(IB)%*% t(F)
    mu0 <- mu
    if(!is.null(mu)) mu0 <- mu%*%t(IB)%*%t(F)
    rbind(c(muO=mu0,invChol=Cov2invChol(cov)))
  }

# Notice that S is a covariance matrice, so parameters associated with it
# are not totally free. You should consider using the function invParCor for parametrization,
# as in the following construction helper functions

RAM.Fstd <- function(nl,no) {cbind(matrix(0,no,nl),diag(no))}
RAM.S <- function(logstdL=0,logstdO=0,nl=length(logstdL),no=length(logstdO),
                  zcorL=rep(0,nl*(nl-1)/2),zcorO=rep(0,no*(no-1)/2),zcorLO=rep(0,no*nl))
  {
    S <- diag(nl+no)
    S[lower.tri(S)] <- invParCor(c(zcorL,zcorLO,zcorO))
    S[upper.tri(S)] <- t(S)[upper.tri(S)]
    Std <- diag(exp(c(logstdL,logstdO)))
    Std%*%S%*%Std
  }
RAM.B <- function(LL=0,LO=0,OL=0,OO=0)
  { rbind(cbind(LL,LO),cbind(OL,OO))} # Note: LO are for {L <--- O} paths and OL for {O <--- L} paths.

# Simple basic verifications of RAM matrices properties

check.RAM <- function(F=NULL,S=NULL,B=NULL,mu=NULL)
  {
    if(!is.matrix(F)) stop("F should be a matrix.\n")
    if(!is.matrix(S)) stop("S should be a matrix.\n")
    if(!is.matrix(B)) stop("S should be a matrix.\n")
    dF=dim(F)
    dS=dim(S)
    dB=dim(B)
    if(diff(dS)!=0) stop("S should be a square-matrix.\n")
    if(diff(dB)!=0) stop("B should be a square-matrix.\n")
    if(any(dS!=dB))  stop("S and B should have the same dimensions.\n")
    if(dF[2]!=dB[2])stop("F should have the same number of columns as S and B.\n")
    varS <- diag(S)
    if(!all(varS > 0)) stop("Diagonal elements of S should be > 0\n")
    corS <- diag(sqrt(1/varS)) %*% S %*% diag(sqrt(1/varS))
    if(!all(corS <=1 & corS >= -1))
      stop("S should define a covariance matrix such that -1 < cor < 1 \n")

    if(!is.null(mu) && length(mu)!=db[2]) stop("Wrong mu length,\n")
    c(nl=dF[2]-dF[1],no=dF[1])
  }

# LOB specification of SEM is much simpler than RAM and does not require inversion
# of a matrix in its computions. It correspond more closely to standard factor analysis
# model X = B Z + E with Cov(Z) = L and Cov(E) = O with Cov(Z,E) = {0}
# It also requires three matrices:

# 1 B  m by n matrix of factor loadings of latent variables on observed variables

# 2 L  n by n symmetric variances-covariances matrix of the [L]atent variables (factors)

# 3 O  m by m symmetric variances-covariances matrix of the residuals for the [O]bserved variables.

# So, the resulting covariances matrix for the observed variables is:
#      B * L *B' + O 
# notice that O and L are covariances matrices, so parameters associated with them
# are not totally free. You should consider using the function invParCor for parametrization.

LOB.muInvChol <- function(L=NULL,O=NULL,B=NULL,mu=NULL)
  {
    dB <- dim(B)
    cov <- B %*% L %*% t(B) + O
    mu0 <- mu
    if(!is.null(mu)) mu0 <-  mu[,1:dB[2]] %*% t(B) +  mu[,(dB[2]+1):sum(dB)]
    rbind(c(muO=mu0,invChol=Cov2invChol(cov)))
  }

# Simple constructors for LOB model matrices

LOB.B <- function(OL,nl=1,no=1) { matrix(OL,no,nl)} # just organize the loadings parameters

LOB.COV <- function(
                    =0,n=length(logstd),  # simply construct a covariances matrix from
                    zcor=rep(0,n*(n-1)/2))      # log(std.dev) parameters and ParCor(cor) parameters.
                                                # Used for L and O matrices.
  {
    S <- diag(n)
    S[lower.tri(S)] <- invParCor(zcor)
    S[upper.tri(S)] <- S[lower.tri(S)]
    Std <- diag(exp(c(logstd)))
    Std%*%S%*%Std
  }


# Simple basic verifications of LOB matrices properties

check.LOB <- function(L=NULL,O=NULL,B=NULL,mu=NULL)
  {
    if(!is.matrix(O)) stop("O should be a matrix.\n")
    if(!is.matrix(L)) stop("L should be a matrix.\n")
    if(!is.matrix(B)) stop("S should be a matrix.\n")
    dO=dim(O)
    dL=dim(L)
    dB=dim(B)
    if(diff(dL)!=0) stop("L should be a square-matrix.\n")
    varL <- diag(L)
    if(!all(varL > 0)) stop("Diagonal elements of L should be > 0\n")
    corL <- diag(sqrt(1/varL)) %*% L %*% diag(sqrt(1/varL))
    if(!all(corL <=1 & corL >= -1))
      stop("L should define a covariance matrix such that -1 < cor < 1 \n")
    if(diff(dO)!=0) stop("O should be a square-matrix.\n")
    varO <- diag(O)
    if(!all(varO > 0)) stop("Diagonal elements of O should be > 0\n")
    corO <- diag(sqrt(1/varO)) %*% O %*% diag(sqrt(1/varO))
    if(!all(corO <=1 & corO >= -1))
      stop("O should define a covariance matrix such that -1 < cor < 1 \n")
    if(any(dB!=c(dO[1],dL[2])))
      stop("B should have the same number of rows as O\n\t and the same number of columns as L.\n")
    if(!is.null(mu) && length(mu)!=sum(dB)) stop("Wrong mu length,\n")
    c(nl=dL[1],no=dO[1])
  }

