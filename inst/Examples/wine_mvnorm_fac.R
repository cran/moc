## This example shows how to build a mixture of factor analyzers in MOC.
## We use the wine data of Forina, M. et al, PARVUS
## available at UCI repository (http://www.ics.uci.edu/~mlearn/MLRepository.html)
## and in package=gclus or at http://www.maths.uq.edu.au/~gjm/DATA/mmdata.html

wine.file <- system.file("Examples","wine.RData",package="moc",mustWork=TRUE)
load(wine.file)

require(moc)    ## load the essential

## We need the multivariate normal distribution, assuming
## that the extra parameter contains the parameters corresponding
## to the inverse of the cholesky decomposition of the correlation matrix.

mvnorm <- function(x,mu,sigma,extra)
{ 
  y <- (x-mu)/sigma
  ss <- diag(rep(1,dim(x)[[2]]))
  lind <- upper.tri(ss)
  for ( i in 1:dim(x)[1]) {
      ss[lind] <- extra[i,]
      y[i,] <- t(ss)%*%y[i,]
  }
  apply(dnorm(y)/sigma,1,prod,na.rm=TRUE)
}


mvnorm.fac <- 
function(x,mu,sigma,extra)
{ 
  y <- (x-mu)/sigma
  ss <- diag(rep(1,dim(x)[[2]]))
  lind <- upper.tri(ss,diag=TRUE)   ## for the factorial structure we also need the diagonal elements
  de <- cumsum(1:dim(x)[[2]])
  for ( i in 1:dim(x)[1]) {
      ss[lind] <- extra[i,]
      y[i,] <- t(ss)%*%y[i,]
  }
  apply(dnorm(y)/sigma*extra[,de],1,prod,na.rm=TRUE)
}

mu.mvnorm <- function(p) {rbind(p)}  ## base mean function for mvnorm

sig.mvnorm <- function(p) {rbind(exp(p))}  ## base standard deviation function for mvnorm

extra.mvnorm <- function(p) {rbind(p)}  ## base extra parameter function for mvnorm

extra.mvnorm.fac <- function(p,fact) {r <- rbind(apply(matrix(p,ncol=fact),1,function(x) x/sqrt(sum(x^2)+1)))
                             Mr <- t(r)%*%r;                               ## base extra parameter function for mvnorm with
                             diag(Mr) <- 1                                 ## factorial correlation structure
                             rbind((solve(chol(Mr))[row(Mr)<=col(Mr)]))}                                                                              
## Note that mixture models are quite good to approximate multivariate
## distributions with a mixing of locally independent marginal distributions,
## but the number of groups to achieve this can be large especially
## when the number of variables is large as in this case.
## We then often prefer to drop the local independance constraint and
## model a multivariate structure in each group directly, even in this
## case the number of parameters becomes large ( for 13 variables the
## covariances matrix has 91 parameters ). So using a factorial
## structure help decrease the number of parameters.

## First we simply compute a single group with 1 factor structure:

ll1 <- unclass(loadings(factanal(wine[,2:14],factor=1,rotation="none")))

wine.1f.1g <- moc(wine[,2:14],density=mvnorm.fac,joint=TRUE,
              gmu=list(G1=function(p) mu.mvnorm(p[1:13])),
              gshape=list(G1=function(p) sig.mvnorm(p[1:13])),
              gextra=list(G1=function(p) extra.mvnorm.fac(p[1:13],fact=1)),
check.length=FALSE,pgmu=sapply(wine[,2:14],mean),pgshape=log(sapply(wine[,2:14],sd)),pgextra=ll1/sqrt(1-ll1^2))

## Note that the preceeding model already has 39 parameters, 13 for the correlation matrix
## with only 178 observations, but when we increase the number of groups to 3 we have
## already 39 parameters for the correlation matrix solely and 110 parameters total
## which is already quite large with respect to the sample size.
## Though the mixture groups need not correspond to any predefined (a priori) known groups 
## we first try to make a correspondance with known cultivars (wine[,1]=wine$Class).

ll2 <- c(unclass(loadings(factanal(wine[,2:14],factor=1,rotation="none",subset=wine[,1]==1))),  ##for starting values
         unclass(loadings(factanal(wine[,2:14],factor=1,rotation="none",subset=wine[,1]==2))),
         unclass(loadings(factanal(wine[,2:14],factor=1,rotation="none",subset=wine[,1]==3))))

## Be patient it can take a long time to estimate.

wine.1f.3g <- moc(wine[,2:14],density=mvnorm.fac,joint=TRUE,groups=3,
               gmu=list(G1=function(p) mu.mvnorm(p[1:13]),G2=function(p) mu.mvnorm(p[14:26]),G3=function(p) mu.mvnorm(p[27:39])),
               gshape=list(G1=function(p) sig.mvnorm(p[1:13]),G2=function(p) sig.mvnorm(p[14:26]),G2=function(p) sig.mvnorm(p[27:39])),
               gextra=list(G1=function(p) extra.mvnorm.fac(p[1:13],fact=1),G2=function(p) extra.mvnorm.fac(p[14:26],fact=1),
               G3=function(p) extra.mvnorm.fac(p[27:39],fact=1)) ,
               check.length=FALSE,
               pgmu=unlist(by(wine[,2:14],wine[,1],sapply,mean)),
               pgshape=log(unlist(by(wine[,2:14],wine[,1],sapply,sd))),
               pgextra= ll2/sqrt(1-ll2^2),pgmix=c(0.2,-0.2),gradtol=1e-3)


## We can then estimates the mixture probabilities for the known cultivar classes

class.mix <- function(p) {t(apply(cbind(wine[,1]),1, function(x) inv.glogit(c(p[1]+p[2]*x,p[3]+p[4]*x))))}

wine.1f.3g.class <- moc(wine[,2:14],density=mvnorm.fac,joint=TRUE,groups=3,
               gmu=list(G1=function(p) mu.mvnorm(p[1:13]),
                        G2=function(p) mu.mvnorm(p[14:26]),
                        G3=function(p) mu.mvnorm(p[27:39])),
               gshape=list(G1=function(p) sig.mvnorm(p[1:13]),
                           G2=function(p) sig.mvnorm(p[14:26]),
                           G3=function(p) sig.mvnorm(p[27:39])),
               gextra=list(G1=function(p) extra.mvnorm.fac(p[1:13],fact=1),
                           G2=function(p) extra.mvnorm.fac(p[14:26],fact=1),
                           G3=function(p) extra.mvnorm.fac(p[27:39],fact=1)) ,
               check.length=FALSE,gmixture=class.mix,
               pgmu=wine.1f.3g$coef[1:39],pgshape=wine.1f.3g$coef[40:78],
               pgextra=wine.1f.3g$coef[78:116],
               pgmix=c(wine.1f.3g$coef[117],0,wine.1f.3g$coef[118],0),gradtol=1e-3)

