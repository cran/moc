\name{moc}
\alias{moc}
\title{Fit a General Nonlinear Multivariate Mixture Model}

\description{
    \code{moc} fits user-specified mixture models with one,two
    and three parameters distributions to multivariate data that can be of
    discrete or continuous type and a mix of both types. The likelihood for the
    vector of observations or repeated measurements for
    subject \eqn{i} has the form

    \deqn{f( Y_i = y_i |\, Z_i = z_i, X_i = x_i)  = \sum_k P( G_i = k |\, z_i)
    \,h( y_i |\, G_i = k, x_i)}{f( Y[i] = y[i] | z[i], x[i]) =
    \Sum_k P( G[i] = k | z[i]) h( y[i] | G[i] = k, x[i])}
 
    The \eqn{G_i}{G[i]} represent the mixture groups and \eqn{h()}
    the conditional joint density of the observations \eqn{Y_i}{Y[i]}
    given the covariates \eqn{(Z_i,X_i)}{(Z[i],X[i])}.
    The user supplies either the joint or marginal conditional
    density(ies) of the components of \eqn{Y_i}{Y[i]}. In the latter case,
    the joint conditional density is constructed by taking the product of
    the marginal densities ( assuming conditional independence ).
  }
  
\usage{

moc(y, density=NULL, joint=FALSE, groups=1,
gmu=NULL, gshape=NULL, gextra=NULL, gmixture=inv.glogit, expected = NULL,
pgmu=NULL, pgshape=NULL, pgextra=NULL, pgmix=NULL,
check.length=TRUE, scale.weight=FALSE, wt=1, data=NULL,
ndigit=10, gradtol=0.0001, steptol=gradtol,iterlim=100,print.level=2,\dots)
}

\arguments{
  \item{y}{ A matrix or data frame giving the vectors of observations
  (of length \eqn{nvar}) for each subject.}

  \item{density}{ A function returning the conditional joint or marginal
  density of the observations and calling the location, shape and extra functions.}

  \item{joint}{ Specify if the density gives the joint or common marginal
    density of the vector of observations. When using a joint density
    remenber that this density will receive its parameters as matrices.}

  \item{groups}{ Number of mixtures.}

  \item{gmu, gshape, gextra}{ User-specified \code{lists} of functions
  returning the location, shape and extra density parameters for each
  mixture group and observation as a function of the parameters
  \emph{pgmu}, \emph{pgshape}, \emph{pgextra} and covariates.}

  \item{gmixture}{ A user-specified function of \emph{pgmix}, giving the regression 
  function of the mixture probabilities. The default is the
  inverse generalized logit with respect to the first group.}

  \item{expected}{ A \code{list} of functions returning the expected response value with
  respect to all the parameters \code{c(pgmu,pgshape,pgextra,pgmix)} for
  each mixture groups. Defaults to \code{gmu}.}

  \item{pgmu, pgshape, pgextra, pgmix}{ Vector of initial estimates for the
  parameters of the location, shape, extra and mixture functions.}


  \item{wt}{ Vector of subjects sampling weights. Actually the program uses
  standard sample-weighted log-Likelihood assuming fixed weights. }

  \item{scale.weight}{ Logical value specifying if the vector of weights 
  \emph{wt} should be rescaled to sum to the sample size.}

\item{check.length}{ Logical value specifying weither to check the rows
  length returned by the functions in \code{gextra} against the number
  of variables in \eqn{y}. Especially usefull when the \code{density}
  requires more paramaters then the number of variables like covariance
  parameters for multivariate normal.}

  \item{data}{ An optional \emph{data frame} or \emph{list} containing
  some or all variables and functions required to fit the model. It is
  \code{attach}ed at run time and subsequently by \code{moc}
  \code{methods} applied to the fitted model. Be aware that
  \code{objects} in \emph{.GlobalEnv} will mask any \code{object} with
  the same name in \emph{data}.
  }
  
  \item{ndigit, gradtol, steptol, iterlim, print.level, \ldots}{
    Arguments controlling \code{nlm}.}
}

\details{
    The procedure minimizes the resulting \eqn{-\log(Likelihood)} without
    constraints, the parameters are all assumed to be real numbers.
    Thus the user should supply appropriate link functions and parameterize
    the density and parameters functions accordingly ( see the examples
    ). By default missing values in the response variables \eqn{y} are assumed
    to be missing at random, that is the \eqn{Likelihood} for the subset
    of valid observations is just the marginal \eqn{Likelihood} for this
    subset in each mixture. Specific treatment of missing values in the
    response variables can be achieved by handling them explicitly in
    the functions \code{density},\code{gmixture},\code{gmu},\code{gshape}
    and \code{gextra}. The function \code{density} can return \code{NA}
    and yields the default treatment of missing values in the response.
    The functions
    \code{gmixture},\code{gmu},\code{gshape} and \code{gextra},
    can not return \code{NA} thus missing values in the covariates
    should be treated explicitly by these functions.

    The \code{lists} of functions \code{gmu}, \code{gshape},
    \code{gextra} returns the location, shape and extra parameters
    to the density for each observation and mixture group as a function
    of \emph{pgmu}, \emph{pgshape} and \emph{pgextra} and covariates.
    Each function should return a vector of length \eqn{nvar} or a
    matrix of such vectors. The first function in the \code{list} is
    for the first group, the second function for the second group and so
    on. The functions in the same \code{list} share the same parameters
    but the different \code{lists} have different parameters
    ( see the examples ). 

    Setting the \code{attributes} \emph{parameters} for functions
    \code{gmu},\code{gshape},\code{gextra} and \code{gmixture} will
    generate parameter labels in the printout of the object.
    
    The deviance residuals, fitted values and posterior probabilities
    are obtained through the use of the methods \code{residuals},\code{fitted}
    and \code{post}.
}


\value{
    A list of class \emph{moc} is returned that contains all of the
    relevant information calculated, including error code generated by
    \code{nlm}.
    The printed output includes  \eqn{-2 \log(Likelihood)},
    the corresponding df , AIC, BIC, entropy and ICL-BIC
    ( entropy corrected BIC ), mean mixture probabilities,
    mean expected and observed values for each mixture group,
    the maximum likelihood estimates and standard errors.}

\note{
    The \code{expected} function is used to compute the fitted values
    and response residuals (not deviance). It is especially useful when
    the expected value differs from the location parameters as for censored
    normal or zero inflated Poisson distributions.
  
    The method of fixed sample-weight provides design-consistent parameters
    estimates. However, for the moment the program does not provide any
    methods to include sampling variances resulting from weights
    estimation. If the user wants to incorporate weights estimation
    sampling variances it could be achieved, for example, by including
    \code{moc} model estimation in a \emph{jackknife} loop.

    Be aware that degrees of freedom (df) for mixture models are usually
    useless ( if not meaningless ) and likelihood-ratio of
    \emph{apparently} nested models often does not converge to a
    Chi-Square with corresponding df.
  }   


  \references{
McLachlan, Geoffrey and Peel, David (2000) \emph{Finite mixture models},Wiley-Interscience, New York.

Roeder, K., Lynch, K. and Nagin, D. (1999) \emph{Modeling Uncertainty in
  Latent Class Membership: A Case Study in Criminology}, J. Amer. Statist. Assoc., \bold{94}, pp. 766--776.

Lindsay, Bruce G. and Roeder, K. (1992) \emph{Residual diagnostics for mixture models}, J. Amer. Statist. Assoc., \bold{87},
     pp. 785--794.
   }
   
\author{Bernard Boulerice <\email{Bernard.Boulerice@umontreal.ca}>}

\seealso{\code{\link{print.moc}},\code{\link{plot.moc}},\code{\link{residuals.moc}},
  \code{\link{plot.residuals.moc}},\code{\link{fitted.moc}},\code{\link{post.moc}},\code{\link{AIC.moc}},
    \code{\link{logLik.moc}},\code{\link{obsfit.moc}},\code{\link{nlm}}}

\examples{

data(moc.dat)

cnorm.dat<-list()   #This is used as a container for functions and data

# Censored Normal (marginal density)

cnorm.dat$cnorm<-function(x,mu,sig,min,max)
{mi<-(x==min)*1
ma<-(x==max)*1
mi*pnorm((min-mu)/sig)+ma*(1-pnorm((max-mu)/sig))+
(1-mi-ma)*dnorm((x-mu)/sig)/sig}

# For this data set the range of the dependent variables is [0,14]

cnorm.dat$cnorm1<-function(x,mu,sig,...) {cnorm(x,mu,sig,0,14)}

# We have 4 observations

cnorm.dat$gmu1<- list(
   Group1=function(pmu) {t(1)\%*\%rep(pmu[1],4)},
   Group2=function(pmu) {t(1)\%*\%rep(pmu[2],4)},
   Group3=function(pmu) {t(1)\%*\%rep(pmu[3],4)})

attr(cnorm.dat$gmu1,"parameters")<-c("  cons1","  cons2","  cons3")

# Expected value of a general censored normal

cnorm.dat$cmean<-function(mu,sig,min,max) {
max-(max-mu)*pnorm((max-mu)/sig)+(min-mu)*pnorm((min-mu)/sig)-
sig*(dnorm((max-mu)/sig)-dnorm((min-mu)/sig)) }

# Homogeneous variances

cnorm.dat$gshape1<- list(
   Group1=function(psh) {t(1)\%*\%rep(exp(psh[1]),4)},
   Group2=function(psh) {t(1)\%*\%rep(exp(psh[1]),4)},
   Group3=function(psh) {t(1)\%*\%rep(exp(psh[1]),4)})

attr(cnorm.dat$gshape1,"parameters")<-c("  log(std.dev)")

cnorm.dat$cmean1<- list(
   Group1=function(p) {cmean(gmu1[[1]](p[1:3]),gshape1[[1]](p[4]),0,14) },
   Group2=function(p) {cmean(gmu1[[2]](p[1:3]),gshape1[[2]](p[4]),0,14) },
   Group3=function(p) {cmean(gmu1[[3]](p[1:3]),gshape1[[3]](p[4]),0,14) })

moc1<-
moc(moc.dat[,1:4],density=cnorm1,groups=3,gmu=gmu1,gshape=gshape1,expected=cmean1,
pgmu=c(0.5, 2, 5),pgshape=c(0.7),pgmix=c(-0.6, -2.0), data=cnorm.dat,gradtol=1E-4)

\dontrun{
# Heterogeneous variances across mixture groups

cnorm.dat$gshape2<-list(
   Group1=function(psh) {t(1)\%*\%rep(exp(psh[1]),4)},
   Group2=function(psh) {t(1)\%*\%rep(exp(psh[2]),4)},
   Group3=function(psh) {t(1)\%*\%rep(exp(psh[3]),4)})

cnorm.dat$cmean2<-list(
   Group1=function(p) {cmean(gmu1[[1]](p[1:3]),gshape2[[1]](p[4:6]),0,14) },
   Group2=function(p) {cmean(gmu1[[2]](p[1:3]),gshape2[[2]](p[4:6]),0,14) },
   Group3=function(p) {cmean(gmu1[[3]](p[1:3]),gshape2[[3]](p[4:6]),0,14) })

moc2<-
moc(moc.dat[,1:4],density=cnorm1,groups=3,gmu=gmu1,gshape=gshape2,expected=cmean2,
pgmu=moc1$coef[1:3],pgshape=c(rep(moc1$coef[4],3)),pgmix=moc1$coef[5:6]
,data=cnorm.dat,gradtol=1E-4)
}
# Heterogeneous variances across time

cnorm.dat$gshape3<-list(
   Group1=function(psh) {exp(t(1)\%*\%psh[1:4])},
   Group2=function(psh) {exp(t(1)\%*\%psh[1:4])},
   Group3=function(psh) {exp(t(1)\%*\%psh[1:4])})

cnorm.dat$cmean3<-list(
   Group1=function(p) {cmean(gmu1[[1]](p[1:3]),gshape3[[1]](p[4:7]),0,14)},
   Group2=function(p) {cmean(gmu1[[2]](p[1:3]),gshape3[[2]](p[4:7]),0,14)},
   Group3=function(p) {cmean(gmu1[[3]](p[1:3]),gshape3[[3]](p[4:7]),0,14)})

moc3<-
moc(moc.dat[,1:4],density=cnorm1,groups=3,gmu=gmu1,gshape=gshape3,expected=cmean3,
pgmu=moc1$coef[1:3],pgshape=c(rep(moc1$coef[4],4)),pgmix=moc1$coef[5:6],
data=cnorm.dat,gradtol=1E-4)

cnorm.dat$ages<-cbind(1.7,3,4.2,5.6)

\dontrun{
# Last group is a linear function of time

cnorm.dat$gmu1t<-list(
   Group1=function(pmu) {pmu[1]*ages^0},
   Group2=function(pmu) {pmu[2]+pmu[3]*ages},
   Group3=function(pmu) {pmu[4]*ages^0})

cnorm.dat$cmean1t<-list(
   Group1=function(p) {cmean(gmu1t[[1]](p[1:4]),gshape1[[1]](p[5]),0,14)},
   Group2=function(p) {cmean(gmu1t[[2]](p[1:4]),gshape1[[2]](p[5]),0,14)},
   Group3=function(p) {cmean(gmu1t[[3]](p[1:4]),gshape1[[3]](p[5]),0,14)})

moc4<-
moc(moc.dat[,1:4],density=cnorm1,groups=3,gmu=gmu1t,gshape=gshape1,expected=cmean1t,
pgmu=append(moc1$coef[1:3],0,after=2),pgshape=c(moc1$coef[4]),pgmix=moc1$coef[5:6],
data=cnorm.dat,gradtol=1E-4)


# Zero inflated Poisson log-linear in time for the third group
# Be carefull dpois requires integer values

zip<- function(x,la,shape=1,extra)
{ mix<- exp(extra)/(1+exp(extra))
  mix*(x==0)+(1-mix)*dpois(x,la) }
}

ages<-cbind(1.7,3,4.2,5.6)

gmup<-list(
   Group1=function(pmu) {exp(pmu[1]*ages^0)},
   Group2=function(pmu) {exp(pmu[2]+pmu[3]*ages)},
   Group3=function(pmu) {exp(pmu[4]*ages^0)})

\dontrun{
zipfit<-list(
   Group1=function(p) { gmup[[1]](p)/(1+exp(p[5]))},
   Group2=function(p) { gmup[[2]](p)/(1+exp(p[5]))},
   Group3=function(p) { gmup[[3]](p)/(1+exp(p[5]))})

gextrap<-list(
   Group1=function(pxt) {t(1)\%*\%rep(pxt[1],4)},
   Group2=function(pxt) {t(1)\%*\%rep(pxt[1],4)},
   Group3=function(pxt) {t(1)\%*\%rep(pxt[1],4)})

moc5<-
moc(moc.dat[,1:4],density=zip,groups=3,gmu=gmup,gextra=gextrap,expected = zipfit,
pgmu=c(-0.6, 0.64,0, 1.6),pgextra=c(-3),pgmix=c(-0.7, -2), gradtol=1E-4)
}

# Standard Poisson with mixture depending on time independent
# dichotomous covariate
# Be aware that dpoiss require integer values

dumm<-moc.dat[,5]-1
gmixt<-function(pm){
mix<-cbind(1,dumm)\%*\%matrix(pm[1:4],2,2)
cbind(1,exp(mix))/(1+apply(exp(mix),1,sum))}

poiss<-function(x,la,...) {dpois(x,la)}

moc6<-
moc(moc.dat[,1:4],density=poiss,groups=3,gmu=gmup,gmixture=gmixt,
pgmu=c(-0.6,0.64, 0, 1.6),pgmix=c(-0.2,-1, -1 ,-2),gradtol=1E-4)

obsfit.moc(moc6,along=dumm)
\dontrun{
 plot(moc6,against=ages,main="MOC profiles",xlab="age",ylab="Y")
 plot(residuals(moc6))
}
}
\keyword{models}
\keyword{nonlinear}
\keyword{cluster}
\keyword{classif}
\keyword{multivariate}
