\name{moc}
\alias{moc}
\title{General Nonlinear Mixture of Curves.}

\description{
  \code{moc} fits user-specified mixture of curves models to one,two
  and three parameters distributions. The likelihood for the
  vector of observation over time or repeated measurements of
  a subject takes the form

  \deqn{f( Y_i = y_i | Z_i = z_i, W_i = w_i)  = \sum_k P( C_i = k | z_i) h( y_i | C_i = k, w_i)}{f( Y[i] = y[i] | Z[i] = z[i], W[i] = w[i])  = Sum_k P( C[i] = k | z[i]) h( y[i] | C[i] = k, w[i])}

 
  The \eqn{C_i}{C[i]} represents the mixture groups ( or latent classes ) and \eqn{h()}
  the conditional joint density of \eqn{Y_i}{Y[i]} given the covariates.
  The user supplies either the joint or marginal conditional
  density(ies) of the components of \eqn{Y_i}{Y[i]}. In the latter case,
  the joint conditional density is constructed by taking the product of
  the marginal densities. Censoring should be  handled by the supplied density. 
  The user can also define the mixture probability function.

  The procedure minimizes the resulting \eqn{-\log likelihood} without
  constraints, the parameters are all assumed to be real numbers.
  Thus the user should supply appropriate link functions and parameterize
  the density accordingly ( see the examples ).

  The printed output includes  \eqn{-2 \log likelihood},
  the corresponding df, AIC and BIC, mean mixture probabilities,
  mean expected and observed latent curves,
  the maximum likelihood estimates and standard errors.

  The deviance residuals, fitted values and posterior probabilities
  are obtained through the use of the methods \code{residuals},\code{fitted}
  and \code{post}. The user has the option of weighting the residuals
  according to the posterior mixture probabilities.
}
\usage{

       moc(y, density=NULL, joint=FALSE, groups=1,
       gmu=NULL, gshape=NULL, gextra=NULL, gmixture=glogit, expected = NULL,
       pgmu=NULL, pgshape=NULL, pgextra=NULL, pgmix=NULL, wt=1,
       ndigit=10, gradtol=0.0001, steptol=gradtol,iterlim=100,...)

}
\arguments{
  \item{y}{ A matrix giving the vector of observations for each subject.}

  \item{density}{ A function giving the conditional joint or marginal
  density of the observations and calling the location, shape and extra functions.}

  \item{joint}{ Specify if the density gives the joint or common marginal
  density of the vector of observations.}

  \item{groups}{ Number of mixtures.}

  \item{gmu}{ A user-specified function of `pgmu', giving the regression 
  equation of the location parameter for each time and latent groups.
  The function should return a vector of length times*groups when
  the locations are the same across individuals or a matrix of
  n such vectors.}

  \item{gshape}{ A user-specified function of `pgshape', giving the regression
  equation of the dispersion or shape parameter for each 
  time and latent groups. The function should return a vector of 
  length times*groups or a matrix of n such vectors. 
  This function is usually parametrized such that the parameters 
  are the log of the dispersion.}

  \item{gextra}{ A user-specified function of `pgextra', giving the regression
  equation of the extra parameter for each time and latent groups.
  The function should return a vector of length times*groups or a matrix of n such vectors.}

  \item{gmixture}{ A user-specified function `pgmix', giving the regression 
  function of the mixture probabilities. The function should return 
  a vector of length groups or a matrix of n such vectors. 
  The default is the the inverse generalized logit.}

  \item{expected}{ A function returning the expected response value with
    respect to the parameters. This function will be used to compute the
    fitted values and response residuals (not deviance). By default, `gmu' will
    be taken. It is especially useful for cases where the location
    function dosn't correspond to the expected value as for censored
    normal or zero inflated poisson distributions. }

  \item{pgmu}{ Vector of initial estimates for parameters of the location function.}

  \item{pgshape}{ Vector of initial estimates for parameters of the shape function.}

  \item{pgextra}{ Vector of initial estimates for parameters of the extra function.}

  \item{pgmix}{ Vector of initial estimates for parameters of the mixture function.}

  \item{wt}{ Weight vector.}

   \item{ndigit,gradtol,steptol,iterlim,...}{ Arguments controlling `nlm'.}
}

\value{
     A list of class `moc' is returned that contains all of the
     relevant information calculated, including error codes.}

\references{
McLachlan, Geoffrey and Peel, David (2000) \emph{Finite mixture models},Wiley-Interscience, New York.

Roeder, K., Lynch, K. and Nagin, D. (1999) \emph{Modeling Uncertainty in
  Latent Class Membership: A Case Study in Criminology}, J. Amer. Statist. Assoc., \bold{94}, pp. 766--776.

Lindsay, Bruce G. and Roeder, K. (1992) \emph{Residual diagnostics for mixture models}, J. Amer. Statist. Assoc., \bold{87},
     pp. 785--794.
   }
   
\author{Bernard Boulerice <\email{Bernard.Boulerice@umontreal.ca}>}

\seealso{\code{\link{residuals.moc}},\code{\link{fitted.moc}},\code{\link{post.moc}},\code{\link{AIC.moc}},\code{\link{logLik.moc}},\code{\link{nlm}}}

\examples{

data(hyp)

# Censored Normal

cnorm<-function(x,mu,sig,min,max)
{mi<-(x==min)*1
ma<-(x==max)*1
mi*pnorm((min-mu)/sig)+ma*(1-pnorm((max-mu)/sig))+(1-mi-ma)*dnorm((x-mu)/sig)/sig}

cmean<-function(mu,sig,min,max) {
max-(max-mu)*pnorm((max-mu)/sig)+(min-mu)*pnorm((min-mu)/sig)-sig*(dnorm((max-mu)/sig)-dnorm((min-mu)/sig)) }

cnorm1<-function(x,mu,sig,...) {cnorm(x,mu,sig,0,10)}
gmu1<- function(p) {rbind(c(rep(p[1],3),rep(p[2],3),rep(p[3],3)))}
# Homogeneous variances

gshape1<- function(p) {rbind(c(rep(exp(p[1]),9)))}
cmean1<-function(p) { cmean(gmu1(p[1:3]),gshape1(p[4]),0,10) }

moc1<-
moc(hyp[,2:4],density=cnorm1,groups=3,gmu=gmu1,gshape=gshape1,expected=cmean1,
pgmu=c(2.1, 4.5, 7.8),pgshape=c(0.56),pgmix=c(0.6, -1.0), gradtol=1E-4)

\dontrun{
# Heterogeneous variances across latent groups

gshape2<-function(p) {rbind(c(rep(exp(p[1]),3),rep(exp(p[2]),3),rep(exp(p[3]),3)))}

moc2<-
moc(hyp[,2:4],density=cnorm1,groups=3,gmu=gmu1,gshape=gshape2,
pgmu=moc1$coef[1:3],pgshape=c(rep(moc1$coef[4],3)),pgmix=moc1$coef[5:6],gradtol=1E-4)
}
# Heterogeneous variances across time

gshape3<-function(p) {rbind(c(rep(c(exp(p[1]),exp(p[2]),exp(p[3])),3)))}

moc3<-
moc(hyp[,2:4],density=cnorm1,groups=3,gmu=gmu1,gshape=gshape3,
pgmu=moc1$coef[1:3],pgshape=c(rep(moc1$coef[4],3)),pgmix=moc1$coef[5:6],gradtol=1E-4)

times<-c(1.7,3,4.2)

\dontrun{
# Last class is a linear function of time

gmu1t<-function(p) {rbind(c(rep(p[1],3),rep(p[2],3),p[3]+p[4]*times))}

moc4<-
moc(hyp[,2:4],density=cnorm1,groups=3,gmu=gmu1t,gshape=gshape1,
pgmu=c(moc1$coeff[1:3],0),pgshape=c(moc1$coef[4]),pgmix=moc1$coef[5:6],gradtol=1E-4)


# Zero inflated Poisson log-linear in time for the third class

zip<- function(x,la,shape=1,extra)
{ mix<- exp(extra)/(1+exp(extra))
  mix*(x==0)+(1-mix)*dpois(x,la) }
}
gmup<-function(p) {rbind(c(rep(exp(p[1]),3),rep(exp(p[2]),3),exp(p[3]+p[4]*times)))}
\dontrun{
zipfit<-function(p) { gmup(p)/(1+exp(p[5]))}

gextrap<-function(p) {rbind(c(rep(p[1],9)))}

moc5<-
moc(hyp[,2:4],density=zip,groups=3,gmu=gmup,gextra=gextrap,expected = zipfit,
pgmu=c(0.6, 1.5, 1.9, 0),pgextra=c(-9),pgmix=c(1, -0.5), gradtol=1E-4)
}

# Standard Poisson with mixture depending on time independant covariate

dumm<-hyp[,1]-1
gmixt<-function(pm)
    {mix<-cbind(pm[1]+pm[2]*dumm,pm[3]+pm[4]*dumm)
cbind(1,exp(mix))/(1+apply(exp(mix),1,sum))}

poiss<-function(x,la,...) {dpois(x,la)}

moc6<-
moc(hyp[,2:4],density=poiss,groups=3,gmu=gmup,gmixture=gmixt,
pgmu=c(0.6, 1.5, 1.9, 0),pgmix=c(1.6,-1,0.4,-1),gradtol=1E-4)
}
\keyword{models}
