\name{residuals.moc}
\alias{residuals.moc}
\alias{fitted.moc}
\alias{post.moc}
\alias{print.moc}
\alias{AIC.moc}
\alias{logLik.moc}
\alias{obsfit.moc}
\title{Methods for Mixture of Curves models.}

\description{

  \code{print.moc} prints major information contained in a fitted \code{moc}
  object. The \code{attributes} \emph{parameters} of the functions
  \code{gmu},\code{gshape},\code{gextra} and \code{gmixture} will be
  used to label the output whenever they are set.
   
  \code{residuals.moc} computes \eqn{response = y - expected} or\cr
  \eqn{deviance = \sqrt(2 * wt *
    (\log(density(y,y,shape,extra)/density(y,mu,shape,extra)))) * sign(response)}
  \crresiduals of an object of class \code{moc}. These residuals are then
  optionally weighted by the posterior mixture probabilities,
  globally ( with \emph{post} ) or within each group ( in that case
  \emph{post} is divided by its mean for each group ).

  \code{post.moc} computes the posterior mixture probabilities of each subject.

  \code{fitted.moc} computes the expected  values for each observation of a \code{moc} object.

  \code{obsfit.moc} computes the mean posterior probabilities and the mean
  of a user specified function of the expected and observed values,
  splitted with respect to the specified variable. 
  
  \code{AIC.moc} generates a table of \eqn{\log likelihood}, AIC, BIC ,
  ICL-BIC values along with the degrees of freedom of multiple \code{moc}
  objects. The computed values are \eqn{-2 \log Like + k * npar}.
  For BIC \eqn{k = \log(nobs)}. Setting k = "BIC", will produce a table
  with BIC , \eqn{entropy = \sum post * \log(post)} which is an
  indicator of mixtures separation, df and \eqn{ICL-BIC = BIC + 2 * entropy}
  which is an entropy corrected BIC, see McLachlan, G. and Peel, D. (2000).

  \code{logLik} returns on object of class logLik containing the
  \eqn{\log likelihood},degrees of freedom and number of observations.

}
\usage{

       print.moc(x,...)

       residuals.moc(object,...,type="deviance",post.weight=TRUE,within=FALSE)

       post.moc(object,...)

       fitted.moc(object,...)

       obsfit.moc(object,along=NULL,FUN=function(x) x)

       AIC.moc(object,...,k=2)

       logLik.moc(object,...)

}
\arguments{
  \item{x,object,...}{ Objects of class `moc'.}

  \item{type}{ Either deviance ( the default ) or response. Deviance
    residuals are defined as properly scaled difference in the log
    likelihood at the observed and fitted value. Response residuals are
    simply the difference between the observed and expected values.}

  \item{post.weight}{ Specify if the residuals must be weighted by the
    posterior mixture probabilities. Weighthing is preferable, it is the
  default.}

  \item{within}{ Specify if the weights are rescaled within each
    mixture group. Global residuals are preferable to detect influencial
  data,wrong number of groups and differences between groups. Rescaled
  residuals are more useful when plotted against some variables or time
  to detect misspecified regression function or curve.}

  \item{along}{ Splitting variable. }

  \item{FUN}{ User defined function to apply to observed and expected values.}
  
  \item{k}{ can be any real number or the string "BIC".}

}

\references{

  McLachlan, G. and Peel, D. (2000) \emph{Finite mixture
    models},Wiley-Interscience, New York.
  
  Lindsay, B. G. and Roeder, K. (1992) \emph{Residual diagnostics for mixture models}, J. Amer. Statist. Assoc., \bold{87},
     pp. 785--794.
   }
   
\author{Bernard Boulerice <\email{Bernard.Boulerice@umontreal.ca}>}

\seealso{\code{\link{moc}}}

\keyword{methods}
