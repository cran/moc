\name{residuals.moc}
\alias{residuals.moc}
\alias{fitted.moc}
\alias{post.moc}
\alias{print.moc}
\alias{AIC.moc}
\alias{logLik.moc}
\title{Methods for Mixture of Curves.}

\description{

  \code{print.moc} computes the posterior mixture probabilities of each subject.

  \code{residuals.moc} computes the response or deviance residuals of an
  object of class `moc'. These can be weighted by the posterior mixture
  probabilities within each group or globally.

  \code{post.moc} computes the posterior mixture probabilities of each subject.

  \code{fitted.moc} computes the expected  values for each observation of a `moc' object.

  \code{AIC.moc} generates a table of \eqn{\log likelihood}, AIC or BIC values
  along with the degrees of freedom of multiple `moc' objects. The
  computed values are \eqn{-2 \log Like + k * npar}. For the BIC \eqn{k = \log(nobs)}.

  \code{logLik} returns on object of class logLik containing the \eqn{\log likelihood}, degrees of freedom and number of observations.

}
\usage{

       print.moc(x,...)

       residuals.moc(object,...,type="deviance",post.weight=TRUE,within=TRUE)

       post.moc(object,...)

       fitted.moc(object,...)

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
    mixture group.}

  \item{k}{ can be any real number or the string "BIC".}

}

\references{

Lindsay, Bruce G. and Roeder, K. (1992) \emph{Residual diagnostics for mixture models}, J. Amer. Statist. Assoc., \bold{87},
     pp. 785--794.
   }
   
\author{Bernard Boulerice <\email{Bernard.Boulerice@umontreal.ca}>}

\seealso{\code{\link{moc}}}

\keyword{methods}
