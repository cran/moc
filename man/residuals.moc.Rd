\name{residuals.moc}
\alias{residuals.moc}
\alias{post.moc}
\alias{post}

\title{Posterior and residuals methods for MOC models.}

\description{
   
   \code{post} is a generic method for computing posterior probabilities of a fitted model.
  
   \code{post.moc} computes the posterior mixture probabilities of each subject.

   \code{residuals.moc} computes response or deviance residuals.
   The residuals are optionally weighted by the posterior mixture probabilities,
   globally ( with \emph{post} ) or within each group ( in that case
   \emph{post} is divided by its mean for each group ).

}

\usage{

   post(object,\dots)
   
   \method{post}{moc}(object,\dots)

   \method{residuals}{moc}(object,\dots,type="deviance",post.weight=TRUE,within=FALSE)

}

\arguments{

  \item{object}{Object of class \code{moc}.}
 
  \item{type}{Type of residuals: either deviance ( the default ) or
  response.}
  
  \item{post.weight}{ Specify if the residuals must be weighted by the
  posterior mixture probabilities. Weighting is preferable, it is the
  default.}

  \item{within}{ Specify if the posterior weights are rescaled within each
    mixture group.}
  
  \item{\dots}{Unused.}
}

\details{
  Response residuals are simply the difference between the
  observed and expected values, \deqn{response = y - expected}

  Deviance residuals are defined as properly scaled
  difference in the log likelihood at the observed and fitted value.
  \deqn{deviance = \sqrt{2 * wt * \log{\frac{\rm{density}(y,y,shape,extra)}{\rm{density}(y,mu,shape,extra)}}} * sign(response)}{deviance = \sqrt(2 * wt * (\log(density(y,y,shape,extra)/density(y,mu,shape,extra)))) * sign(response)}

  Globally weighted residuals are preferable to detect influential
  data, wrong number of groups and differences between groups.
  Rescaled weight residuals are more useful when plotted against some
  variables or time to detect misspecified regression function or profiles.
}  

\value{ \code{residuals.moc} returns an array of class
  \code{residuals.moc} and \code{residuals} with attributes \emph{type},
  \emph{post.weight} and \emph{within}. All these methods return their
  values invisibly.}
  
\references{

  McLachlan, G. and Peel, D. (2000) \emph{Finite mixture
  models},Wiley-Interscience, New York.
  
  
  Lindsay, B. G. and Roeder, K. (1992) \emph{Residual diagnostics for
  mixture models}, J. Amer. Statist. Assoc., \bold{87}, pp. 785--794.
}
   
\author{Bernard Boulerice <\email{Bernard.Boulerice@umontreal.ca}>}

\seealso{\code{\link{moc}},\code{\link{plot.moc}},\code{\link{print.moc}},
  \code{\link{AIC.moc}}}


\keyword{methods}
\keyword{utilities}
