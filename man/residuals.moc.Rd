\name{residuals.moc}
\alias{residuals.moc}
\alias{plot.residuals.moc}
\alias{post.moc}

\title{Residuals methods for Mixture of Curves models.}

\description{
   
  \code{residuals.moc} computes response or deviance residuals.
  The residuals are optionally weighted by the posterior mixture probabilities,
  globally ( with \emph{post} ) or within each group ( in that case
  \emph{post} is divided by its mean for each group ).

  \code{plot.residuals.moc} nicely plots and object of type
  \code{residuals.moc}.
  
  \code{post.moc} computes the posterior mixture probabilities of each subject.

}

\usage{

       \method{residuals}{moc}(object,\dots,type="deviance",post.weight=TRUE,within=FALSE)

       \method{plot}{residuals.moc}(x,against="Index",sunflower=FALSE,\dots)

       \method{post}{moc}(object,\dots)

}

\arguments{

  \item{object}{Object of class \code{moc}.}
 
  \item{type}{Type of residuals: either deviance ( the default ) or
  response.}
  
  \item{post.weight}{ Specify if the residuals must be weighted by the
  posterior mixture probabilities. Weighthing is preferable, it is the
  default.}

  \item{within}{ Specify if the posterior weights are rescaled within each
  mixture group. }

  \item{x}{An object of class \code{residuals.moc}.}

  \item{against}{A variable against which to plot the residuals or
    the strings
    \describe{
    \item{Index:}{ The default, use the index of the residuals array.}
    \item{Observation:}{Use the column (variable) index of the response matrix.}
    \item{Subject:}{Use the row (subject) index of the response matrix.}
  }}

  \item{sunflower,\dots}{Specify if a sunflower or standard plot is requested.
  Other arguments to be passed to \code{plot}.}
  
}

\details{
  Response residuals are simply the difference between the
  observed and expected values, \deqn{response = y - expected}

  Deviance residuals are defined as properly scaled
  difference in the log likelihood at the observed and fitted value.
  \deqn{deviance = \sqrt(2 * wt * (\log(density(y,y,shape,extra)/density(y,mu,shape,extra)))) * sign(response)}

  Globally weighted residuals are preferable to detect influencial
  data, wrong number of groups and differences between groups.
  Rescaled weigth residuals are more useful when plotted against some
  variables or time to detect misspecified regression function or curve.
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

\seealso{\code{\link{moc}}}

\keyword{methods}
\keyword{utilities}
