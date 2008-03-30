\name{confint.moc}
\alias{confint.moc}
\alias{density.moc}
\alias{profiles.postCI}

\title{Parameter and profiles confidence intervals and likelihood
  profiling of MOC models.}

\description{
  \code{confint.moc} computes confidence intervals (CI) of specified
  function of the parameters based on crude Wald asymptotics. More
  precise CI for the original parameters are obtained through profiling
  of the likelihood function (that is evaluation of the likelihood over
  a wide range of values in the parameter space). When profiling is
  requested the deviance over different parameters' values is also
  returned.
  
  \code{profiles.postCI} computes data values for which the empirical
  probability of observing such subject values, given mixture group,
  lies between the confidence bounds (see details).

  \code{density.moc} computes the estimated mixture density at data
  points along some factor and optionally plot it.
}

\usage{
\method{confint}{moc}(object, parm = list(), level = 0.95,
        profiling = c("none","simple","complete"), \dots)

profiles.postCI(object, data = NULL, level = 0.95,
                interpolate = TRUE)

\method{density}{moc}(x, var = NULL, along = NULL,
            plot = c("none","pp-plot","density","pq-plot"),
            type = "l", \dots)
}

\arguments{
  
  \item{object, x}{ A fitted \code{moc} object.}
  
  \item{parm}{ A list of formulas beginning with \code{~} or expressions of the
    parameters denoted \eqn{p1,p2,\ldots} for which confidence intervals
    are requested.

    For example, \code{parm = list(~p1,~exp(p2)/(1+exp(p1+p3)))}.}
  
  \item{level}{ Alpha level in the interval \eqn{(0,1)} for the confidence
    bounds \eqn{[(1-level)/2,(1+level)/2]}. In \code{profiles.postCI}, you
    can also directly specify the bounds like \code{level = c(0.02,0.98)}.} 
  
  \item{profiling}{ A string that specifies the desired type of likelihood profiling.
    This can be one of
    \describe{
      \item{none:}{(default) no profiling.}
      \item{simple:}{evaluate the likelihood on a grid of values, one
	parameters at a time holding the other parameters fixed.}
      \item{complete:}{fix one parameter at a time over a grid of values
	and re-estimate the other parameters.}}}
      
  \item{data}{ An optional \code{data.frame}, \code{matrix} or \code{vector}
    of length \code{nsubject} containing the values for which the
    confidence limits are requested. The default is to take the original
    response profile.}
  
  \item{interpolate}{ A logical value indicating whether interpolation of
    data values must be performed to achieve the probability limits. When
    "FALSE", the data points with the probabilities nearest to the given
    bounds are taken (thus using the corresponding step function).}
  
  \item{var}{ A vector of integer values specifying the response
    variables for density evaluation.}
  
  \item{along}{ A factor used to split the density estimator.}
  
  \item{plot}{ A string that specifies the kind of desired plot. Allowed
    values are
    \describe{
      \item{none:}{the default.}
      \item{pp-plot:}{plot the estimated mixture cumulative
	probability (CDF) against its empirical counterpart.}
      \item{pq-plot:}{plot the estimated CDF against the quintiles.}
      \item{density:}{plot the estimated mixture distribution at
	observed data points.}}}
  
  \item{type}{ The type of lines in the plot, see \code{plot} for
    details.}
  
  \item{\dots}{ Used in \code{density.moc} to pass arguments
    directly to the plotting function. In \code{confint.moc} iterlim
    will be passed to \code{update.moc} and offscal will change the 
    profiling parameters search range.}
}
\value{
  \code{confint.moc} returns a list containing a list of arrays
  with likelihood deviance for each parameters configuration of the
  requested profiling, a function \code{ellip} corresponding to the
  asymptotic elliptic distance
  \deqn{ellip(p)=(p-\hat{p})\hat{\Sigma}^{-1}(p-\hat{p})}{ellip(p) =
    (p - p_max) S^-1 (p - p_max)}
  where \eqn{\hat{p}}{p_max} is the maximum likelihood estimator of the
  parameters and \eqn{\hat{\Sigma}}{S} its asymptotic covariance
  matrix. It also returns univariate, joint conditional and likelihood
  rejection confidence intervals for the parameters (when profiling
  has been requested).
  
  \code{profiles.postCI} returns a list of array with upper and lower
  bounds on \emph{data} profiles for each mixture group.
  
  \code{density.moc} returns nothing when a plot is requested, otherwise
  an array with mixture density estimate and data points is returned.
}

\details{
  The methods included here primarily exploit the empirical estimators
  of the conditional expectation given mixture group for some
  appropriately chosen function of the data \eqn{g()}, that is
  \deqn{\hat{g}_k=\frac{\sum_i wt_i\,\hat{\tau}_{i,k}\,g(y_i)}{\sum_i
      wt_i\,\hat{\tau}_{i,k}}.}{g_k = \Sum_i (wt[i] * post[i,k] * g(y[i])) /
    \Sum_i (wt[i] * post[i,k]).}
  Profiles confidence intervals and density estimates are defined
  by choosing \eqn{g()} as the indicator function over appropriate sets.
  See \code{\link{print.moc}} and \code{\link{residuals.moc}}.
}

\note{Mixture models are powerful tools to capture and describe the
  variability present in the data. The methods \code{profiles.postCI}
  and \code{density.moc} are especially intended to this purpose.
  The method \code{confint.moc} is quite different in essence since it
  is used to describe parameters' uncertainty that depends on sampling
  scheme and size, estimation method and goodness-of-fit of the model.
  A model with small parameters' uncertainty can poorly describe data
  variability while a model with large parameters' uncertainty can be
  very good at describing data variability.
}

\seealso{\code{\link{moc}}, \code{\link{print.moc}},
  \code{\link{residuals.moc}}, \code{\link{post.moc}},
  \code{\link{loglike.moc}}, \code{\link{profilesplot}}
}

\author{Bernard Boulerice <\email{Bernard.Boulerice@sympatico.ca}>}
\keyword{models}
\keyword{methods}
\keyword{multivariate}
\keyword{htest}
