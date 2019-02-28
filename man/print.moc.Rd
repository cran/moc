\name{print.moc}
\alias{print.moc}
\alias{fitted.moc}
\alias{obsfit.moc}
\alias{coef.moc}

\title{Summary methods for fitted MOC models}

\description{
  
  \code{print.moc} prints information contained in a fitted \code{moc}
  object. The \code{attributes} \emph{parameters} of the functions
  \code{gmu}, \code{gshape}, \code{gextra} and \code{gmixture} will be
  used to label the output.

  \code{coef.moc} returns the coefficients (estimated parameters) of a
  fitted \code{moc} object.
  
  \code{fitted.moc} computes the expected values for each observation
  of a \code{moc} object using its \code{expected} function.
  
  \code{obsfit.moc} computes and prints the mean posterior
  probabilities and the posterior means of a user specified function of
  the expected and observed values, separated with respect
  to the specified variable.
}

\usage{

\method{print}{moc}(x, digits = 5, expand = TRUE, transpose = FALSE, \dots)

\method{coef}{moc}(object, split=FALSE, \dots)

\method{fitted}{moc}(object, \dots)

obsfit.moc(object, along = list(cons = rep(1, object$nsubject)),
           FUN = function(x) x)

 }

 \arguments{
   
   \item{x, object}{ Objects of class \code{moc}.}

   \item{split}{If split is TRUE, returns a list with elements
     corresponding to mu, shape, extra and mixture parameters.}
   
   \item{digits}{ Number of digits to be printed.}

   \item{expand}{ Expand density, gmu, gshape, gextra, gmixture
     function body in the print.}

   \item{transpose}{ Transpose fitted.mean and observed.mean in the print.}
   
   \item{along}{ Splitting variable. }
   
   \item{FUN}{ User defined function to apply to observed and expected values.}
   
   \item{\ldots}{ Unused.}
   
 }
 
 \details{ \code{obsfit.moc} will first compute the posterior probabilities
   for all subjects in each mixture using \code{\link{post.moc}} and
   then the weighted posterior mean probabilities 
   \deqn{\hat{\bar{\tau}}_k = \frac{\sum_i wt_i\,\hat{\tau}_{i,k}}
     {\sum_i wt_i}}{\Sum_i (wt[i] * post[i,k]) / \Sum_i wt[i]}
   The weighted posterior means of a function \eqn{g()} of the data
   (which are the empirical estimators of the conditional expectation given
   mixture group) are computed as
   \deqn{\frac{\sum_i wt_i\,\hat{\tau}_{i,k}\,g(y_i)}{\sum_i
       wt_i\,\hat{\tau}_{i,k}}}{\Sum_i (wt[i] * post[i,k] * g(y[i])) /
     \Sum_i (wt[i] * post[i,k])}
   where both sums are taken over index of valid data \eqn{y_i}{y[i]}.}

 \value{ All these methods return their results invisibly.}

 \author{Bernard Boulerice <\email{bernard.boulerice.bb@gmail.com}>}
 
 \seealso{\code{\link{moc}}, \code{\link{residuals.moc}}, \code{\link{post.moc}},
   \code{\link{plot.moc}}, \code{\link{AIC.moc}}}
 
 \keyword{methods}
 \keyword{print}
 \keyword{utilities}
