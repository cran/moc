\name{print.moc}
\alias{print.moc}
\alias{fitted.moc}
\alias{obsfit.moc}
\title{Fitted MOC model}

\description{

   \code{print.moc} prints information contained in a fitted \code{moc}
   object. The \code{attributes} \emph{parameters} of the functions
   \code{gmu},\code{gshape},\code{gextra} and \code{gmixture} will be
   used to label the output.

   \code{fitted.moc} computes the expected values for each observation
   of a \code{moc} object.

   \code{obsfit.moc} computes and prints the mean posterior
   probabilities and the posterior means of a user specified function of
   the expected and observed values, splitted with respect
   to the specified variable.
 }

 \usage{

   \method{print}{moc}(x,digits=5,\dots)

   \method{fitted}{moc}(object,\dots)

   obsfit.moc(object,along=NULL,FUN=function(x) x)

 }

 \arguments{

   \item{x,object}{ Objects of class \code{moc}.}

   \item{digits}{ Number of digits to be printed.}
   
  \item{along}{ Splitting variable. }

   \item{FUN}{ User defined function to apply to observed and expected values.}

   \item{\ldots}{Unused.}
 
 }

 \details{ \code{Obsfit.moc} will first compute the posterior probabilities
   \deqn{\hat{\tau}_{i,k} = \frac{\hat{P}( G_i = k | z_i)
   \hat{h}( y_i | G_i = k, x_i)}{\sum_k \hat{P}( G_i = k | z_i)
   \hat{h}( y_i | G_i = k, x_i)}}{tau[i,k] = P( G[i] = k | z[i])
   h( y[i] | G[i] = k, x[i]) /  \Sum_k (P( G[i] = k | z[i])
   h( y[i] | G[i] = k, x[i]) )}
 using \code{post.moc} and then the weighted posterior mean probabilities
 \deqn{\hat{\bar{\tau}}_k = \frac{\sum_i wt_i\,\hat{\tau}_{i,k}}
   {\sum_i wt_i}}{\Sum_i (wt[i] * tau[i,k]) / \Sum_i wt[i]}
 The weighted posterior means of a function \eqn{f()} of the data are computed as
 \deqn{\frac{\sum_i wt_i\,\hat{\tau}_{i,k}\,f(y_i)}{\sum_i
     wt_i\,\hat{\tau}_{i,k}}}{\Sum_i (wt[i] * tau[i,k] * f(y[i])) /
   \Sum_i (wt[i] * tau[i,k])}
where both sums are taken over index of valid data \eqn{y_i}{y[i]}.}


 \value{ All these methods return their results invisibly.}

 \author{Bernard Boulerice <\email{Bernard.Boulerice@umontreal.ca}>}

 \seealso{\code{\link{moc}},\code{\link{residuals.moc}},\code{\link{plot.moc}},
 \code{\link{AIC.moc}}}

\keyword{methods}
\keyword{print}
\keyword{utilities}
