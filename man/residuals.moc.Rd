\name{residuals.moc}
\alias{residuals.moc}
\alias{post.moc}
\alias{post}
\alias{npmle.gradient}

\title{Posterior and residuals methods for MOC models}

\description{
  \code{post} is a generic method for computing posterior probabilities of a fitted model.
  
  \code{post.moc} computes the posterior mixture probabilities for each
  subject of a fitted \code{moc} model.
  
  \code{residuals.moc} computes response, deviance, gradient and mixture
  residuals. The residuals are optionally weighted by the posterior
  mixture probabilities, globally ( with \code{post} ) or within each
  group ( in that case \code{post} is divided by its mean for each group
  ).
  
  \code{npmle.gradient} computes the components and average of the
  mixture gradient function at some specified parameters values
  \deqn{\frac{1}{n}\sum_{i=1}^n\Bigl{(}\frac{\hat{h}_j(y_i|\,x_i)}
    {\sum_k\hat{\pi}_k(z_i,x_i)\,\hat{h}_k(y_i|\,x_i)} - 1\Bigr{)}}{
    1/n Sum_i { h_k(y[i] | x[i] ) / Sum_k p_k(z[i],x[i]) h_k( y[i] | x[i]) - 1}}
  or the empirical mixture gradient which is
  \deqn{\frac{1}{n}\sum_{i=1}^n\Bigl{(}\frac{d\hat{F}(y_i)}{\sum_k\hat{\pi}_k(z_i,x_i)\,
      \hat{h}_k(y_i|\,x_i)} - 1\Bigr{)}}{
    1/n Sum_i { dF(y[i]) / Sum_k p_k(z[i],x[i]) h_k( y[i] | x[i]) - 1}} 
  where \eqn{d\hat{F}()}{dF()} is the empirical measure.
}

\usage{

post(object, \dots)
   
\method{post}{moc}(object, \dots)

\method{residuals}{moc}(object, \dots,
          type = c("deviance","response","mixture","gradient"),
          post.weight = TRUE, within = FALSE)

npmle.gradient(object, parm = object$coef, gradient = TRUE,
               average = FALSE)
}

\arguments{
  
  \item{object}{ Object of class \code{moc}.}
  
  \item{type}{Type of residuals: either "deviance" ( the default ),
    "response", "gradient" and "mixture" (see description and details).}
  
  \item{post.weight}{ Specify if the residuals must be weighted by the
    posterior mixture probabilities. Weighting is preferable, it is the
    default.}
  
  \item{within}{ Boolean that specify if the posterior weights are
    rescaled within each mixture group.}
  
  \item{parm}{ Parameters values at which evaluation of the
    gradient takes place.}
  
  \item{gradient}{ Boolean specifying if you require the density for each
    mixture group (\emph{TRUE}) or the empirical measure (\emph{FALSE}) in
    the numerator of \code{npmle.gradient} ( see the description above ).}
  
  
  \item{average}{ Boolean that specify if \code{npmle.gradient} should return
    the individual components or the average.}
  
  \item{\ldots}{ Unused.}
}

\details{
  Posterior probabilities are the conditional probabilities of mixture
  groups given the subject response observations and are computed using the
  formula:
  \deqn{\hat{\tau}_{i,k} = \frac{\hat{\pi}_k(z_i,x_i)\,
      \hat{h}_k( y_i | x_i)}{\sum_k \hat{\pi}_k(z_i, x_i)\,
      \hat{h}_k( y_i | x_i)}}{post[i,k] = p_k( z[i], x[i])
    h_k( y[i] | x[i]) /  \Sum_k (p_k(z[i], x[i])
    h_k( y[i] | x[i]) )}
  
  Response residuals are simply the difference between the
  observed and expected values, \deqn{response = y - expected}
  
  Deviance residuals are defined as properly scaled
  difference in the log likelihood at the observed and fitted value.
  \deqn{deviance = \Bigl{\lbrace}2\cdot wt\cdot
      \log{\Bigl{(}\frac{\mathrm{density}(y,y,shape,extra)}
	{\mathrm{density}(y,mu,shape,extra)}\Bigr{)}}\Bigr{\rbrace}^{1/2} \cdot \mathrm{sign}(response)}{
    deviance = \sqrt(2*wt*(\log(density(y,y,shape,extra)/density(y,mu,shape,extra)))) * sign(response)}
  
  The \code{npmle.gradient} function is primarily intended to compute the
  components which are used to define the mixture and gradient residuals in
  \code{residuals.moc}. Mixture residuals uses the empirical gradient
  components while gradient residuals uses the mixture gradient
  components defined in the description above. The average of the first
  form above which we call the gradient can be used to check some
  minimal requirement about a a proposed solution given with parm: if it
  is a maximum likelihood (not necessarily a non-parametric maximum
  likelihood) the average should be \eqn{0}.
  See Lindsay, B.G. (1983) for details about the use of the gradient
  function for finite mixture and non-parametric mixture and Lindsay,
  B.G. and Roeder, K. (1992) for residuals diagnostics specific to mixture
  models.
 
  Globally weighted residuals are preferable to detect influential
  data, wrong number of groups and differences between groups.
  Rescaled weight residuals are more useful when plotted against some
  variables or variable index to detect misspecified regression function
  or density.
}  

\value{ \code{residuals.moc} returns an array of class
  \code{residuals.moc} and \code{residuals} with attributes \emph{type},
  \emph{post.weight} and \emph{within}. All these methods return their
  values invisibly.}

\references{
  McLachlan, G. and Peel, D. (2000) \emph{Finite mixture
    models}, Wiley-Interscience, New York.
  
  Lindsay, B. G. and Roeder, K. (1992) \emph{Residual diagnostics for
    mixture models}, Journal of the American Statistical Association,
  \bold{87}, pp. 785--794.  
}
   
\author{Bernard Boulerice <\email{Bernard.Boulerice@sympatico.ca}>}

\seealso{\code{\link{moc}}, \code{\link{plot.moc}}, \code{\link{print.moc}},
  \code{\link{AIC.moc}}, \code{\link{plot.residuals.moc}}}


\keyword{methods}
\keyword{models}
\keyword{utilities}
