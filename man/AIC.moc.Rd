\name{AIC.moc}
\alias{AIC.moc}
\alias{logLik.moc}
\alias{loglike.moc}
\alias{entropy}
\alias{entropy.default}
\alias{entropy.moc}

\title{Information criterions for MOC models.}

\description{
  \code{AIC.moc} generates a table of \eqn{\log(Likelihood)}, AIC, BIC ,
  ICL-BIC and entropy values along with the degrees of freedom of multiple \code{moc}
  objects. 

  \code{logLik} returns on object of class \code{logLik} containing the
  \eqn{\log(Likelihood)}, degrees of freedom and number of observations.
  
  \code{loglike.moc} computes the \eqn{\log(Likelihood)} of a \code{moc}
  object evaluated at the supplied parameters values, contrary to
  \code{logLik} above which uses the estimated parameter values. It gives
  the option to re-evaluate the model in which case the supplied
  parameter values are used as new starting values.
  
  \code{entropy} is a generic method to compute the entropy of
  sets of probabilities.

  The entropy of a set of \eqn{k} probabilities
  (\eqn{\pi_1,\ldots,\pi_k}{p_1,\ldots,p_k}) is computed as
  \eqn{entropy = - \sum_{i=1}^k(\pi_i \log(\pi_i) )}{entropy = - Sum_i( p_i * \log(p_i) )},
  it reaches its minimum of \eqn{0} when one of the \eqn{\pi_i=1}{p_i=1}
  (minimum uncertainty) and its maximum of \eqn{\log(k)} when all
  probabilities are equal to \eqn{\pi_i=1/k}{p_i = 1/k} (maximum uncertainty). 
  Standardized entropy is just \eqn{entropy/\log(k)} which lies in the
  interval \eqn{[0,1]}. The total and mean mixture entropy are the
  weighted sum and mean of the mixture probabilities entropy of all
  subjects. These are computed for both the prior ( without knowledge of
  the response patterns ) and the posterior mixture probabilities ( with
  knowledge of the responses).
  
  The default method \code{entropy.default} compute entropy and
  standardized entropy of a set of probabilities.
  
  \code{entropy.moc} generates a table containing weighted total and
  mean standardized entropy of prior and posterior mixture probabilities
  of \code{moc} models.  
}

\usage{

\method{AIC}{moc}(object, \dots, k = 2)

\method{logLik}{moc}(object, \dots)

loglike.moc(object, parm = object$coef, evaluate = FALSE)

\method{entropy}{moc}(object, \dots)

}

\arguments{
  
  \item{object, \ldots}{ Objects of class \code{moc}.}
  
  \item{k}{ Can be any real number or the string "BIC".}
  
  \item{parm}{ Parameters values at which the \eqn{\log(Likelihood)} is
    evaluated.}
    
  \item{evaluate}{ Boolean indicating whether re-evaluation of the model
    is desired. If \code{TRUE} \code{parm} will be used as new starting values.}
}

\details{
  The computed value in \code{AIC.moc} is \eqn{-2\cdot \log(Likelihood) +
    k\cdot npar}{-2*\log(Likelihood) + k*npar}.
  Specific treatment is carried for BIC
  (\eqn{k=\log(nsubject\cdot nvar)}{k = \log(nsubject*nvar)}),
  AIC (\eqn{k = 2}) and \eqn{\log(Likelihood)} (\eqn{k = 0}).
  Setting k = "BIC", will produce a table with BIC, mixture posterior
  \eqn{entropy =  - \sum_{i,k}( wt_i\cdot \hat{\tau}_{i,k}\,
    \log(\hat{\tau}_{i,k}) )}{entropy = - \Sum_i_k( wt[i] * post[i,k] * \log(post[i,k]) )}
  which is an indicator of mixture separation, df and
  \eqn{ICL-BIC = BIC + 2\cdot entropy}{ICL-BIC = BIC + 2 * entropy} which is
  an entropy corrected BIC, see McLachlan, G. and Peel, D. (2000) and
  Biernacki, C. et al. (2000).
}

\value{\code{AIC.moc} returns a data frame with the relevant
  information for one or more \code{moc} objects.
  
  The likelihood methods works on a single \code{moc} object:
  \code{logLik.moc} returns an object of class \code{logLik} with
  attributes \emph{df}, \emph{nobs} and \emph{moc.name} while
  \code{loglike.moc} returns a matrix containing \eqn{\log(Likelihood)}
  and corresponding estimated parameters with attributes \emph{moc.name} and
  \emph{parameters}.
  
  \code{entropy.moc} returns a \code{data.frame} with number of groups,
  total and mean standardized prior and posterior entropy of multiple
  \code{moc} objects. The percentage of reduction from prior to
  posterior entropy within a model is also supplied.} 

\note{
  Be aware that degrees of freedom (df) for mixture models are usually
  useless ( if not meaningless ) and likelihood-ratio of
  \emph{apparently} nested models often doesn't converge to a
  Chi-Square with corresponding df.
}   


\references{
  McLachlan, G. and Peel, D. (2000) \emph{Finite
    mixture models}, Wiley-Interscience, New York.

  Biernacki, C., Celeux, G., Govaert, G. (2000) \emph{Assessing a
    Mixture Model with the Integrated Completed Likelihood}, IEEE
  Transaction on Pattern Analysis and Machine Learning, \bold{22}, pp. 719--725.
}

\author{Bernard Boulerice <\email{Bernard.Boulerice@sympatico.ca}>}

\seealso{\code{\link{moc}}, \code{\link{confint.moc}}, \code{\link{profiles.postCI}},
 \code{\link{entropyplot.moc}}, \code{\link{npmle.gradient}}}

\keyword{methods}
\keyword{utilities}
