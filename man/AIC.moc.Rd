\name{AIC.moc}
\alias{AIC.moc}
\alias{logLik.moc}

\title{Information Criterion for Mixture of Curves models.}

\description{

  \code{AIC.moc} generates a table of \eqn{\log likelihood}, AIC, BIC ,
  ICL-BIC and entropy values along with the degrees of freedom of multiple \code{moc}
  objects. 

  \code{logLik} returns on object of class logLik containing the
  \eqn{\log likelihood},degrees of freedom and number of observations.

}

\usage{

       \method{AIC}{moc}(object,\dots,k=2)

       \method{logLik}{moc}(object,\dots)
}

\arguments{

  \item{object,\dots}{Objects of class \code{moc}.}

  \item{k}{ can be any real number or the string "BIC".}
}

\details{
  The computed value is \eqn{-2 \log Like + k * npar}.
  Specific treatment is carried for BIC (\eqn{k = \log(nsubject*nvariables)}), AIC
  (\eqn{k = 2}) and log likelihood (\eqn{k = 0}).
  Setting k = "BIC", will produce a table with BIC ,
  \eqn{entropy = \sum( post * \log(post) )} which is an indicator of
  mixture separation, df and \eqn{ICL-BIC = BIC + 2 * entropy} which is
  an entropy corrected BIC, see McLachlan, G. and Peel, D. (2000).
}

\value{ A data frame with the relevant information for one or more
  objects is returned .}

\note{

  Be aware that degrees of freedom (df) for mixture models are usually
  useless ( if not meaningless ) and likelihood-ratio of
  \emph{apparently} nested models often doesn't converge to a
  Chi-Square with corresponding df.
}   

  
\references{

  McLachlan, G. and Peel, D. (2000) \emph{Finite mixture
    models},Wiley-Interscience, New York.
}

\author{Bernard Boulerice <\email{Bernard.Boulerice@umontreal.ca}>}

\seealso{\code{\link{moc}}}

\keyword{methods}
\keyword{utilities}
