\name{utils.moc}
\alias{inv.glogit}
\alias{glogit}
\title{MOC utility functions}
\description{ Functions to compute generalized logit and inverse logit
  with respect to a reference group.
}
\usage{
   inv.glogit(gmix,ref=1)
   glogit(p,ref=1)
}
\arguments{
  \item{p}{vector of probabilities to transform in generalized logit
  with ref as a reference group.}

  \item{gmix}{vector of generalized logit with respect to ref to
  transform in probabilities.}

  \item{ref}{reference group.}
  }
  
\value{A vector of probabilities or generalized logit.}

\details{These functions do not perform any check on their arguments.}

\author{Bernard Boulerice <\email{Bernard.Boulerice@umontreal.ca}>}
\keyword{utilities}
