\name{utils.moc}
\alias{inv.glogit}
\alias{glogit}
\alias{mix.colors.moc}
\alias{mocUtils}

\title{MOC utility functions.}
\description{ Functions to compute generalized logit and inverse logit
  with respect to a reference group.
  
  \code{mix.colors.moc} computes subject colors useful for plotting by
  mixing group base colors according to the subject posterior
  probabilities.

  \code{mocUtils} lists or returns an environment (which can be
  attached) containing the utility functions from \emph{Utils}
  subdirectory of \emph{moc} package.
  
}

\usage{
inv.glogit(gmix, ref = 1)
glogit(p, ref = 1)

mix.colors.moc(object, group.colors = rainbow(object$groups))

mocUtils(filename)
}
\arguments{
  \item{object}{A fitted moc object.}
  
  \item{p}{ Vector of probabilities to transform in generalized logit
    (log-odds) with ref as a reference group.}

  \item{gmix}{ Vector of generalized logit with respect to ref to
    transform in probabilities.}
  
  \item{ref}{ Reference group.}
  
  \item{group.colors}{ The groups base colors to be mixed in proportion
    corresponding to posterior probabilities.}

  \item{filename}{ The filename containing the source of utility functions or
    empty for a list of available source codes.}
}

\value{A vector of probabilities or generalized logit.
  \code{mix.colors.moc} invisibly returns a vector of subject mixed colors in
  hexadecimal \emph{RGB} format.}

\details{These functions do not perform any check on their arguments.}

\note{More utility functions are available in the \emph{Utils} directory
  of the \code{moc} distribution.
}
\author{Bernard Boulerice <\email{bernard.boulerice.bb@gmail.com}>}
\keyword{utilities}
