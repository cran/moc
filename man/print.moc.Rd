\name{print.moc}
\alias{print.moc}
\alias{plot.moc}
\alias{fitted.moc}
\alias{obsfit.moc}
\title{Fitted Mixture of Curves model}

\description{

   \code{print.moc} prints information contained in a fitted \code{moc}
   object. The \code{attributes} \emph{parameters} of the functions
   \code{gmu},\code{gshape},\code{gextra} and \code{gmixture} will be
   used to label the output.

   \code{plot.moc} plots the fitted and observed values of a fitted
   \emph{moc} model.

   \code{fitted.moc} computes the expected values for each observation of a \code{moc} object.

   \code{obsfit.moc} computes and prints the mean posterior probabilities and the mean
   of a user specified function of the expected and observed values,
   splitted with respect to the specified variable.

 }

 \usage{

   \method{print}{moc}(x,digits=5,\dots)

   \method{plot}{moc}(x,times=1:x$ntimes,main="",xlab="",ylab="",\dots)

   \method{fitted}{moc}(object,\dots)

   obsfit.moc(object,along=NULL,FUN=function(x) x)

 }

 \arguments{

   \item{x,object}{ Objects of class \code{moc}.}

   \item{digits}{ Number of digits to be printed.}
   
   \item{times}{ x axis for plotting the curves. }

   \item{main, xlab, ylab, \dots}{ Arguments to be passed to \code{matplot}.}

   \item{along}{ Splitting variable. }

   \item{FUN}{ User defined function to apply to observed and expected values.}
  
 }

 \value{ All these methods return their results invisibly.}

 \author{Bernard Boulerice <\email{Bernard.Boulerice@umontreal.ca}>}

 \seealso{\code{\link{moc}}}

\keyword{methods}
\keyword{print}
\keyword{hplot}
\keyword{utilities}
