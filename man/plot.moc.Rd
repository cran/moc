\name{plot.moc}
\alias{plot.moc}
\alias{profilesplot}
\alias{profilesplot.moc}
\alias{plot.residuals.moc}
\alias{mix.colors.moc}

\title{Plotting methods for MOC models.}


\description{

  \code{plot.moc} plots the fitted and observed values of a fitted
  \emph{moc} model.

  \code{profilesplot} is a generic method to plot subject profiles of fitted model.
   
  \code{profilesplot.moc} plots the response profiles of each subject
  of fitted \code{moc} object with colors that are a posterior
  probability mix of group colors.

  \code{plot.residuals.moc} nicely plots and object of type
  \code{residuals.moc}.

  \code{mix.colors.moc} computes subject colors by mixing group base colors
  accordind to the subject posterior probabilities.
}

\usage{

  \method{plot}{moc}(x,against=1:x$ntimes,main="",xlab="",ylab="",
                      prob.legend=TRUE,scale=FALSE,group.colors=rainbow(x$groups),\dots)

  \method{plot}{residuals.moc}(x,against="Index",groups=1:dim(x)[3],sunflower=FALSE,\dots)

  profilesplot(x,\dots)
  
  \method{profilesplot}{moc}(x,against=1:x$ntimes,main="",xlab="",ylab="",
                     col.legend=TRUE,scale=FALSE,
                     group.colors=rainbow(x$groups),\dots)

  mix.colors.moc(object,group.colors=rainbow(object$groups))
  
 }

\arguments{

  \item{x,object}{ Objects of class \code{moc} or \code{residuals.moc}.}
 
  \item{against}{ x axis for plotting the profiles.
      A variable against which to plot the residuals or the strings
    \describe{
  \item{Index:}{ The default, use the index of the residuals array.}
  \item{Observation:}{Use the column (variable) index of the response matrix.}
  \item{Subject:}{Use the row (subject) index of the response matrix.}
  }}

  \item{main, xlab, ylab, \dots}{ Arguments to be passed to \code{plot},
    \code{matplot}.}

  \item{prob.legend, col.legend}{ Add a mixture probabilities and
      color legend to the plot.}

  \item{scale}{ Specify if each variable should be scaled.(see \code{scale})}

  \item{group.colors}{ The groups colors to be used in the plot.}

  \item{groups}{ Specify for which groups \code{residuals.moc} plot is
    requested.}
   
  \item{sunflower}{Specify if a sunflower or standard plot is requested.}
}

\value{ \code{plot.moc} invisibly returns a \code{list} containing the
  plotted values and scaling information. \code{plot.residuals.moc}
  invisibly returns the plotted residual values. \code{mix.colors.moc} 
  invisibly returns subject mixed colors in hexadecimal RGB values.}
  
\author{Bernard Boulerice <\email{Bernard.Boulerice@umontreal.ca}>}

\seealso{\code{\link{moc}},\code{\link{residuals.moc}},\code{\link{print.moc}},
  \code{\link{AIC.moc}}}

\keyword{methods}
\keyword{hplot}

  

  
  
