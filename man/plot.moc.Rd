\name{plot.moc}
\alias{plot.moc}
\alias{profilesplot}
\alias{profilesplot.moc}
\alias{plot.residuals.moc}
\alias{entropyplot}
\alias{entropyplot.moc}

\title{Plotting methods for MOC models.}

\description{
  
  \code{plot.moc} plots the fitted and observed values of a fitted
  \code{moc} model.
  
  \code{profilesplot} is a generic method to plot subject profiles of
  fitted model.
  
  \code{profilesplot.moc} plots the subject response profiles (variables
  or posterior) of each subject of a fitted \code{moc} object with
  colors that are a mix of group colors (see \code{mix.colors.moc}). 
  
  \code{plot.residuals.moc} nicely plots and object of type
  \code{residuals.moc}.
  
  \code{entropyplot} is a generic method to plot subject entropy.
  
  \code{entropyplot.moc} plots the prior and posterior subject entropy
  (standardized or not) corresponding to a \code{moc} model. The subject
  entropy are first sorted by prior then by posterior within prior.
}

\usage{

\method{plot}{moc}(x, against = 1:x$nvar,main = paste(substitute(x)),
     xlab = "", ylab = "", prob.legend = TRUE, scale = FALSE,
     group.colors = rainbow(x$groups), \dots)

\method{plot}{residuals.moc}(x, against = "Index", groups = 1:dim(x)[3],
     sunflower = FALSE, group.colors = NULL, \dots)

profilesplot(x, \dots)
  
\method{profilesplot}{moc}(x, against = 1:x$nvar, main = NULL,
             xlab = "", ylab = "", col.legend = TRUE,
             scale = FALSE, group.colors = rainbow(x$groups),
             type = "subject", \dots)

entropyplot(x, \dots)

\method{entropyplot}{moc}(x, main = NULL, std = TRUE, lwd = 1.5,
            col = c("red3","green3","gray95"),
            shade.gr.col = gray(1-0.5*(0:(x$groups-1))/(x$groups-1)),
            legend = TRUE, \dots)
}

\arguments{
  \item{x}{ Objects of class \code{moc} or \code{residuals.moc}.}
 
  \item{against}{ x axis for plotting the profiles.
    A variable against which to plot the residuals or the strings
    \describe{
      \item{Index:}{ The default, use the index of the residuals array.}
      \item{Observation:}{Use the column (variable) index of the response matrix.}
      \item{Subject:}{Use the row (subject) index of the response matrix.}
  }}

  \item{main, xlab, ylab, \ldots}{ Arguments to be passed to \code{plot},
    \code{matplot}.}
  
  \item{prob.legend, col.legend, legend}{ Add mixture probabilities,
    color legend and lines to the plot.}
  
  \item{sunflower}{ Specify if a sunflower or standard plot is requested.}
  
  \item{scale}{ Specify if each variable should be scaled.(see \code{\link{scale}})}
  
  \item{groups}{ Specify for which groups \code{residuals.moc} plot is
    requested.}

  \item{type}{ A string giving the type of profiles to plot
    \describe{
      \item{subject:}{(default) for subject profiles.}
      \item{variable:}{for pairs plot of response variables.}
      \item{posterior:}{for pairs of posterior probabilities. }
  }}
  
  \item{std}{ If TRUE standardized entropy are plotted.}
  
  \item{lwd, col}{ \code{entropyplot.moc} will plot the prior and posterior
    entropy with lines of width \eqn{lwd} and colors \code{col[1:2]}, the
    third color is used for the area between the two curves.}

  \item{shade.gr.col}{ \code{entropyplot.moc} will mix these colors with
    \code{mix.colors.moc} and paint the area between the two curves
    accordingly. }
  
  \item{group.colors}{ The groups base colors.}
}

\value{ \code{plot.moc} invisibly returns a \code{list} containing the
  plotted values and scaling information.

  \code{plot.residuals.moc} invisibly returns the plotted residual values.}

\author{Bernard Boulerice <\email{Bernard.Boulerice@sympatico.ca}>}

\seealso{\code{\link{moc}}, \code{\link{residuals.moc}}, \code{\link{print.moc}},
  \code{\link{AIC.moc}}}

\keyword{methods}
\keyword{hplot}

  

  
  
