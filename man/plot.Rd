\name{plot}

\alias{plot.groc}
\alias{plot.rocbands}
\alias{plot.cdroc}

\title{
Plot an ROC curve object}

\description{
This function plots a 'groc', 'rocbands' or 'cdroc' object.
}

\usage{
\S3method{plot}{groc}(x, lwd = 2, xlab = "False-Positive Rate",
      ylab = "True-Positive Rate", main = "ROC curve", ...)
\S3method{plot}{rocbands}(x, type='s', lwd=2, xlim=c(0,1), ylim=c(0,1),
     xlab="False-Positive Rate", ylab="True-Positive Rate",
     main=paste("ROC curve \n (", obj$method, " confidence bands)", sep=""),
     col='aquamarine3', col.inside="azure2", col.frontier="azure3",
     lwd.frontier=2, ...)
\S3method{plot}{cdroc}(x, type='s', lwd=3, xlab='1 - Specificity',
     ylab='Sensitivity', xaxs='i', yaxs='i',
     main=paste("ROC curve at time", obj$predict.time), ...)
}

\arguments{
  \item{x }{
a 'groc', 'rocbands' or 'cdroc' object from the \code{gROC}, \code{ROCbands} or \code{cdROC} respectively.}
  \item{type }{
what type of plot should be drawn.}
  \item{lwd }{
the line width to be used for ROC curve estimate, a positive number. See \code{par}.}
  \item{col }{
the color to be used for ROC curve estimate. See \code{par}.}
  \item{lwd.frontier }{
the line width to be used for ROC curve confidence bands estimate.
}
  \item{col.inside, col.frontier }{
the color to be used for ROC curve confidence bands estimate (\code{col.frontier}) and for the area inside (\code{col.inside}).}
  \item{xlim, ylim }{
numeric vectors of length 2, giving the x and y coordinates ranges. See \code{plot.window}.}
  \item{xlab, ylab }{
a title for the x and y axis, respectively. See \code{title}.}
  \item{xaxs, yaxs }{
the style of axis interval calculation to be used for the x and y axis, respectively. See \code{par}.}
  \item{main }{
an overall title for the plot. See \code{title}.}
  \item{\dots }{
further arguments to be passed to methods, such as graphical parameters. See \code{\link[graphics]{par}}.}
}

\value{
These functions return a plot of the object they were passed.
}

\examples{
# Data generation
set.seed(123)
X <- c(rnorm(45), rnorm(30,2,1.5))
D <- c(rep(0,45), rep(1,30))

# Plot an ROC curve
grocobj <- gROC(X,D)
plot(grocobj)

# Plot ROC curve confidence bands
rocbandsobj <- ROCbands(grocobj)
plot(rocbandsobj)

# Plot cumulative/dynamic ROC curve
set.seed(123)
stime <- rchisq(50,3)
status <- sample(c(rep(1,40), rep(0,10)))
marker <- max(stime) - stime + rnorm(50,0,2)
cdrocobj <- cdROC(stime, status, marker, 2.8, ci=TRUE)
plot(cdrocobj)
}


\keyword{ plot}
