\name{print}

\alias{print.groc}
\alias{print.rocbands}
\alias{print.cdroc}

\title{
Print an ROC curve object}

\description{
This function prints a 'groc', 'rocbands' or 'cdroc' object.
}

\usage{
\S3method{print}{groc}(x, ...)
\S3method{print}{rocbands}(x, ...)
\S3method{print}{cdroc}(x, ...)
}

\arguments{
  \item{x }{
a 'groc', 'rocbands' or 'cdroc' object from the \code{gROC}, \code{ROCbands} or \code{cdROC} respectively.}
  \item{\dots}{
     further arguments to be passed to other methods. Ignored.
}
}

\value{
These functions return information about the object they were passed.
}

\seealso{
\code{\link{gROC}}, \code{\link{ROCbands}}, \code{\link{cdROC}}
}

\examples{
# Data generation
set.seed(123)
X <- c(rnorm(45), rnorm(30,2,1.5))
D <- c(rep(0,45), rep(1,30))

# Print a groc object
grocobj <- gROC(X,D)
print(grocobj)

# Print an rocbands object
grocobj <- ROCbands(grocobj)
print(grocobj)

# Print a cdroc object
set.seed(123)
stime <- rchisq(50,3)
status <- sample(c(rep(1,45), rep(0,5)))
marker <- max(stime) - stime + rnorm(50,0,2)
cdrocobj <- cdROC(stime, status, marker, 3, ci=TRUE)
print(cdrocobj)
}

\keyword{ print}% use one of  RShowDoc("KEYWORDS")
