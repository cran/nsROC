\name{checkROC}

\alias{checkROC}

\title{Check data to compute an ROC curve
}

\description{
This internal function checks if the data introduced for building the curve is correct or not. It shows if there are some missing marker or response values and whether there are less or more than two levels at the response vector, D.
It also splits the data into two groups: controls and cases, depending on the corresponding value in the response vector.
}

\usage{
checkROC(X,D)
}

\arguments{
  \item{X }{
vector of (bio)marker values. It should be numeric.}
  \item{D }{
vector of response values. It should contain at least two different levels.}
}

\details{
The code will not run and an error will be showed in these cases:
\itemize{\item{}{\code{X} or \code{D} is missing,}
\item{}{\code{X} or \code{D} is NULL or full of NA's,}
\item{}{\code{X} is not a numeric vector,}
\item{}{\code{D} has less than two different values, and/or}
\item{}{\code{X} and \code{D} have different lengths.}
}
If the response vector has more than two different levels, only the two first ones are considered as controls and cases, respectively.

If the user does not agree with the codification, it can be changed modifying the order of the levels using the \code{factor} function; for instance, \code{factor(D, levels=c("1","0"))}.}

\value{
If the marker and response vectors are correct a list with the following fields is returned:
\item{levels }{levels in \code{D}. The two first ones are the labels of \code{D} considered as controls and cases, respectively.}
\item{controls }{marker values for controls.}
\item{cases }{marker values for cases.}
\item{n0 }{number of controls.}
\item{n1 }{number of cases.}
\item{X }{marker values corresponding to controls and cases (in this order).}
\item{D }{response vector consisting of \code{n0} repetitions of \code{levels[1]} and \code{n1} repetitions of \code{levels[2]} (in this order).}
}

\examples{
# Basic example with full information
set.seed(123)
X <- c(rnorm(45), rnorm(30,2,1.5))
D <- c(rep(0,45), rep(1,30))
checkROC(X,D)

# Example with some missing values and more than two levels
X <- replace(c(rnorm(25), rnorm(30,2,1.5), rnorm(20,-3,1)), seq(1,75,5), NA)
D <- replace(c(rep(0,25), rep(1,30), rep(2,20)), seq(1,75,11), NA)
checkROC(X,D)
}
