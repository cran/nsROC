\name{gROC}

\alias{gROC}
\alias{gROC.default}

\title{ROC curve estimation (generalization included)
}

\description{
This function estimates the Receiver Operating Characteristic curve and returns a 'groc' object. This object can be \code{printed} or \code{plotted}. It is the main function of the package and it includes the ROC curve generalization for non-monotone relationships (see \emph{References} below).}

\usage{
gROC(X, D, ...)
\S3method{gROC}{default}(X, D, side=c("auto","right","left","both"), Ni=1000,
     plot.roc=FALSE, plot.density=FALSE, ...)
}

\arguments{
  \item{X }{
vector of (bio)marker values.}
  \item{D }{
vector of response values.}
  \item{side }{
type of ROC curve. One of "auto", "right", "left" or "both". If the user does not specify the method "auto" is considered, hence one of "right" or "left" is assigned automatically according to the Wilcoxon statistic. If the estimate is lower than \eqn{\frac{number \; of \; controls \times number \; of \; cases}{2}}, the right-sided ROC curve is computed, otherwise the left one is considered.}
  \item{Ni }{
number of subintervals of the unit interval (FPR values) considered to calculate the curve. Default: 1000.}
  \item{plot.roc }{
if TRUE, ROC curve estimate considered is displayed.
}
  \item{plot.density }{
if TRUE, density estimates for the (bio)marker in the positive and negative subjects are displayed in the same plot.
}
  \item{\dots }{
additional arguments for \code{gROC}. Ignored.
}
}

\details{
First of all, the data introduced is checked by the \code{ROCcheck} function.

If \code{side} is not specified, one of "right" or "left" options is considered according to the comparison of the Wilcoxon test estimate and \eqn{\frac{number \; of \; controls \times number \; of \; cases}{2}}. In this case, Wilcoxon rank test is performed in order to test the alternative hypothesis \eqn{median(controls)<median(cases)} and the resulting p-value is shown but side selection is just based in the aforementioned comparison, without taking into account the p-value of the test.

The general ROC curve, \eqn{R_g(.)} is estimated using \{0, 1e-3, 2e-3, ..., 1\} as a grid for \eqn{\gamma} in the unit interval.

If both \code{plot.density} and \code{plot.roc} are \code{TRUE} they are displayed in the same window.
}

\value{
A list of class "groc" with the following content:
\item{levels }{levels in \code{D}. The two first ones are the labels of \code{D} considered controls and cases, respectively.}
\item{controls }{marker values for controls.}
\item{cases }{marker values for cases.}
\item{side }{direction of the comparison between controls and cases. One of "right", "left" or "both". If \code{side} is not specified (or it is \code{auto}) in input arguments, the output will be the one considered by the function (according to the comparison specified above).}
\item{pvalue.wilcox }{p-value of Wilcoxon test performed to compare cases and controls. Alternative hypothesis: \eqn{median(controls)<median(cases)}.}
\item{points }{vector of cut-off points of the (bio)marker considered to estimate left or right-sided ROC curves. It corresponds to the sorted marker-values, including \eqn{-\infty} and \eqn{\infty}.}
\item{pairpoints }{matrix whose rows correspond to each pair of cut-off points (xl,xu) such that xl<xu of the (bio)marker considered to estimate general ROC curves.}
\item{specificities }{vector of specificities associated to \code{points} if left or right-sided ROC curves are considered and those ones corresponding to \code{pairpoints} in case of general ROC curves. It could contain repeated values in both cases.}
\item{sensitivities }{vector of sensitivities associated to \code{points} if left or right-sided ROC curves are considered and those ones corresponding to \code{pairpoints} in case of general ROC curves. It could contain repeated values in both cases.}
\item{coordinates }{matrix whose first column is constituted by the unique values of FPR (1-\code{specificities}) and its second column is formed by the corresponding values of TPR (\code{sensitivities}). In the case of general ROC curves, several TPR values could be associated with the same FPR value.}
\item{points.coordinates }{coordinates of the points (FPR, TPR) where the ROC curve estimate has a step in case of right or left-sided curves. The first column corresponds to the values of the \code{points} vector providing these coordinates.}
\item{pairpoints.coordinates }{coordinates of the points (FPR, TPR) where the ROC curve estimate has a step in case of general curves. The first and second columns correspond to the values of the \code{pairpoints} matrix providing these coordinates.}
\item{index }{in general ROC curves, a vector displaying which row(s) of the \code{pairpoints} matrix correspond(s) to \code{pairpoints.coordinates}.}
\item{area }{area under the ROC curve estimate by trapezoidal rule.}
\item{Ni }{number of subintervals of the unit interval considered to build the curve.}
\item{ROC.t }{vector of values of \eqn{R(t)} for each \eqn{t \in} \{0, 1/\code{Ni}, 2/\code{Ni}, ..., 1\}.}
}

\references{
Martinez-Camblor P., Corral N., Rey C., Pascual J., Cernuda-Morollon E., 2014, ROC curve generalization for non-monotone relationships, \emph{Statistical Methods in Medical Research}, \bold{26}(1), 113-123.
}

\examples{
# Basic example (side="auto") -> Output side is "right"
set.seed(123)
X <- c(rnorm(45), rnorm(30,2,1.5))
D <- c(rep(0,45), rep(1,30))
gROC(X,D)

# Basic example (side="auto") -> Output side is "left"
X <- c(rnorm(45), rnorm(30,-2,1.5))
D <- c(rep(0,45), rep(1,30))
gROC(X,D)

# General ROC curve example
X <- c(rnorm(45), rnorm(30,1,4))
D <- c(rep(0,45), rep(1,30))
gROC(X, D, side="both")

# Plot density estimates and ROC curve in the same plot
X <- c(rnorm(45), rnorm(30,2,1.5))
D <- c(rep(0,45), rep(1,30))
gROC(X, D, plot.roc=TRUE, plot.density=TRUE)
}

\keyword{ generalROC}
