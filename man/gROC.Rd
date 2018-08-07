\name{gROC}

\alias{gROC}
\alias{gROC.default}

\title{ROC curve estimation (generalization included)
}

\description{
This function estimates the Receiver Operating Characteristic curve and returns a 'groc' object. This object can be \code{printed} or \code{plotted}. It is the main function of the package and it includes the ROC curve generalization for non-monotone relationships (see \emph{References} below).}

\usage{
gROC(X, D, ...)
\S3method{gROC}{default}(X, D, side=c("right", "left", "both", "both2", "auto"),
          Ni = NULL,  plot.roc = FALSE, plot.density = FALSE,
          pval.auc = FALSE, B = 500, ...)
}

\arguments{
  \item{X }{
vector of (bio)marker values.}
  \item{D }{
vector of response values.}
  \item{side }{
type of ROC curve. One of "right", "left", "both", "both2" or "auto". If the user does not specify the method, "right" is considered, i.e., the right-sided ROC curve is estimated. If "auto", one of "right" or "left" is assigned automatically according to the Wilcoxon statistic. If the estimate is lower than \eqn{\frac{number \; of \; controls \times number \; of \; cases}{2}}, the right-sided ROC curve is computed, otherwise the left one is considered. If "both", the gROC curve is estimated in the usual way (controls inside the interval) and if "both2", the opposite direction is considered.}
  \item{Ni }{
number of subintervals of the unit interval (FPR values) considered to calculate the curve. Default: NULL (in this case the fast algorithm considering as many FPR values as number of controls is considered).}
  \item{plot.roc }{
if TRUE, ROC curve estimate considered is displayed.
}
  \item{plot.density }{
if TRUE, density estimates for the (bio)marker in the positive and negative subjects are displayed in the same plot.
}
  \item{pval.auc }{
if TRUE, a permutation test to test if the AUC differs from $0.5$ is performed. Available if Ni is NULL. Default: FALSE.
}
  \item{B }{
number of permutations used for testing. Default: 500.
}
\item{\dots }{
additional arguments for \code{gROC}. Ignored.
}
}

\details{
First of all, the data introduced is checked by the \code{ROCcheck} function.

If \code{side} is not specified, one of "right" or "left" options is considered according to the comparison of the Wilcoxon test estimate and \eqn{\frac{number \; of \; controls \times number \; of \; cases}{2}}. In this case, Wilcoxon rank test is performed in order to test the alternative hypothesis \eqn{median(controls)<median(cases)} and the resulting p-value is shown but side selection is just based in the aforementioned comparison, without taking into account the p-value of the test.

If \code{Ni} is \code{NULL}, the general ROC curve, \eqn{R_g(.)} is estimated considering every different pair (t, \eqn{\gamma} \eqn{\cdot} t) available on data. Otherwise it is estimated using \{0, 1e-3, 2e-3, ..., 1\} as a grid for \eqn{\gamma} in the unit interval.

If both \code{plot.density} and \code{plot.roc} are \code{TRUE} they are displayed in the same window.
}

\value{
A list of class "groc" with the following content:
\item{levels }{levels in \code{D}. The two first ones are the labels of \code{D} considered controls and cases, respectively.}
\item{controls }{marker values for controls.}
\item{cases }{marker values for cases.}
\item{side }{direction of the comparison between controls and cases. One of "right", "left", "both" or "both2". If \code{side} is \code{auto} in input arguments, the output will be the one considered by the function (according to the comparison specified above).}
\item{pvalue.wilcox }{p-value of Wilcoxon test performed to compare cases and controls. Alternative hypothesis: \eqn{median(controls)<median(cases)}.}
\item{t }{vector of values of t: \{0, 1/\code{Ni}, 2/\code{Ni}, ..., 1\}.}
\item{roc }{vector of values of \eqn{R(t)} for each \eqn{t}.}
\item{auc }{area under the ROC curve estimate by trapezoidal rule.}
\item{pval.auc }{p-value of the permutation test over the AUC.}
\item{Paucs }{different permutation AUCs displayed if the hypothesis test is performed.}
\item{points.coordinates }{coordinates of the points (FPR, TPR) where the ROC curve estimate has a step in case of right or left-sided curves. The first column corresponds to the values of the \code{points} vector providing these coordinates.}
\item{pairpoints.coordinates }{coordinates of the points (FPR, TPR) where the ROC curve estimate has a step in case of general curves. The first and second columns correspond to the values of the \code{pairpoints} matrix providing these coordinates.}
\item{param }{a logical value indicating if the estimation procedure whether parametric or not.}
\item{Ni }{number of subintervals of the unit interval considered to build the curve.}
\item{points }{if \code{Ni} != \code{NULL}, vector of cut-off points of the (bio)marker considered to estimate left or right-sided ROC curves. It corresponds to the sorted marker-values, including \eqn{-\infty} and \eqn{\infty}.}
\item{pairpoints }{if \code{Ni} != \code{NULL}, matrix whose rows correspond to each pair of cut-off points (xl,xu) such that xl<xu of the (bio)marker considered to estimate general ROC curves.}
\item{specificities }{if \code{Ni} != \code{NULL}, vector of specificities associated to \code{points} if left or right-sided ROC curves are considered and those ones corresponding to \code{pairpoints} in case of general ROC curves. It could contain repeated values in both cases.}
\item{sensitivities }{if \code{Ni} != \code{NULL}, vector of sensitivities associated to \code{points} if left or right-sided ROC curves are considered and those ones corresponding to \code{pairpoints} in case of general ROC curves. It could contain repeated values in both cases.}
\item{coordinates }{if \code{Ni} != \code{NULL}, matrix whose first column is constituted by the unique values of FPR (1-\code{specificities}) and its second column is formed by the corresponding values of TPR (\code{sensitivities}). In the case of general ROC curves, several TPR values could be associated with the same FPR value.}
\item{index }{if \code{Ni} != \code{NULL}, in general ROC curves, a vector displaying which row(s) of the \code{pairpoints} matrix correspond(s) to \code{pairpoints.coordinates}.}
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
