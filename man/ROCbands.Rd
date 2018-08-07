\name{ROCbands}

\alias{ROCbands}
\alias{ROCbands.default}

\title{Confidence bands for ROC curves}

\description{
This function computes and plots confidence bands for ROC curves (both left/right-sided and general one) using three different procedures. Particularly, one parametric approach assuming the binormal model (\emph{Demidenko}) and two non-parametric techniques (\emph{Jensen et al.} and \emph{Martinez-Camblor et al.}). See \emph{References} below.}

\usage{
ROCbands(groc, ...)
\S3method{ROCbands}{default}(groc, method = c("PSN", "JMS", "DEK"), conf.level = 0.95,
          B = 500, bootstrap.bar = TRUE, alpha1 = NULL, s = 1, a.J = NULL, b.J = NULL,
          plot.bands = FALSE, plot.var = FALSE, seed = 123, ...)
}

\arguments{
  \item{groc }{
a 'groc' object from the \code{gROC} function.}
  \item{method }{
method used to compute the confidence bands. One of "PSN" (\emph{Martinez-Camblor et al.}), "JMS" (\emph{Jensen et al.}) or "DEK" (\emph{Demidenko}).}
  \item{conf.level }{
the width of the confidence band as a number in (0,1). Default: 0.95, resulting in a 95\% confidence band.
}
  \item{B }{
number of bootstrap replicates. Default: 500 (only used in "PSN" and "JMS" methods).
}
  \item{bootstrap.bar }{
if TRUE, a bar showing bootstrap replication progress is displayed.
}
  \item{alpha1 }{
\eqn{\alpha_1} in "PSN" approach a number in (0,1) affecting the width between the lower band and the ROC curve estimate. Default: NULL, the one which minimizes the theoretical area between lower and upper bands is considered.
}
  \item{s }{
scale parameter used to compute the smoothed kernel distribution functions in "PSN" method. The bandwidth \eqn{h = s \cdot min(m,n)^{1/5} \cdot \hat{\sigma}} where \eqn{m} and \eqn{n} stand by the number of controls and cases, respectively, is considered. Default: 1.}
  \item{a.J, b.J }{
extremes of interval in (0,1) in which compute the regional confidence bands by "JMS" methodology. Default: (1/Ni, 1 - 1/Ni.).
}
  \item{plot.bands }{
if TRUE, confidence bands at level \code{conf.level} are displayed.
}
  \item{plot.var }{
if TRUE, a plot of \eqn{\sigma_n^{*,1}(t)} with \eqn{t} in [0,1] (if "PSN" method is selected) or \eqn{Var(\Psi(p))} with \eqn{p} in (\code{a.J}, \code{b.J}) (if "JMS" method is selected) is displayed.
}
  \item{seed }{
seed used to compute the bootstrap controls and cases samples in "PSN" method or Brownian Bridges in "JMS" method.
}
  \item{\dots }{
additional arguments for \code{ROCbands}. Ignored.
}
}

\details{
\itemize{
\item{}{\emph{Martinez-Camblor et al. methodology} - "PSN" \code{method}}

The \code{theoretical.area} is computed as \eqn{(c_{\alpha_1} - c_{\alpha_2}) n^{-1/2} \int \sigma_n^*(t) dt} where \eqn{\sigma_n^*(t)} is the standard deviation estimate of \eqn{\sqrt{n} [\hat{R}(\omega, .) - R(.)]} and \eqn{n} is the cases sample size.

Due to computation can take some time depending on the number of bootstrap replicates considered, a progress bar is shown.

Confidence bands are truncated in the following way: on one hand, if the lower band is lower than 0 or higher than 0.95 it is forced to be 0 or 0.95, respectively; on the other hand, if the upper band is higher than 1 or lower than 0.05 it is forced to be 1 or 0.05, respectively.

\item{}{\emph{Jensen et al. methodology} - "JMS" \code{method}}

\eqn{K^\alpha_{a,b}} denote the upper \eqn{\alpha/2}-quantile of the distribution of \eqn{\sup_{a \le p \le b} \frac{|\Psi(p)|}{\sqrt{Var \Psi(p)}}} where \eqn{(a,b)} is the interval in which the regional confidence bands are calculated and \eqn{\Psi(.)} is the limiting process of the stochastic process \eqn{\Delta_N = \sqrt{N} [\hat{R}(\omega, .) - R(.)]} with \eqn{N} being the total sample size.

Extremes of the interval (\code{a.J}, \code{b.J}) used in order to display the regional confidence bands must be divisors of \code{Ni} in the interval \eqn{[0,1]}.

Confidence bands are truncated in a similar way as in "PSN" method in order not to have bands lower than 0 or higher than 1.

\item{}{\emph{Demidenko methodology} - "DEK" \code{method}}

Demidenko ROC curve estimate does not correspond to the empirical one due to the fact that the (bio)marker values in controls and cases are supposed to come from a normal distribution is exploited.
}
}

\value{
A list of class 'rocbands' with the following content:
\item{method }{method used to compute the confidence bands. One of "PSN" (Martinez-Camblor et al.), "JMS" (Jensen et al.) or "DEK" (Demidenko).}
\item{conf.level }{the width of the confidence band as a number in (0,1).}
\item{B }{number of bootstrap replicates used in "PSN" and "JMS" methods.}
\item{L, U }{vectors containing the values of lower and upper bands, respectively, for each \eqn{t \in} \{0, 1/Ni, 2/Ni, ..., 1\}. In case of "JMS" method \code{p} is considered as \eqn{t}.}
\item{practical.area }{area between lower and upper bands (\code{L} and \code{U}) computed by trapezoidal rule.}
\item{Ni }{number of subintervals of the unit interval considered to build the curve.}
\item{ROC.t }{vector of values of \eqn{R(t)} for each \eqn{t \in} \{0, 1/Ni, 2/Ni, ..., 1\}.}

If the \code{method} is "PSN":
\item{s }{scale parameter used to compute the smoothed kernel distribution functions.}
\item{alpha1, alpha2 }{if the \code{alpha1} input argument is not specified, \eqn{\alpha_1} and \eqn{\alpha_2} values which minimize area between bands are automatically computed. If \code{alpha1} is chosen by the user, \code{alpha2} is computed by \code{alpha1} = (1 - \code{conf.level}) - \code{alpha1}.}
\item{fixed.alpha1 }{if TRUE, \code{alpha1} has been fixed by the user.}
\item{c1, c2 }{\eqn{c_{\alpha_1}} and \eqn{c_{\alpha_2}} resulting from the algorithm to compute confidence bands.}
\item{ROC.B }{matrix of size \code{Ni}+1, \code{B} whose columns contain the ROC curve estimate for each bootstrap sample.}
\item{sd.PSN }{vector \eqn{\sigma_n^*(t)} which is the estimate of the standard deviation of the empirical process considered.}
\item{theoretical.area }{theoretical area between confidence bands by trapezoidal rule.}

If the \code{method} is "JMS":
\item{a.J, b.J }{extremes of the interval in which the regional confidence bands have been computed.}
\item{p }{vector of FPR points considered in the interval (\code{a.J}, \code{b.J}).}
\item{smoothROC.p }{smooth ROC curve estimate for each value of \code{p}.}
\item{K.alpha }{value of \eqn{K_{a,b}^\alpha} computed to calculate confidence bands (see Details above).}
\item{var.JMS }{value of \eqn{Var(\Psi(p))} estimated from the formula given by Hsieh and Turnbull (see \emph{Jensen et al.} in \emph{References}).}

If the \code{method} is "DEK":
\item{DEK.fpr, DEK.tpr }{values of FPR and TPR computed to calculate the Demidenko confidence bands taking into account that it is a binormal technique.}
}

\note{
Brownian bridges needed to estimate \eqn{\Psi(.)} in "JMS" method are computed using the \code{BBridge} function in the \code{sde} package.

It should be noted that both the "PSN" and "JMS" methods are non-parametric, while the "DEK" approach is designed assuming the binormal model, so it is not convenient to use this method when distribution assumptions are not fulfilled. Furthermore, both the "JMS" and "DEK" methodologies are implemented just for the right-sided ROC curve. If \code{side} is \code{left} or \code{both} only the "PSN" method provides confidence bands.
}

\references{
Martinez-Camblor P., Perez-Fernandez S., Corral N., 2016, Efficient nonparametric confidence bands for receiver operating-characteristic curves, \emph{Statistical Methods in Medical Research}, DOI: 10.1177/0962280216672490.

Jensen K., Muller H-H., Schafer H., 2000, Regional confidence bands for ROC curves, \emph{Statistical in Medicine}, \bold{19}, 493-509.

Demidenko E., 2012, Confidence intervals and bands for the binormal ROC curve, \emph{Journal of Applied Statistics}, \bold{39}(1), 67-79.
}

\examples{
# Basic example
set.seed(123)
X <- c(rnorm(45), rnorm(30,2,1.5))
D <- c(rep(0,45), rep(1,30))
groc.obj <- gROC(X,D)

# PSN confidence bands with conf.level=0.95
ROCbands(groc.obj)
\donttest{# Plot standard deviation estimate of the curve and confidence bands in the same window
ROCbands(groc.obj, plot.bands=TRUE, plot.var=TRUE)
# PSN confidence bands with alpha1 fixed (alpha1=0.025)
ROCbands(groc.obj, alpha1=0.025)

# JMS confidence bands in (0.2,0.7) interval
ROCbands(groc.obj, method="JMS", a.J=0.2, b.J=0.7)
# Plot variance estimate of the curve and confidence bands in the same window
ROCbands(groc.obj, method="JMS", a.J=0.2, b.J=0.7, plot.bands=TRUE, plot.var=TRUE)

# DEK confidence bands with conf.level=0.99
ROCbands(groc.obj, method="DEK", conf.level=0.99)
}
}

\keyword{ confidence}
