\name{compareROCindep}

\alias{compareROCindep}
\alias{compareROCindep.default}

\title{Comparison of k independent ROC curves
}
\description{
This function compares k ROC curves from independent data. Different statistics can be considered in order to perform the comparison: those ones included in \emph{Martinez-Camblor et al. (2011)} based on distances, the \emph{Venkatraman (2000)} methodology for comparing curves for continuous unpaired data and one based in AUC (area under the curve) comparison. See \emph{References} below.}

\usage{
compareROCindep(X, G, D, ...)
\S3method{compareROCindep}{default}(X, G, D, statistic=c("L1","L2","CR","other","VK","AUC"),
                FUN.stat.int=function(roc.i, roc){mean(abs(roc.i - roc))},
                FUN.stat.cons=function(n.cases, n.controls){sqrt(n.cases)},
                side=c("right","left"), Ni=1000, raw=FALSE, perm=500,
                seed=123, plot.roc=TRUE, type='s', lwd=3,
                lwd.curves=rep(2,length(table(G))), lty=1,
                lty.curves=rep(1,length(table(G))), col='black',
                col.curves=rainbow(length(table(G))), cex.lab=1.2,
                legend=c(sapply(1:length(table(G)),function(i){
                eval(bquote(expression(hat(R)[.(i)](t))))}),
                expression(hat(R)(t))), legend.position='bottomright',
                legend.inset=0.03, cex.legend=1, ...)
}

\arguments{
  \item{X }{
vector of (bio)marker values.}
  \item{G }{
vector of group identifier values (it should have as levels as independent samples to compare).}
  \item{D }{
the vector of response values.}
  \item{statistic }{
the statistic used in order to compare the curves. One of "L1" (\eqn{L_1}-measure), "L2" (\eqn{L_2}-measure), "CR" (Cramer-von Mises), "other" (another statistic defined by \eqn{\sum_{i=1}^k} \code{FUN.stat.cons} \eqn{\cdot} \code{FUN.stat.int}), "VK" (Venkatraman) or "AUC" (area under the curve).}
  \item{FUN.stat.int }{
a function of two variables, \code{roc.i} and \code{roc} standing for ROC curve estimate for the \eqn{i-th} sample and mean ROC curve estimate along the k samples, respectively. This function represents the integral to consider in case of \code{statistic="other"}.
}
  \item{FUN.stat.cons }{
a function of two variables, \code{n.cases} and \code{n.controls} standing for the cases and controls sample size, respectively. This function represents the constant to multiply \code{FUN.stat.int} above in case of \code{statistic="other"}.}
  \item{side }{
type of ROC curve. One of "right" or "left". If \code{method="VK"} only right-sided could be considered.}
  \item{Ni }{
number of subintervals of the unit interval (FPR values) considered to calculate the curve. Default: 1000.}
  \item{raw }{
if TRUE, raw data is considered; if FALSE, data is ranked and a method to break ties in the permutations is considered (see \emph{Venkatraman (2000)} in \emph{References}). Default: FALSE.}
  \item{perm }{
number of permutations. Default: 500.}
  \item{seed }{
seed considered to generate the permutations (for reproducibility). Default: 123.
}
  \item{plot.roc }{
if TRUE, a plot including ROC curve estimates for the k samples and the mean of all of them is displayed.}
  \item{type }{
what type of plot should be drawn.}
  \item{lwd }{
the line width to be used for mean ROC curve estimate.}
  \item{lwd.curves }{
a vector with the line widths to be used for ROC curve estimates of each sample.}
  \item{lty }{
the line type to be used for mean ROC curve estimate.}
  \item{lty.curves }{
a vector with the line types to be used for ROC curve estimates of each sample.}
  \item{col }{
the color to be used for mean ROC curve estimate.}
  \item{col.curves }{
a vector with the colors to be used for ROC curve estimates of each sample.}
  \item{cex.lab }{
the magnification to be used for x and y labels relative to the current setting of \code{cex}.}
  \item{legend }{
a character or expression vector to appear in the legend.}
  \item{legend.position, legend.inset, cex.legend }{
the position of the legend, the inset distance from the margins as a fraction of the plot region when legend is placed, and the character expansion factor relative to current \code{par("cex")}, respectively.}
  \item{...}{
another graphical parameters to be passed.
}
}

\details{
If the Venkatraman statistic is chosen in order to compare left-sided ROC curves, an error will be displayed and it will not work. The Venkatraman methodology is just implemented for right-sided ROC curves.

If \code{raw=FALSE} the data will be ranked in each sample using the \code{rank} function with \code{ties.method='first'} option. Furthermore, the permutation samples possible ties will be broken using \code{ties.method='random'} option.

The statistic is defined by \eqn{\sum_{i=1}^k} \code{statistic.cons} \eqn{\cdot} \code{statistic.int} where \code{statistic.cons} = \code{FUN.stat.cons}('number of cases in the i-th sample', 'number of controls in the i-th sample') and \code{statistic.int} = \code{FUN.stat.int}('ROC curve estimate from the i-th sample', 'mean ROC curve estimate along the k samples'). It is usual to consider the function \code{FUN.stat.int} as an integral of a distance between \eqn{\hat{R}_i(t)} and \eqn{\hat{R}(t)} where \eqn{\hat{R}(t) := k^{-1} \sum_{i=1}^k \hat{R}_i(t)}.

The statistics implemented are defined by the following \code{FUN.stat.cons} and \code{FUN.stat.int} functions:
\itemize{
\item{}{\code{statistic="L1"}:}

\code{FUN.stat.int(roc.i, roc) = mean(abs(roc.i - roc))}

\code{FUN.stat.cons(n.cases, n.controls) = sqrt(n.cases)}

\item{}{\code{statistic="L2"}:}

\code{FUN.stat.int(roc.i, roc) = mean((roc.i - roc)^2)}

\code{FUN.stat.cons(n.cases, n.controls) = n.cases}

\item{}{\code{statistic="CR"}:}

\code{FUN.stat.int(roc.i, roc) = mean((roc.i[seq(2,2*Ni+1,2)] -}

\code{      roc[seq(2,2*Ni+1,2)])^2 * (roc[seq(3,2*Ni+1,2)] - roc[seq(1,2*Ni-1,2)]))}.

\code{FUN.stat.cons(n.cases, n.controls) = n.cases}

In order to use this statistic, the ROC curves have been estimated in a grid with \code{2*Ni} subintervals of the unit interval.
}

The permutation method proposed in \emph{Venkatraman (2000)} is used in order to generate the \code{perm} samples in all methodologies (i.e., any \code{statistic}).

In case of \code{statistic="VK"} the Venkatraman methodology (see \emph{References} below) is computed to calculate the statistic. If \eqn{k>2} the statistic value is the sum of the statistic values of each pair such that \eqn{i < j}.

In case of \code{statistic="AUC"}, the statistic considered is \eqn{k^{-1} \sum_{i=1}^k | \widehat{AUC}_i - \widehat{AUC} |} where \eqn{\hat{AUC}} is the mean of \eqn{\hat{AUC}_i} along the k samples.
}

\value{
\item{n.controls }{
vector of number of controls in each sample.
}
\item{n.cases }{
vector of number of cases in each sample.
}
\item{controls.k }{
a vector of all controls along the k samples, ordered by sample.
}
\item{cases.k }{
a vector of all cases along the k samples, ordered by sample.
}
\item{statistic }{
the value of the test statistic.
}
\item{stat.perm }{
a vector of statistic values for permutations.
}
\item{p.value }{
the p-value for the test.
}
}

\references{
Venkatraman E.S., 2000, A permutation test to compare receiver operating characteristic curves, \emph{Biometrics}, \bold{56}, 1134-1138.

Martinez-Camblor P., Carleos C., Corral N., 2011, Powerful nonparametric statistics to compare k independent ROC curves, \emph{Journal of Applied Statistics}, \bold{38}(7), 1317-1332.
}

\examples{
set.seed(123)
X1 <- c(rnorm(45), rnorm(30,2,1.5))
D1 <- c(rep(0,45), rep(1,30))
X2 <- c(rnorm(45), rnorm(38,3,1.5))
D2 <- c(rep(0,45), rep(1,38))
X3 <- c(rnorm(30), rnorm(42,3,1))
D3 <- c(rep(0,30), rep(1,42))
X <- c(X1, X2, X3)
D <- c(D1, D2, D3)
G <- c(rep(1,75), rep(2,83), rep(3,72))

\donttest{# Default method: L1 statistic proposed in Martinez-Camblor
output <- compareROCindep(X, G, D)

# Venkatraman statistic
output1 <- compareROCindep(X, G, D, statistic="VK")

# DeLong AUC comparison methodology
output2 <- compareROCindep(X, G, D, statistic="AUC")
}
}

\keyword{ comparison}
