\name{compareROCdep}

\alias{compareROCdep}
\alias{compareROCdep.default}

\title{Comparison of k paired ROC curves}

\description{
This function compares k ROC curves from dependent data. Different statistics can be considered in order to perform the comparison: those ones included in \emph{Martinez-Camblor et al. (2013)} based on general distances between functions, the \emph{Venkatraman et al. (1996)} methodology for comparing diagnostic the accuracy of the k markers based on data from a paired design and the \emph{DeLong et al. (1988)} one based on the AUC (area under the curve) comparison. Two different methods could be considered to approximate the distribution function of the statistic: the procedure proposed by \emph{Venkatraman et al. (1996)} (based on permutated samples) or the one introduced by \emph{Martinez-Camblor et al. (2012)} (based on bootstrap samples). See \emph{References} below.
}

\usage{
compareROCdep(X, D, ...)
\S3method{compareROCdep}{default}(X, D, method=c("general.bootstrap","permutation","auc"),
              statistic=c("KS","L1","L2","CR","VK","other"),
              FUN.dist=function(g){max(abs(g))}, side=c("right","left"),
              Ni=1000, B=500, perm=500, seed=123, h.fun=function(H,x){
              H*sd(x)*length(x)^{-1/3}}, H=1, plot.roc=TRUE, type='s', lwd=3,
              lwd.curves=rep(2,ncol(X)), lty=1, lty.curves=rep(1,ncol(X)),
              col='black',col.curves=rainbow(ncol(X)), cex.lab=1.2,
              legend=c(sapply(1:ncol(X), function(i){eval(bquote(expression(
              hat(R)[.(i)](t))))}), expression(hat(R)(t))),
              legend.position='bottomright', legend.inset=0.03,
              cex.legend=1, ...)

}

\arguments{
  \item{X }{
a matrix of k columns in which each column is the vector of (bio)marker values corresponding to each sample.}
  \item{D }{
the vector of response values.}
  \item{method }{
the method used to approximate the statistic distribution. One of "general.bootstrap" (\emph{Martinez-Camblor et al. (2012)}), "permutation" (\emph{Venkatraman et al. (1996)}) or "auc" (\emph{DeLong et al. (1988)}).}
  \item{statistic }{
the statistic used to compare the curves. One of "KS" (Kolmogorov-Smirnov criteria), "L1" (\eqn{L_1}-measure), "L2" (\eqn{L_2}-measure), "CR" (Cramer-von Mises), "other" (another statistic defined by the \code{FUN.dist} input parameter), "VK" (Venkatraman) or "AUC" (area under the curve).}
  \item{FUN.dist }{
the distance considered as a function of one variable. If \code{statistic="other"} the statistic considered is \eqn{\sum_{i=1}^k} \code{FUN.dist}(\eqn{\sqrt{n_1}(\hat{R}_i(t) - \hat{R}(t)}) where \eqn{n_1} is the number of cases, \eqn{\hat{R}_i(t)} is the ROC curve estimate from the i-th sample and \eqn{\hat{R}(t) := k^{-1} \sum_{i=1}^k \hat{R}_i(t)}.}
  \item{side }{
type of ROC curve. One of "right" or "left". If \code{method="VK"} only right-sided could be considered.}
  \item{Ni }{
number of subintervals of the unit interval (FPR values) considered to calculate the curve. Default: 1000.}
  \item{B }{
number of bootstrap samples if \code{method="general.bootstrap"}. Default: 500.}
  \item{perm }{
number of permutations if \code{method="permutation"}. Default: 500.}
  \item{seed }{
seed considered to generate the permutations (for reproducibility). Default: 123.
}
  \item{h.fun }{
a function defining the bandwidth calculus used to generate the bootstrap samples if \code{method="general.bootstrap"}. It has two arguments: the first one referred to the \code{H} value and the second one, \code{x}, referred to the sample. Default: \code{function(H,x){H*sd(x)*length(x)^{-1/3}}}.}
  \item{H }{
the value used to compute \code{h.fun}, that is, the bandwidth. Default: 1.
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
the position of the legend, the inset distance from the margins as a fraction of the plot region when legend is placed and the character expansion factor relative to current \code{par("cex")}, respectively.}
  \item{...}{
another graphical parameters to be passed.
}
}

\details{
First of all, the data introduced is checked and those subjects with some missing information (marker or response value(s)) are removed. Data from a paired design should have the same length along the samples. If this is not fulfilled the code will not run and an error will be showed.

If the Venkatraman statistic is chosen in order to compare left-sided ROC curves, an error will be displayed and it will not work. The Venkatraman methodology is just implemented for right-sided ROC curves. Furthermore, for this statistics, \code{method="permutation"} is automatically assigned.

The statistic is defined by \eqn{\sum_{i=1}^k} \code{FUN.dist}(\eqn{\sqrt{n_1} \cdot (\hat{R}_i(t) - \hat{R}(t))}) where \code{FUN.dist} stands by the distance function, \eqn{n_1} is the number of cases, \eqn{\hat{R}_i(t)} is the ROC curve estimate from the i-th sample and \eqn{\hat{R}(t) := k^{-1} \sum_{i=1}^k \hat{R}_i(t)}.

The statistics implemented are defined by the following \code{FUN.dist} functions:
\itemize{
\item{}{\code{statistic="KS"}:}

\code{FUN.dist(g) = max(abs(g))}

\item{}{\code{statistic="L1"}:}

\code{FUN.dist(g) = mean(abs(g))}

\item{}{\code{statistic="L2"}:}

\code{FUN.dist(g) = mean(g^2)}

\item{}{\code{statistic="CR"}:}

\code{FUN.dist.CR(g,h) = sum(g[-length(g)]^2*(h[-1]-h[-length(h)]))}

\emph{Cramer von-Mises} statistic is defined by \eqn{\sum_{i=1}^k} \code{FUN.dist.CR}(\eqn{\sqrt{n_1} \cdot (\hat{R}_i(t) - \hat{R}(t))}, \eqn{\hat{R}(t)})
}

In case of \code{statistic="VK"} the Venkatraman methodology (see \emph{References} below) is computed to calculate the statistic. If \eqn{k>2} the statistic value is the sum of statistic values of each pair such that \eqn{i < j}.

If \code{method="general.bootstrap"} it is necessary to have a bandwidth in order to compute the bootstrap samples from the smoothed (the gaussian kernel is considered) multivariate empirical distribution functions referred to controls and cases. This bandwidth is defined by the \code{h.FUN} function whose parameters are a bandwidth constant parameter defined by the user, \code{H}, and the sample (cases or controls values of the marker) considered, \code{x}.

If \code{method="auc"}, the methodology proposed by \emph{DeLong et al.} is implemented. This option is slower because of the Mann-Whitney statistic inside requires \eqn{number~of~cases \cdot number~of~controls} comparisons. In this case, \code{statistic} returns the value of the Mann-Whitney statistic estimate and \code{test.statistic} the final test statistic estimate (formula (5) in the paper) which follows a chi-square distribution.
}

\value{
\item{n.controls }{
the number of controls.
}
\item{n.cases }{
the number of cases.
}
\item{controls.k }{
a matrix whose columns are the controls along the k samples.
}
\item{cases.k }{
a matrix whose columns are the cases along the k samples.
}
\item{statistic }{
the value of the test statistic.
}
\item{stat.boot }{
a vector of statistic values for bootstrap replicates if \code{method="general.bootstrap"}.
}
\item{stat.perm }{
a vector of statistic values for permutations if \code{method="permutation"}.
}
\item{test.statistic }{
statistic estimate given in formula (5) of \emph{DeLong et al. (1988)} (See \emph{References} below) if \code{method="auc"}.
}
\item{p.value }{
the p-value for the test.
}
}

\references{
Venkatraman E.S., Begg C.B., 1996, A distribution-free procedure for comparing receiver operating characteristic curves from a paired experiment, \emph{Biometrika}, \bold{83}(4), 835-848.

Martinez-Camblor P., Corral, N., 2012, A general bootstrap algorithm for hypothesis testing, \emph{Journal of Statistical Planning and Inference}, \bold{142}, 589-600.

Martinez-Camblor P., Carleos C., Corral N., 2013, General nonparametric ROC curve comparison, \emph{Journal of the Korean Statistical Society}, \bold{42}(1), 71-81.

DeLong E.R., DeLong D.M., Clarke-Pearson D.L., 1988, Comparing the areas under two or more correlated receiver operating characteristic curves: a nonparametric approach, \emph{Biometrics}, \bold{44}, 837-845.
}

\examples{
n0 <- 45; n1 <- 60
set.seed(123)
D <- c(rep(0,n0), rep(1,n1))

library(mvtnorm)
rho.12 <- 1/4; rho.13 <- 1/4; rho.23 <- 0.5
sd.controls <- c(1,1,1)
sd.cases <- c(1,1,1)
var.controls <- sd.controls\%*\%t(sd.controls)
var.cases <- sd.cases\%*\%t(sd.cases)
sigma.controls <- var.controls*matrix(c(1,rho.12,rho.13,rho.12,1,rho.23,rho.13,rho.23,1),3,3)
sigma.cases <- var.cases*matrix(c(1,rho.12,rho.13,rho.12,1,rho.23,rho.13,rho.23,1),3,3)
controls <- rmvnorm(n0, mean=rep(0,3), sigma=sigma.controls)
cases <- rmvnorm(n1, mean=rep(1.19,3), sigma=sigma.cases)
marker.samples <- rbind(controls,cases)

\donttest{# Default method: KS statistic proposed in Martinez-Camblor by general bootstrap
output <- compareROCdep(marker.samples, D)

# L1 statistic proposed in Martinez-Camblor by general bootstrap
output1 <- compareROCdep(marker.samples, D, statistic="L1")

# CR statistic proposed in Martinez-Camblor by permutation method
output2 <- compareROCdep(marker.samples, D, method="permutation", statistic="CR")

# Venkatraman statistic
output3 <- compareROCdep(marker.samples, D, statistic="VK")

# DeLong AUC comparison methodology
output4 <- compareROCdep(marker.samples, D, method="auc")
}
}

\keyword{ comparison}
