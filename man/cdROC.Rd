\name{cdROC}

\alias{cdROC}
\alias{cdROC.default}

\title{Cumulative/dynamic ROC curve estimate
}

\description{
This function estimates a time-dependent ROC curve following the cumulative/dynamic approach and returns a 'cdroc' object. This object can be \code{printed} or \code{plotted}. To deal with the right censored problem different statistics can be considered: those ones proposed by \emph{Martinez-Camblor et al. (2016)} based on the hazard Cox regression model (semiparametric) or the Kaplan-Meier estimator (non-parametric); and the one included in \emph{Li et al. (2016)} based on the kernel-weighted Kaplan-Meier method. See \emph{References} below.
}

\usage{
cdROC(stime, status, marker, predict.time, ...)
\S3method{cdROC}{default}(stime, status, marker, predict.time, method=c('Cox', 'KM', 'wKM'),
      kernel=c('normal', 'Epanechnikov', 'other'), h=1,
      kernel.fun = function(x,xi,h){u <- (x-xi)/h; 1/(2*h)*(abs(u) <= 1)},
      ci=FALSE, boot.n=100, conf.level=0.95, seed=2032, ...)
}

\arguments{
  \item{stime }{
vector of observed times.}
  \item{status }{
vector of status (takes the value 0 if the subject is censored and 1 otherwise).}
  \item{marker }{
vector of (bio)marker values.}
  \item{predict.time }{
considered time point (scalar).}
  \item{method }{
procedure used to estimate the probability. One of "Cox" (method based on Cox regression), "KM" (method based on Kaplan-Meier estimator) or "wKM" (method based on kernel-weighted Kaplan-Meier estimator).}
  \item{kernel }{
procedure used to calculate the kernel function. One of "normal", "Epanechnikov" or "other".  Only considered if \code{method='wKM'}.}
  \item{h }{
bandwith used to calculate the kernel function. Only considered if \code{method='wKM'}.}
  \item{kernel.fun }{
if \code{method='wKM'} and \code{kernel='other'}, function used to calculate the kernel function. It has three input parameters: \code{x}=vector, \code{xi}=value around which the kernel weight should be computed, \code{h}=bandwidth. Default: Uniform kernel.}
  \item{ci }{
if TRUE, a confidence interval for the area under the curve is computed.}
  \item{boot.n }{
number of bootstrap replicates considered to build the confidence interval. Default: 100.}
  \item{conf.level }{
the width of the confidence band as a number in (0,1). Default: 0.95, resulting in a 95\% confidence band.}
  \item{seed }{
seed considered to generate bootstrap replicates (for reproducibility).
}
  \item{\dots }{
additional arguments for \code{cdROC}. Ignored.
}
}

\details{
Assuming that larger values of the marker are associated with higher probabilities of occurrence of the event, the cumulative sensitivity and the dynamic specificity are defined by:

\eqn{Se^C(x,t) = P(}\code{marker} \eqn{> x | }\code{stime} \eqn{\le t)}   and   \eqn{Sp^D(x,t) = P(}\code{marker} \eqn{\le x |} \code{stime} \eqn{> t)}.

The resulting ROC curve is known as the cumulative/dynamic ROC curve, \eqn{R_t^{C/D}}, where \eqn{t = } \code{predict.time}.

Data censored before \eqn{t} is the major handicap with regard to the estimation of the time-dependent ROC curve. In order to estimate the probability of surviving beyond \eqn{t} for the \eqn{i}-th subject, \eqn{\hat{P}_i}, three different methods are considered:

\itemize{
\item{}{A semiparametric one, using a proportional hazard Cox regression model:}

The hazard function is estimated by \eqn{\lambda(t) = \lambda_0(t) \cdot exp(\beta \cdot X)} where \eqn{X} denotes the marker.

The probability is estimated by \eqn{\hat{P}_i = \frac{\hat{S}(t | X = x_i)}{\hat{S}(z_i | X = x_i)}} where \eqn{z_i} stands for the observed time of the \eqn{i}-th subject and \eqn{\hat{S}} is the survival function estimated from the Cox regression model.

\item{}{A non-parametric one, using the Kaplan-Meier estimator directly:}

The probability is estimated by \eqn{\hat{P}_i = \frac{\hat{S}(t)}{\hat{S}(z_i)}} where \eqn{z_i} stands for the observed time of the \eqn{i}-th subject and \eqn{\hat{S}} is the survival function estimated by the Kaplan-Meier method referred to those subjects satisfying \eqn{X \le x_i}.

\item{}{A non-parametric one, using the kernel-weighted Kaplan-Meier estimator:}

The survival function is estimated by \eqn{\hat{S}(t | X = x_i) = \prod_{s \leq t} \left[ 1- \frac{\sum_{j=1}^n K_h(x_j,x_i) I(z_j = s) status_j}{\sum_{j=1}^n K_h(x_j,x_i) I(z_j = s)} \right]} where \eqn{z_j} stands for the observed time of the \eqn{j}-th subject, \eqn{I} is the indicator function and \eqn{status_j} takes the value 0 if the \eqn{j}-th subject is censored and 1 otherwise.

Two different methods can be considered in order to define the kernel function, \eqn{K_h(x_j,x_i)}:
\itemize{
\item{}{\code{kernel='normal'}: }

\eqn{K_h(x_j,x_i) =  \frac{1}{h \sqrt{2 \pi}} exp\{ - \frac{(x_j - x_i)^2}{2 h^2} \}}

\item{}{\code{kernel='Epanechnikov'}: }

\eqn{K_h(x_j,x_i) = \frac{3}{4h} \left( 1 - \frac{x_j - x_i}{h}  \right) I(|x_j - x_i| \le h)}
}

where \eqn{h} is the bandwidth considered for kernel weights.

If the user decide to use another kernel function, \code{kernel='other'}, it should be defined by the \code{kernel.fun} input parameter, which has three parameters following this order: \code{x} is a vector, \code{xi} is the value around which the kernel weight should be computed and \code{h} is the bandwidth.

The probability is estimated by \eqn{\hat{P}_i = \frac{\hat{S}(t | X = x_i)}{\hat{S}(z_i | X = x_i)}} where \eqn{z_i} stands for the observed time of the \eqn{i}-th subject and \eqn{\hat{S}} is the survival function estimated by the kernel-weighted Kaplan-Meier method considered above.
}

}

\value{
A list of class 'cdroc' with the following content:
  \item{TP }{
vector of sensitivities (true positive rates).}
  \item{TN }{
vector of specificities (true negative rates).}
  \item{cutPoints }{
vector of thresholds considered for the (bio)marker. It coincides with the \code{marker} vector adding \eqn{min(}\code{marker}\eqn{)-1} and \eqn{max(}\code{marker}\eqn{)+1}.}
  \item{auc }{
area under the curve estimate by trapezoidal rule.}
  \item{ci }{
if TRUE, a confidence interval for the area under the curve has been computed.}
  \item{boot.n }{
number of bootstrap replicates considered to build the confidence interval. Default: 100.}
  \item{conf.level }{
the width of the confidence band as a number in (0,1). Default: 0.95, resulting in a 95\% confidence band.}
  \item{seed }{
seed considered to generate bootstrap replicates (for reproducibility).
}
  \item{meanAuc }{
bootstrap area under the curve estimate (mean along bootstrap replicates).
}
  \item{ciAuc}{
bootstrap confidence interval for the area under the curve.
}
  \item{aucs}{
vector of bootstrap area under the curve estimates.
}
  \item{stime }{
vector of observed times.}
  \item{status }{
vector of status (takes the value 0 if the subject is censored and 1 otherwise).}
  \item{marker }{
vector of (bio)marker values.}
  \item{predict.time }{
considered time point (scalar).}
  \item{method }{
procedure used in order to estimate the probability.}
  \item{kernel }{
procedure used to calculate the kernel function. Only considered if \code{method='wKM'}.}
  \item{h }{
bandwith used to calculate the kernel function. Only considered if \code{method='wKM'}.}

}

\note{
\code{survfit} and \code{Surv} functions in \code{survival} package are used in order to estimate the survival functions in both methodologies. Additionally, \code{coxph} from the same package is used to fit the Cox proportional hazard regression model in the semiparametric approach.
}

\references{
Martinez-Camblor P., F-Bayon G., Perez-Fernandez S., 2016, Cumulative/dynamic ROC curve estimation, \emph{Journal of Statistical Computation and Simulation}, \bold{86}(17), 3582-3594.

Li L., Greene T., Hu B., 2016, A simple method to estimate the time-dependent receiver operating characteristic curve and the area under the curve with right censored data, \emph{Statistical Methods in Medical Research}, DOI: 10.1177/0962280216680239.
}


\examples{
# Basic example. Data
set.seed(123)
stime <- rchisq(50,3)
status <- sample(c(rep(1,40), rep(0,10)))
marker <- max(stime) - stime + rnorm(50,0,2)

# Cumulative/dynamic ROC curve estimate at time 2.8 (Cox method is used) with 0.95 confidence
# interval for the area under the curve
cdROC(stime, status, marker, 2.8, ci=TRUE)

# Cumulative/dynamic ROC curve estimate at time 3.1 (Kaplan-Meier method is used)
cdROC(stime, status, marker, 3.1, method="KM")

# Cumulative/dynamic ROC curve estimate at time 3 (kernel-weighted Kaplan-Meier method with
# gaussian kernel and bandwidth 1 is used)
cdROC(stime, status, marker, 3, method="wKM")

# Cumulative/dynamic ROC curve estimate at time 3 (kernel-weighted Kaplan-Meier method with
# biweight kernel and bandwidth equals to 2 is used)
cdROC(stime, status, marker, 3, method="wKM", kernel="other", h=2,
      kernel.fun = function(x,xi,h){u <- (x-xi)/h; 15/(16*h)*(1-u^2)^2*(abs(u)<=1)})
}

\keyword{ time-dependent}
