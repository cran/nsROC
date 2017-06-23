\name{metaROC}

\alias{metaROC}
\alias{metaROC.default}

\title{Non-parametric ROC curve estimate for meta-analysis}

\description{
This function performs meta-analytic studies of diagnostic tests for both the fixed and random-effects models. In particular it reports a fully non-parametric ROC curve estimate when data come from a meta-analysis study using the information of all cut-off points available in the selected original studies. The approach considered is the one proposed by \emph{Martinez-Camblor et al. (2017)} based on weighting each individual interpolated ROC curve. See \emph{References} below.}

\usage{
metaROC(data, ...)
\S3method{metaROC}{default}(data, Ni=1000, model=c("fixed-effects","random-effects"),
        plot.Author=FALSE, plot.bands=TRUE, plot.inter.var=FALSE,
        cex.Author=0.7, lwd.Author=12, col.curve='blue',
        col.bands='light blue', alpha.trans=0.5, col.border='blue', ...)
}

\arguments{
  \item{data}{
  a data frame containing at least the following variables (with these names):
  \itemize{
\item{}{\code{Author}: a vector assigning different numbers to each paper/author.}
\item{}{\code{TP}: true positives.}
\item{}{\code{FP}: false positives.}
\item{}{\code{TN}: true negatives.}
\item{}{\code{FN}: false negatives.}
}
}
  \item{Ni}{
  number of points of the unit interval (FPR values) considered to calculate the curve. Default: 1000.}
  \item{model}{
  the meta-analysis model used to estimate the ROC curve. One of "fixed-effects" (it only considers the within-study variability) or "random-effects" (it takes into account the variability between the studies).
}
  \item{plot.Author}{
  if TRUE, a plot including ROC curve estimates (by linear interpolation) for each paper under study is displayed.
}
  \item{plot.bands}{
  if TRUE, confidence interval estimate for the curve is added to the plot of the ROC curve estimate.
}
  \item{plot.inter.var}{
  if TRUE, a plot including inter-study variability estimate is displayed on an additional window.
}
  \item{cex.Author}{
  the magnification to be used to display the paper/author points labels relative to the current setting of \code{cex}.
}
  \item{lwd.Author}{
  the size to be used for the paper/author points.
}
  \item{col.curve}{
  the color to be used for the (summary) ROC curve estimate. Default: blue.
}
  \item{col.bands}{
  the color to be used for the confidence interval of ROC curve estimate. Default: light blue.
}
  \item{alpha.trans}{
  proportion of opacity to be used for the confidence interval of ROC curve estimate. A number in the unit interval where 0 means transparent. Default: 0.5.
}
  \item{col.border}{
  the color to be used for the border of confidence interval of ROC curve estimate. Default: blue.
}
  \item{...}{
another graphical parameters to be passed.
}
}

\details{
The slight modification considered to ensure the monotonicity of the summary ROC curve estimate is the following
\eqn{sRA(t) = max(sup_{z \in [0,t]} sRA(z), RA(t)).}

Some basic information about the model used and the results obtained are printed.
}

\value{
  \item{data}{
  the data-frame considered ordered by Author-FPR-TPR and including the following variables:    \itemize{
\item{}{\code{n}: positive subjects sample size.}
\item{}{\code{m}: negative subjects sample size.}
\item{}{\code{FPR}: false positive rate.}
\item{}{\code{TPR}: true positive rate.}
}
}
  \item{t}{
  values of the unit interval (FPR values) considered to calculate the curve.
}
  \item{model}{
  the meta-analysis model used to estimate the ROC curve. One of "fixed-effects" (it only considers the within-study variability) or "random-effects" (it takes into account the variability between the studies).
}
  \item{sRA}{
  non-parametric summary ROC curve estimate following the \code{model} considered with a slight modification to ensure the monotonicity. This is the one reported in graphics.
}
  \item{RA}{
  non-parametric summary ROC curve estimate following the \code{model} without modifications.
}
  \item{se.RA}{
  standard-error of summary ROC curve estimate.
}
  \item{area}{
  area under the summary ROC curve estimate by trapezoidal rule.
}
  \item{youden.index}{
  the optimal specificity and sensitivity (in the Youden index sense).
}
  \item{roc.j}{
  a matrix whose column j contains the estimated ROC curve for the j-th study in each point \code{t} considered.
}
  \item{w.j}{
  a matrix whose column j contains the weights in fixed-effects model for the j-th study in each point \code{t} considered.
}
  \item{w.j.rem}{
  a matrix whose column j contains the weights in random-effects model for the j-th study in each point \code{t} considered.
}
  \item{inter.var}{
  inter-study variability estimate in each point \code{t} considered. Only computed if \code{model} is \code{"random-effects"}.
}
}

\references{
Martinez-Camblor P., 2017, Fully non-parametric receiver operating characteristic curve estimation for random-effects meta-analysis, \emph{Statistical Methods in Medical Research}, \bold{26}(1), 5-20.}

\examples{
data(interleukin6)

# Fixed-effects meta-analysis showing linear interpolations of the papers considered in the graphic
output1 <- metaROC(interleukin6, plot.Author=TRUE)

# Random-effects meta-analysis displaying also a window with a plot of the inter-study
# variability estimate
output2 <- metaROC(interleukin6, model="random-effects", plot.Author=TRUE)
}


\keyword{ meta-analysis}
