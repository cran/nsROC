\name{nsROC-package}

\alias{nsROC-package}
\alias{nsROC}

\docType{package}

\title{nsROC}

\description{
    Tools for estimating Receiver Operating Characteristic (ROC) curves,
    building confidence bands, comparing several curves both for dependent and
    independent data, estimating the cumulative-dynamic ROC curve in presence
    of censored data, and performing meta-analysis studies, among others.
}

\details{
  The basic function of the nsROC package is the \code{\link{gROC}} function. It
  will estimate an ROC curve under one of these considerations: larger values of
  the marker are associated with a higher probability of being positive (right-sided),
  the opposite (left-sided) or when both smaller and larger values of the marker
  are associated with having more probability of being positive (both).

  Confidence bands for an ROC curve estimate resulting of the previous function can
  be computed and displayed by the \code{\link{ROCbands}} function. Three different
  methods are provided to compute them.

  Several paired or unpaired ROC curves can be compared with the \code{\link{compareROCdep}}
  or \code{\link{compareROCindep}} function, respectively. In order to compare ROC curves
  different statistics can be used, and to approximate the distribution of the statistic in
  the paired case both permutation and bootstrap procedures are computed.

  Time-dependent ROC curves can be estimated by the cumulative/dynamic approach using
  the \code{\link{cdROC}} function. In order to deal with the right censored problem three
  different statistics can be considered.

  Meta-analysis of ROC curves following a non-parametric approach can be performed with
  the \code{\link{metaROC}} function. Both the fixed-effects and random-effects model can
  be considered.
}

\section{Abbreviations}{
  The following abbreviations are frequently used in this package:
  \itemize{
    \item ROC: Receiver Operating Characteristic
    \item AUC: Area Under the (ROC) Curve
    \item Sp: Specificity
    \item Se: Sensitivity
    \item TPR: True-Positive Rate
    \item FPR: False-Positive Rate
  }
}

\section{Functions}{
  \tabular{ll}{
    \code{\link{gROC}} \tab ROC curve estimate (generalization included) \cr
    \code{\link{ROCbands}} \tab Confidence bands for ROC curves \cr
    \code{\link{compareROCdep}} \tab Comparison of k paired ROC curves \cr
    \code{\link{compareROCindep}} \tab Comparison of k independent ROC curves \cr
    \code{\link{cdROC}} \tab Cumulative/dynamic ROC curve estimate \cr
    \code{\link{metaROC}} \tab Non-parametric ROC curve estimate for meta-analysis \cr
    \code{\link[=plot.groc]{plot}} \tab Plot an ROC curve \cr
    \code{\link[=plot.rocbands]{plot}} \tab Plot confidence bands for an ROC curve \cr
    \code{\link[=plot.cdroc]{plot}} \tab Plot a time-dependent ROC curve \cr
    \code{\link[=print.groc]{print}} \tab Print a \code{groc} object \cr
    \code{\link[=print.rocbands]{print}} \tab Print a \code{rocbands} object \cr
    \code{\link[=print.cdroc]{print}} \tab Print a \code{cdroc} object \cr
    \code{\link{checkROC}} \tab Check the data to compute an ROC curve (internal function)
  }
}

\section{Dataset}{
  This package comes with a dataset of 9 papers (meta-analysis) with the number of TP
  (true positive), FP (false positive), TN (true negative) and FN (false negative)
  about the use of the Interleukin6 (IL6) as a marker for the early detection of
  neonatal sepsis: \code{\link{interleukin6}}.
}

\section{Installing and using}{
  To install this package:
  \preformatted{
    install.packages("nsROC")
  }

  To load the package:
  \preformatted{
    library(nsROC)
  }
}


\author{
Sonia Perez-Fernandez

Maintainer: Sonia Perez Fernandez <uo217889@uniovi.es>
}

\seealso{
CRAN packages \pkg{sde} and \pkg{survival} employed in this package.
}

\keyword{package}

