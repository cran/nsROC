\name{interleukin6}
\alias{interleukin6}
\docType{data}
\title{
Interleukin 6 (IL6) data}
\description{
This dataset includes the true-positives (TP), false-positives (FP), true-negatives (TN) and false-negatives (FN) reported by 9 different papers which study the use of the Interleukin 6 (IL6) as a marker for the early detection of neonatal sepsis.
}
\usage{interleukin6}
\format{
  A data frame with 19 observations of the following 5 variables.
  \describe{
    \item{\code{Author}}{a vector assigning different numbers to each paper}
    \item{\code{TP}}{vector of true positives}
    \item{\code{FP}}{vector of false positives}
    \item{\code{FN}}{vector of false negatives}
    \item{\code{TN}}{vector of true negatives}
  }
}
\details{
In those papers providing more than one pair of Sensitivity-Specificity all of them are collected.
}

\references{
Martinez-Camblor P., 2017, Fully non-parametric receiver operating characteristic curve estimation for random-effects meta-analysis, \emph{Statistical Methods in Medical Research}, \bold{26}(1), 5-20.
}
\examples{
# Load the dataset
data(interleukin6)

# Plot pairs (FPR, TPR) for each Author

attach(interleukin6)

TPR <- TP/(TP+FN)
FPR <- FP/(FP+TN)
plot(FPR, TPR, xlim=c(0,1), ylim=c(0,1), lwd=10, pch=1, col='gray', xlab="False-Positive Rate",
    ylab="True-Positive Rate", main=paste("ROC curve interpolation"))

S <- unique(Author)
ind <- order(Author, FPR, TPR)
ord.data <- cbind(Author[ind], FPR[ind], TPR[ind])
roc.j <- sapply(S, function(j){
  lines(c(0,ord.data[Author==j,2],1), c(0,ord.data[Author==j,3],1), col='gray')})
for(i in 1:19){text(ord.data[i,2],ord.data[i,3],ord.data[i,1],cex=0.5)}
}

\keyword{datasets}
