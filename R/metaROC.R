metaROC <- function(data, ...) {
  UseMethod("metaROC")
}

metaROC.default <- function(data, Ni=1000, model=c("fixed-effects", "random-effects"), plot.Author=FALSE, plot.bands=TRUE, plot.inter.var=FALSE, cex.Author=0.7, lwd.Author=12, col.curve='blue', col.bands='light blue', alpha.trans=0.5, col.border='blue', ...){

  # Check if input arguments are correct

  model <- match.arg(model)
  if(!is.numeric(Ni) || length(Ni)!=1 || Ni%%1!=0 || Ni <= 0){
    stop("Ni should be a positive integer.")
  }
  if(!is.numeric(alpha.trans) || length(alpha.trans)!=1 || alpha.trans < 0 || alpha.trans > 1){
    stop("alpha.trans should be a number in the unit interval.")
  }
  if(!is.data.frame(data) || sum(names(data)=="Author")!=1 || sum(names(data)=="TP")!=1 || sum(names(data)=="TN")!=1 || sum(names(data)=="FP")!=1 || sum(names(data)=="FN")!=1){
    stop("data should be a data frame with at least 5 variables called 'Author', 'TP', 'TN', 'FP' and 'FN'.")
  }

  n <- data$n <- data$TP+data$FN
  m <- data$m <- data$FP+data$TN
  data$tpr <- data$TP/n
  data$fpr <- data$FP/m

  S <- unique(data$Author)

  cat("Number of papers included in meta-analysis:", length(S), '\n', sep=' ')

  sapply(S, function(i){
    data.S <- data[data$Author==i,]
    n.S <- dim(data.S)[1]
    if(!identical(data.S$n, rep(data.S$n[1], n.S))){
      stop("The number of positives does not coincide for some Author. Check it.")
    }else if(!identical(data.S$m, rep(data.S$m[1], n.S))){
      stop("The number of negatives does not coincide for some Author. Check it.")
    }
  })

  # Plot pairs (False-Positive Rate, True-Positive Rate) for each Author

  if(plot.Author){
    plot(data$fpr, data$tpr, xlim=c(0,1), ylim=c(0,1), lwd=lwd.Author, pch=1, col='gray', xlab="False-Positive Rate", ylab="True-Positive Rate", cex.lab=1.5, main=paste("ROC curve (",model," model)", sep=''), ...)
    axis(1, at=seq(0,1,0.01), labels=F, tck=-0.01)
    axis(1, at=seq(0,1,0.1), labels=F, tck=-0.02)
    axis(2, at=seq(0,1,0.01), labels=F, tck=-0.01)
    axis(2, at=seq(0,1,0.1), labels=F, tck=-0.02)
  }

  index.ord <- order(data$Author, data$fpr, data$tpr)
  data <- data[index.ord,]
  n.points <- dim(data)[1]

  t <- seq(0,1,length=Ni)
  roc.j <- sapply(S, function(j){
    f <- approxfun(c(0,data$fpr[data$Author==j],1), c(0,data$tpr[data$Author==j],1), ties="ordered") # If there exist a point (1,Se) with Se<1, it will not be considered
    if(plot.Author){lines(c(0,t,1), c(0,f(t),1), col='gray')}
    f(t)})
  if(plot.Author){
    for (i in 1:n.points) {
      text(data$fpr[i],data$tpr[i],data$Author[i],cex=cex.Author)
    }
  }

  # Variance of fixed-effects model ROC estimate (for each Author)

  der.roc.j <- matrix(1, Ni, length(S))

  size.j <- sapply(S, function(j){
    M <- m[min(which(data$Author==j))]
    N <- n[min(which(data$Author==j))]
    c(M,N)})

  M <- matrix(rep(size.j[1,],Ni),Ni,length(S),byrow=TRUE)
  N <- matrix(rep(size.j[2,],Ni),Ni,length(S),byrow=TRUE)
  T <- matrix(rep(t,length(S)),Ni,length(S))

  var.j <- (der.roc.j)^2*T*(1-T)/M + (roc.j)*(1-roc.j)/N

  var.j[1,] <- var.j[2,]  # To avoid infinite weights (var.j=0 iff t=0 or t=1)
  var.j[Ni,] <- var.j[(Ni-1),]

  w.j <- 1/var.j
  W <- apply(w.j, 1, sum)

  cat("Model considered:", model, '\n', sep=' ')

  # ROC curve estimate (fixed-effects model)

  RA.fem <- apply(w.j*roc.j, 1, sum)/W
  sRA.fem <- RA.fem
  for(i in 1:(Ni-1)){sRA.fem[i+1] <- max(sRA.fem[i],sRA.fem[i+1])}

  # Plot ROC curve estimate

  if(model=="fixed-effects"){

    se.RA.fem <- W^(-1/2)

    if(!plot.Author){
      plot(t, sRA.fem, 'l', xlim=c(0,1), ylim=c(0,1), lwd=1, col=col.curve, xlab="False-Positive Rate", ylab="True-Positive Rate", cex.lab=1.5, main=paste("ROC curve (",model," model)", sep=''), ...)
      abline(0,1,col='gray', lty=2)
      axis(1, at=seq(0,1,0.01), labels=F, tck=-0.01)
      axis(1, at=seq(0,1,0.1), labels=F, tck=-0.02)
      axis(2, at=seq(0,1,0.01), labels=F, tck=-0.01)
      axis(2, at=seq(0,1,0.1), labels=F, tck=-0.02)
      abline(0, 1, col='gray')
    }

    RA <- RA.fem; sRA <- sRA.fem; se.RA <- se.RA.fem

  }else if(model=="random-effects"){

    inter.var <- apply(w.j*(roc.j - RA.fem)^2, 1, sum)/W  # RA.fem is considered to compute tau^2
    inter.var[1] <- inter.var[2]
    inter.var[Ni] <- inter.var[Ni-1]

    w.j.rem <- 1/(var.j + inter.var)
    W.rem <- apply(w.j.rem, 1, sum)

    # ROC curve estimate (random-effects model)

    RA.rem <- apply(w.j.rem*roc.j, 1, sum)/W.rem
    sRA.rem <- RA.rem
    for(i in 1:(Ni-1)){sRA.rem[i+1] <- max(sRA.rem[i],sRA.rem[i+1])}

    se.RA.rem <- (apply((w.j.rem)^2*var.j, 1, sum)/(W.rem^2))^(1/2)

    if(!plot.Author){
      plot(t, sRA.rem, 'l', xlim=c(0,1), ylim=c(0,1), lwd=1, col=col.curve, xlab="False-Positive Rate", ylab="True-Positive Rate", cex.lab=1.5, main=paste("ROC curve (",model," model)", sep=''), ...)
      abline(0, 1, col='gray', lty=2)
      axis(1, at=seq(0,1,0.01), labels=F, tck=-0.01)
      axis(1, at=seq(0,1,0.1), labels=F, tck=-0.02)
      axis(2, at=seq(0,1,0.01), labels=F, tck=-0.01)
      axis(2, at=seq(0,1,0.1), labels=F, tck=-0.02)
      abline(0, 1, col='gray')
    }

    RA <- RA.rem; sRA <- sRA.rem; se.RA <- se.RA.rem

  }


  makeTransparent <- function(Color, alpha=255){
    newColor <- col2rgb(Color)
    apply(newColor, 2, function(col){
      rgb(red=col[1], green=col[2], blue=col[3], alpha=alpha, maxColorValue=255)
      })
  }

  # Plot ROC interval confidence intervals

  if(plot.bands){
    polygon(c(t,rev(t)),c(sRA-1.96*se.RA,rev(sRA+1.96*se.RA)),col=makeTransparent(col.bands,alpha.trans*255), border=makeTransparent(col.border,alpha.trans*255))
  }

  lines(c(0,t,1), c(0,sRA,1), lwd=2, col=col.curve)
  abline(0, 1, col='gray', lty=2)

  # Area under the curve estimate (AUC)

  auc <- mean(sRA[-1]+sRA[-Ni])/2
  text(0.8, 0.1, paste("AUC =", round(auc,3)))

  cat("The area under the summary ROC curve (AUC) is ", round(auc,3),".\n", sep="")

  # Youden index

  ind.youden <- which.max(1-t+sRA)
  youden.index <- c(1-t[ind.youden], sRA[ind.youden])

  cat("The optimal specificity and sensitivity (in the Youden index sense) for summary ROC curve are ", round(youden.index[1],3)," and ", round(youden.index[2],3),", respectively.\n", sep="")


  # Plot inter-study variability estimate

  if(model=="random-effects" && plot.inter.var){
    dev.new()
    par(mar=c(6,5,4,2))
    plot(t, inter.var, 'l', lwd=2, xlab="t", ylab=expression(tau[M]^2~(t)), cex.lab=1.2, main="Inter-study variability")
  }

  results <- list(data=data, t=t, sRA=sRA, RA=RA, se.RA=se.RA, youden.index=youden.index, auc=auc, roc.j=roc.j, w.j=w.j, model=model)

  if(model=="random-effects"){
    results <- list(data=data, t=t, sRA=sRA, RA=RA, se.RA=se.RA, youden.index=youden.index, auc=auc, roc.j=roc.j, w.j=w.j, w.j.rem=w.j.rem, inter.var=inter.var, model=model)
  }

  results

}
