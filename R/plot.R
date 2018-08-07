plot.groc <- function(x, lwd = 2, xlab = "False-Positive Rate", ylab = "True-Positive Rate", main = "ROC curve", ...){

  obj <- x
  type <- ifelse(obj$param, 'l', 's')

  if(is.null(obj$Ni)) xx <- c(0, obj$t, 1) else xx <- c(0, seq(0, 1, 1/obj$Ni), 1)
  plot(xx, c(0, obj$roc, 1), type = type, xlim = c(0,1), ylim = c(0,1), lwd = lwd, xlab = xlab, ylab = ylab, main = main, ...)
  axis(1, at=seq(0,1,0.01), labels=F, tck=-0.01)
  axis(1, at=seq(0.1,0.9,0.1), labels=F, tck=-0.02)
  axis(2, at=seq(0,1,0.01), labels=F, tck=-0.01)
  axis(2, at=seq(0.1,0.9,0.1), labels=F, tck=-0.02)
  abline(0, 1, lty=2, col="gray")

}


plot.rocbands <- function (x, type = "s", lwd = 2, xlim = c(0, 1), ylim = c(0,1),
                           xlab = "False-Positive Rate", ylab = "True-Positive Rate",
                           main = paste("ROC curve \n (", obj$method, " confidence bands)", sep = ""),
                           col = "aquamarine3", col.inside = "azure2", col.frontier = "azure3",
                           lwd.frontier = 2, ...)
{
  obj <- x
  Ni <- obj$Ni
  if(is.null(Ni)) Ni <- length(obj$controls)
  if (obj$method == "PSN") {
    plot(c(0, seq(0, 1, 1/Ni), 1), c(0, obj$ROC.t, 1),
         type = type, lwd = lwd, xlab = xlab, ylab = ylab,
         main = main, ...)
    axis(1, at = seq(0, 1, 0.01), labels = F, tck = -0.01)
    axis(1, at = seq(0.1, 0.9, 0.1), labels = F, tck = -0.02)
    axis(2, at = seq(0, 1, 0.01), labels = F, tck = -0.01)
    axis(2, at = seq(0.1, 0.9, 0.1), labels = F, tck = -0.02)
    lines(rep(seq(0, 1, 1/Ni), each = 2), c(rbind(obj$L,
                                                  obj$U)), lwd = 2, col = col.inside)
    lines(c(0, seq(0, 1, 1/Ni), 1), c(0, obj$ROC.t, 1),
          type = type, col = col, lwd = lwd)
    lines(seq(0, 1, 1/Ni), obj$L, type = type, col = col.frontier,
          lwd = lwd.frontier)
    lines(seq(0, 1, 1/Ni), obj$U, type = type, col = col.frontier,
          lwd = lwd.frontier)
    abline(0, 1, lty = 2, col = "gray")
  }
  else {
    if (obj$method == "JMS") {
      plot(0, 0, type = type, lwd = lwd, xlab = xlab, ylab = ylab,
           xlim = xlim, ylim = ylim, main = main, ...)
      axis(1, at = seq(0, 1, 0.01), labels = F, tck = -0.01)
      axis(1, at = seq(0.1, 0.9, 0.1), labels = F, tck = -0.02)
      axis(2, at = seq(0, 1, 0.01), labels = F, tck = -0.01)
      axis(2, at = seq(0.1, 0.9, 0.1), labels = F, tck = -0.02)
      abline(0, 1, lty = 1, col = "gray")
      lines(0, 0, type = type, lwd = lwd, col = "white")
      lines(rep(obj$p, each = 2), c(rbind(obj$L, obj$U)),
            lwd = 2, col = col.inside)
      lines(obj$p, obj$U, type = type, col = col.frontier,
            lwd = lwd.frontier)
      lines(obj$p, obj$L, type = type, col = col.frontier,
            lwd = lwd.frontier)
      lines(obj$p, obj$smoothROC.p, type = type, col = col,
            lwd = lwd)
      abline(v = obj$a.J, lty = 2, col = "gray")
      abline(v = obj$b.J, lty = 2, col = "gray")
      abline(0, 1, lty = 2, col = "gray")
    }
    else {
      plot(c(0, sort(obj$DEK.fpr), 1), c(0, sort(obj$DEK.tpr),
                                         1), type = type, lwd = lwd, xlab = xlab, ylab = ylab,
           main = main, ...)
      axis(1, at = seq(0, 1, 0.01), labels = F, tck = -0.01)
      axis(1, at = seq(0.1, 0.9, 0.1), labels = F, tck = -0.02)
      axis(2, at = seq(0, 1, 0.01), labels = F, tck = -0.01)
      axis(2, at = seq(0.1, 0.9, 0.1), labels = F, tck = -0.02)
      abline(0, 1, lty = 1, col = "gray")
      segments(obj$DEK.fpr, pnorm(obj$y[, 1]), obj$DEK.fpr,
               pnorm(obj$y[, 2]), col = col.inside, pch = 20,
               lwd = 2)
      segments(pnorm(obj$x[, 1]), obj$DEK.tpr, pnorm(obj$x[,
                                                           2]), obj$DEK.tpr, col = col.inside, pch = 20,
               lwd = 2)
      points(obj$DEK.fpr, obj$DEK.tpr, col = col.inside,
             pch = 20)
      lines(c(0, sort(obj$DEK.fpr), 1), c(0, sort(obj$DEK.tpr),
                                          1), type = "l", col = col, lwd = lwd)
      lines(obj$inferior.band.points, pch = 19, type = "l",
            col = col.frontier, lwd = lwd.frontier)
      lines(obj$superior.band.points, pch = 19, type = "l",
            col = col.frontier, lwd = lwd.frontier)
      abline(0, 1, lty = 2, col = "gray")
    }
  }
}


plot.cdroc <- function(x, type='s', lwd=3, xlab='1 - Specificity', ylab='Sensitivity', xaxs='i', yaxs='i', main=paste("ROC curve at time", obj$predict.time), ...){

  obj <- x

  fpsorted <- c(0,sort(1 - obj$TN),1)
  tpsorted <- c(0,sort(obj$TP),1)
  plot(fpsorted, tpsorted, type=type, lwd=lwd, xlab=xlab, ylab=ylab, xaxs=xaxs, yaxs=yaxs, main=main, ...)
  abline(0, 1, lty=2)

}
