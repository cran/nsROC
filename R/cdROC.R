# The only required package is "survival" to use the functions "coxph" (to fit proportional hazards regression model) and "survfit" (to create survival curves)

# install.packages("survival")
# library(survival)

cdROC <- function(stime, status, marker, predict.time, ...) {
  UseMethod("cdROC")
}

cdROC.default <- function(stime, status, marker, predict.time, method=c('Cox', 'KM', 'wKM'),
                          kernel=c('normal', 'Epanechnikov', 'other'), h=1,
                          kernel.fun = function(x,xi,h){u <- (x-xi)/h; 1/(2*h)*(abs(u) <= 1)},
                          ci=FALSE, boot.n=100, conf.level=0.95, seed=2032, ...){

  # Check if input arguments are correct

  method <- match.arg(method)
  kernel <- match.arg(kernel)

  if(!is.numeric(stime) || length(stime)<=1){
    stop("stime should be a vector of length higher than 1.")
  }else if(length(status)<=1){
    stop("status should be a vector of length higher than 1.")
  }else if(!is.numeric(marker) || length(marker)<=1){
    stop("marker should be a vector of length higher than 1.")
  }

  if(is.numeric(status)){
    who2 <- !is.na(status)
    temp <- (status == 0 | status == 1)
    status <- ifelse(temp, status, NA)
    if(!all(temp[who2], na.rm = TRUE)){
      warning("Invalid status value, converted to NA.")
    }
  }else{
    stop("Invalid status vector, must be numeric (0 = censored, 1 = not censored).")
  }

  if(length(stime)!=length(marker) | length(marker)!=length(status)){
    stop("stime, status and marker should be vectors with same length.")
  }

  if(!is.numeric(predict.time) || length(predict.time)!=1 || predict.time<0){
    stop("predict.time should be a non-negative number.")
  }else if(max(stime,na.rm=TRUE)<predict.time){
    stop("predict.time should be lower inside interval of times considered.")
  }

  if(method=='wKM' && (is.null(h) ||!is.numeric(h) || length(h)!=1 || h < 0)){
    stop("h should be a non-negative number.")
  }

  if(!is.numeric(boot.n) || length(boot.n)!=1 || boot.n%%1!=0 || boot.n <= 0){
    stop("boot.n should be a positive integer.")
  }

  if(!is.numeric(conf.level) || length(conf.level)!=1){
    stop("conf.level should be a number.")
  }else if(conf.level > 1 && conf.level <= 100){
    conf.level <- conf.level/100
  }
  if(conf.level < 0 || conf.level > 1){
    stop("conf.level should be in the unit interval.")
  }

  # Function which assigns the probability of belonging to the negative group (controls) to those subjects censored before predict.time

  assignProbability <- function(stime, status, marker, predict.time, method=c('Cox', 'KM', 'wKM'),
                                kernel=c('normal', 'Epanechnikov', 'other'), h=1,
                                kernel.fun = function(x,xi,h){u <- (x-xi)/h; 1/(2*h)*(abs(u) <= 1)},
                                undefinedIndices){

    method <- match.arg(method)
    kernel <- match.arg(kernel)
    results <- 1:length(undefinedIndices)

    if(method == 'Cox'){

      fit <- survival::coxph(survival::Surv(stime, status)~marker)
      md <- survival::survfit(fit, newdata=data.frame(cbind(stime,status,marker)))
      for(j in 1:length(undefinedIndices)){
        f <- approxfun(c(min(md$time)-1, md$time, max(md$time)+1), c(1, md$surv[, undefinedIndices[j]], 0))
        results[j] <- f(predict.time) / f(stime[undefinedIndices[j]])
        if(f(stime[undefinedIndices[j]]) == 0){
          results[j] <- 0
        }
      }

    }else if(method=='KM'){

      for(j in 1:length(undefinedIndices)){
        idx <- which(marker <= marker[undefinedIndices[j]])
        fit <- survival::survfit(survival::Surv(stime[idx], status[idx])~1)
        f <- stepfun(fit$time, c(1, fit$surv))
        results[j] <- f(predict.time) / f(stime[undefinedIndices[j]])
        if(f(stime[undefinedIndices[j]]) == 0){
          results[j] <- 0
        }
      }

    }else if(method=='wKM'){

      if(kernel!='other'){
        kernel.option <- function(x,xi,h,kernel){
          u <- (x-xi)/h
          switch(kernel,
                 normal = 1/(h*sqrt(2*pi))*exp(-u^2/2),
                 Epanechnikov = 3/(4*h)*(1-u^2)*(abs(u) <= 1)
          )
        }
        kernel.function <- function(x,xi,h){kernel.option(x,xi,h,kernel=kernel)}
      }else{
        kernel.function <- kernel.fun
      }

      for(j in 1:length(undefinedIndices)){
        weights <- kernel.function(marker, marker[undefinedIndices[j]], h)
        fit <- survival::survfit(survival::Surv(stime, status)~1, weights=weights)
        f <- stepfun(fit$time, c(1, fit$surv))
        results[j] <- f(predict.time) / f(stime[undefinedIndices[j]])
        if(summary(fit, times=stime[undefinedIndices[j]])$surv == 0){
          results[j] <- 0
        }
      }

    }

    return(results)

  }


  # Function which assigns the specificity, sensitivity and AUC estimates besides other aspects of cumulative/dynamic ROC curve estimate

  singleCdroc <- function(stime, status, marker, predict.time, method=c('Cox', 'KM', 'wKM'),
                          kernel=c('normal', 'Epanechnikov', 'other'), h=1,
                          kernel.fun = function(x,xi,h){u <- (x-xi)/h; 1/(2*h)*(abs(u) <= 1)},
                          nsp.auc=1000){

    method <- match.arg(method)
    kernel <- match.arg(kernel)

    if(method!='wKM'){
      kernel <- NULL
      h <- NULL
    }

    positiveIndices <- which(stime <= predict.time & status == 1)
    negativeIndices <- which(stime > predict.time)
    undefinedIndices <- which(stime <= predict.time & status == 0)

    undefinedProb <- NULL
    if(length(undefinedIndices) > 0){
      undefinedProb <- assignProbability(stime, status, marker, predict.time, method, kernel, h, kernel.fun, undefinedIndices)
    }

    cutPoints <- c(min(marker) - 1, unique(sort(marker)), max(marker) + 1)
    nSens <- length(positiveIndices) + sum(1 - undefinedProb)
    nSpec <- length(negativeIndices) + sum(undefinedProb)
    sensitivity <- NULL
    specificity <- NULL
    for(i in 1:length(cutPoints)){
      sensitivity[i] <- (sum(marker[positiveIndices] > cutPoints[i]) + sum(1 - undefinedProb[marker[undefinedIndices] > cutPoints[i]])) / nSens
      specificity[i] <- (sum(marker[negativeIndices] <= cutPoints[i]) + sum(undefinedProb[marker[undefinedIndices] <= cutPoints[i]])) / nSpec
    }

    I <- order(1 - specificity, sensitivity)
    roc.x <- c(0,1 - specificity[I],1)
    roc.y <- c(0,sensitivity[I],1)
    nsp.auc <- length(roc.x)
    auc <- sum((roc.x[-1] - roc.x[-nsp.auc])*(roc.y[-1] + roc.y[-nsp.auc])/2)

    results <- list(TP=sensitivity, TN=specificity, undefinedProb=undefinedProb, cutPoints=cutPoints, auc=auc, predict.time=predict.time, method=method, kernel=kernel, h=h)

    attr(results, 'class') <- 'cdroc'

    return(results)

  }

  # Main function

  if(ci){
    set.seed(seed)
    sampledIndices <- sapply(1:boot.n, function(xx) sample(length(stime), length(stime), replace=TRUE))
    allResults <- lapply(1:boot.n, function(index, sampledIndices){
      currentIndices <- sampledIndices[, index]
      singleCdroc(stime[currentIndices], status[currentIndices], marker[currentIndices], predict.time, method, kernel, h, kernel.fun)
      }, sampledIndices=sampledIndices)

    result <- singleCdroc(stime, status, marker, predict.time, method, kernel, h, kernel.fun)

    allAucs <- sapply(allResults, function(xx) xx$auc)
    result$meanAuc <- mean(allAucs)
    result$ciAuc <- quantile(allAucs, c(1 - conf.level, conf.level))
    result$ci <- TRUE
    result$boot.n <- boot.n
    result$conf.level <- conf.level
    result$seed <- seed
    result$aucs <- allAucs

    return(result)

  }else{

    return(singleCdroc(stime, status, marker, predict.time, method, kernel, h, kernel.fun))

  }

}

