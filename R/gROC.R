gROC <- function(X, D, ...) {
  UseMethod("gROC")
}

gROC.default <- function(X, D, side=c("right", "left", "both", "both2", "auto"), Ni = NULL,  plot.roc = FALSE, plot.density = FALSE, pval.auc = FALSE, B = 500, ...){

  data <- checkROC(X, D)
  levels.names <- levels(as.factor(D))
  # controls <- split(X,D)[[levels.names[1]]]; cases <- split(X,D)[[levels.names[2]]]
  # D <- ifelse(as.factor(D)==levels.names[1], 0, 1); levels <- levels(as.factor(D))
  levels <- data$levels; controls <- data$controls; cases <- data$cases
  m <- n0 <- data$n0; n1 <- data$n1; X <- data$X; D <- data$D

  if(plot.density){
    densityControls <- density(controls)
    densityCases <- density(cases)
    minDensMarker <- min(c(densityControls$x, densityCases$x))
    maxDensMarker <- max(c(densityControls$x, densityCases$x))
    minDens <- min(c(densityControls$y, densityCases$y))
    maxDens <- max(c(densityControls$y, densityCases$y))
    if(plot.roc){par(mfrow = c(1, 2))}else{par(mfrow = c(1, 1))}
    plot(densityControls, col = 3, lwd = 2, xlim = c(minDensMarker, maxDensMarker), ylim = c(minDens, maxDens), xlab = "Marker values", ylab = "", cex.lab = 1.2, main = "Density estimation")
    lines(densityCases, col = 2, lwd = 2)
    legend("topright", legend = c("Controls", "Cases"), inset = 0.03, lwd = c(2, 2), lty = c(1, 1), col = c(3, 2))
  }

  side <- match.arg(side)
  pvalue.wilcox <- NULL
  if(side == "auto"){
    wilcox <- wilcox.test(controls, cases, "less")
    if(wilcox$statistic <= n0 * n1/2){
      side <- "right"
      pvalue.wilcox <- wilcox$p.value
    }else{
      side <- "left"
      pvalue.wilcox <- wilcox$p.value
    }
  }

  Ni <- match.arg(Ni)
  if(is.null(Ni)){

    #m <- length(controls)

    t <- seq(0,1,1/m)
    N <- length(t)

    XX <- sort(c(controls,cases))
    e <- ifelse(length(unique(XX))>1, min(unique(XX)[-1]-unique(XX)[-length(unique(XX))])/2, sqrt(.Machine$double.eps))

    main <- function(side){
      if(side=='right'){
        c <- as.numeric(quantile(controls,1-t,type=3))
        roc <- 1-ecdf(cases)(c)
        results <- list(roc=roc, c=c)
      }
      if(side=='left'){
        c <- as.numeric(quantile(controls,t,type=3)) - e*as.numeric(ecdf(controls)(quantile(controls,t,type=3)) > t)
        roc <- ecdf(cases)(c)
        results <- list(roc=roc, c=c)
      }
      if(side=='both'){
        A <- sapply(1:N, function(i){
          if(i == N){
            roc <- 1; xl <- xu <- max(controls)
          }else{
            gamma <- seq(1,i,1)
            index.gamma.t <- which.max(ecdf(cases)(sort(controls)[gamma]-e) + 1 - ecdf(cases)(sort(controls)[m-i+gamma]))
            gamma.t <- gamma[index.gamma.t]
            xl <-  sort(controls)[gamma.t]; xu <- sort(controls)[m-i+gamma.t]
            roc <- ecdf(cases)(xl-e) + 1 - ecdf(cases)(xu)
          }
          c(roc, xl, xu)
        })
        results <- list(roc=A[1,], xl=A[2,], xu=A[3,])
      }
      if(side=='both2'){
        A <- sapply(1:N, function(i){
          if(i >= m){
            if(i == m){
              xl.opt <- c(min(min(controls)-e, min(cases)), min(controls))
              xu.opt <- c(max(controls), max(max(controls)+e, max(cases)))
              index.opt.t <- which.max(ecdf(cases)(xu.opt-e) - ecdf(cases)(xl.opt))
              xl <- xl.opt[index.opt.t]; xu <- xu.opt[index.opt.t]
            }
            if(i == N){
              xl <- min(min(controls)-e, min(cases)); xu <- max(max(controls)+e, max(cases))
            }
          }else{
            gamma <- seq(1,m-i,1)
            index.gamma.t <- which.max(ecdf(cases)(sort(controls)[gamma+i]-e) - ecdf(cases)(sort(controls)[gamma]))
            gamma.t <- gamma[index.gamma.t]
            xl <-  sort(controls)[gamma.t]; xu <- sort(controls)[gamma.t+i]
          }
          roc <- max(ecdf(cases)(xu-e) - ecdf(cases)(xl), 0)
          c(roc, xl, xu)
        })
        results <- list(roc=A[1,], xl=A[2,], xu=A[3,])
      }
      return(results)
    }

    mainres <- main(side)

    roc <- mainres$roc
    if(side=='right' || side=='left'){
      c <- mainres$c
      points.coordinates <- cbind(c, t, roc)
      colnames(points.coordinates) <- c("c", "FPR", "TPR")
    }
    if(side=='both' || side=='both2'){
      xl <- mainres$xl; xu <- mainres$xu
      pairpoints.coordinates <- cbind(xl, xu,  t, roc)
      colnames(pairpoints.coordinates) <- c("xl", "xu", "FPR", "TPR")
    }

    auc <- sum(roc[-N]*(t[-1] - t[-N]))

    Paucs <- NULL
    if(pval.auc){
      D <- c(rep(0, n0), rep(1, n1))
      cat("Progress bar of permutated iterations for computing pval.auc (B = ", B, ")\n", sep = "")
      bar <- txtProgressBar(min = 0, max = B, style = 3)
      for(i in 1:B){
        setTxtProgressBar(bar, i)
        pD <- sample(D, size = n0 + n1, replace = FALSE)
        Paucs[i] <- gROC(X, pD, side = side)$auc
      }
      close(bar)
      pval.auc <- mean(Paucs > auc)
    }

    if(side=='right' || side=='left'){
      results <- list(levels = levels.names, controls = controls, cases = cases, side = side, pvalue.wilcox = pvalue.wilcox, t = t, roc = roc, auc = auc, pval.auc = pval.auc, Paucs = Paucs, points.coordinates = points.coordinates, param = FALSE, Ni = Ni)
    }
    if(side=='both' || side=='both2'){
      results <- list(levels = levels.names, controls = controls, cases = cases, side = side, pvalue.wilcox = pvalue.wilcox, t = t, roc = roc, auc = auc, pval.auc = pval.auc, Paucs = Paucs, pairpoints.coordinates = pairpoints.coordinates, param = FALSE, Ni = Ni)
    }
    attr(results, "class") <- "groc"

  }else{

    points <- c(-Inf, unique(sort(X)), Inf)
    grid <- seq(0, 1, 1/Ni)
    gamma <- seq(0, 1, 0.001)
    sol <- function(side, controls, cases) {
      switch(side, right = ecdf(1 - ecdf(controls)(cases))(grid), left = ecdf(ecdf(controls)(cases))(grid), both = apply(matrix(1 - ecdf(1 - ecdf(controls)(cases))(1 - gamma %*% t(grid)) + ecdf(1 - ecdf(controls)(cases))((1 - gamma) %*% t(grid)), length(gamma), length(grid)), 2, max))
    }
    ROC.t <- sol(side = side, controls = controls, cases = cases)
    if (side == "right") {
      specificities <- sapply(points, function(c) {
        sum(controls <= c)/length(controls)
      })
      sensitivities <- sapply(points, function(c) {
        sum(cases > c)/length(cases)
      })
    }
    else if (side == "left") {
      specificities <- sapply(points, function(c) {
        sum(controls > c)/length(controls)
      })
      sensitivities <- sapply(points, function(c) {
        sum(cases <= c)/length(cases)
      })
    }
    else if (side == "both") {
      pairpoints <- matrix(0, sum(1:(length(points) - 1)), 2)
      index.pairpoints <- 0
      for (i in 1:(length(points) - 1)) {
        for (j in (i + 1):length(points)) {
          index.pairpoints <- index.pairpoints + 1
          pairpoints[index.pairpoints, ] <- c(points[i], points[j])
        }
      }
      specificities <- sapply(1:index.pairpoints, function(i) {
        round(ecdf(controls)(pairpoints[i, 2]) - ecdf(controls)(pairpoints[i, 1]), 5)
      })
      sensitivities <- sapply(1:index.pairpoints, function(i) {
        round(ecdf(cases)(pairpoints[i, 1]) + 1 - ecdf(cases)(pairpoints[i, 2]), 5)
      })
    }
    else if (side == "both2") {
      stop("side = 'both2' can be considered just if Ni = NULL.")
    }
    FPR <- c(0, sort(1 - specificities), 1)
    TPR <- c(0, sensitivities[order(1 - specificities)], 1)
    sFPR <- unique(FPR[which(duplicated(TPR) == 0)])
    sTPR <- sapply(sFPR, function(x) {
      max(TPR[which(abs(FPR - x) < sqrt(.Machine$double.eps))])
    })
    coordinates <- cbind(sFPR, sTPR)
    colnames(coordinates) <- c("FPR", "TPR")
    if (side != "both") {
      points.coordinates.points <- sapply(1:nrow(coordinates), function(i) {
        index <- which((abs(1 - specificities - sFPR[i]) < sqrt(.Machine$double.eps)) & (abs(sensitivities - sTPR[i]) < sqrt(.Machine$double.eps)))
        points[index]
      })
      points.coordinates <- cbind(points.coordinates.points, coordinates)
      colnames(points.coordinates) <- c("c", "FPR", "TPR")
    }
    else {
      index <- unlist(lapply(1:nrow(coordinates), function(i) {
        index <- which((abs(1 - specificities - sFPR[i]) < sqrt(.Machine$double.eps)) & (abs(sensitivities - sTPR[i]) < sqrt(.Machine$double.eps)))
      }))
      pairpoints.coordinates.points <- matrix(pairpoints[index,], nrow = length(index), ncol = 2)
      pairpoints.coordinates.coor <- matrix(c(1 - specificities[index], sensitivities[index]), ncol = 2)
      pairpoints.coordinates <- cbind(pairpoints.coordinates.points,  pairpoints.coordinates.coor)
      colnames(pairpoints.coordinates) <- c("xl", "xu", "FPR", "TPR")
    }
    area <- mean(ROC.t[-1] + ROC.t[-Ni - 1])/2
    if (side != "both") {
      results <- list(levels = levels, controls = controls, cases = cases, side = side,
                      pvalue.wilcox = pvalue.wilcox, t = grid, roc = ROC.t, auc = area,
                      points.coordinates = points.coordinates, param = FALSE, Ni = Ni,
                      specificities = specificities, sensitivities = sensitivities,
                      points = points, coordinates = coordinates)
    }
    else {
      results <- list(levels = levels, controls = controls, cases = cases, side = side,
                      pvalue.wilcox = pvalue.wilcox, t = grid, roc = ROC.t, auc = area,
                      pairpoints.coordinates = pairpoints.coordinates, param = FALSE, Ni = Ni,
                      specificities = specificities, sensitivities = sensitivities,
                      pairpoints = pairpoints, coordinates = coordinates, index = index)
    }

    attr(results, "class") <- "groc"

  }

  if(plot.roc){
    if(!plot.density){par(mfrow = c(1, 1))}
    plot(results, cex.lab = 1.2)
  }

  return(results)

}
