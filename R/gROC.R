gROC <- function(X, D, ...) {
  UseMethod("gROC")
}

gROC.default <- function(X, D, side=c("auto","right","left","both"), Ni=1000, plot.roc=FALSE, plot.density=FALSE, ...){

	data <- checkROC(X, D)
  levels <- data$levels;	controls <- data$controls; cases <- data$cases; n0 <- data$n0; n1 <- data$n1; X <- data$X; D <- data$D

	# Check if input arguments are correct

		side <- match.arg(side)

		if(!is.numeric(Ni) || length(Ni)!=1 || Ni%%1!=0){
		  stop("Ni should be an integer.")
		}else if(Ni <= 0){
		  stop("Ni should be a positive number.")
		}else if(Ni < n0){
		    warning("It is advisable to consider Ni higher than number of controls.")
		}

	# Plot empirical density functions of controls and cases

		if(plot.density){
			densityControls <- density(controls)
			densityCases <- density(cases)
			minDensMarker <- min(c(densityControls$x,densityCases$x)); maxDensMarker <- max(c(densityControls$x,densityCases$x))
			minDens <- min(c(densityControls$y,densityCases$y)); maxDens <- max(c(densityControls$y,densityCases$y))
			if(plot.roc){
				par(mfrow=c(1,2))
			}else{
			  par(mfrow=c(1,1))
			}
			plot(densityControls, col=3, lwd=2, xlim=c(minDensMarker,maxDensMarker), ylim=c(minDens,maxDens), xlab="Marker values", ylab="", cex.lab=1.2, main="Density estimation")
			lines(densityCases, col=2, lwd=2)
			legend("topright", legend=c("Controls","Cases"), inset=0.03, lwd=c(2,2), lty=c(1,1), col=c(3,2))
		}

	# Calculate cut-off points considered

		# points <- (c(-Inf,unique(sort(X))) + c(unique(sort(X)),Inf))/2
		points <- c(-Inf,unique(sort(X)),Inf)

	# Choose the best option (left or right-side) to estimate the ROC curve in terms of median (only if side=="auto")

  pvalue.wilcox <- NULL
	if(side=="auto"){
	  wilcox <- wilcox.test(controls,cases,"less")
		if(wilcox$statistic <= n0*n1/2){
			side <- "right"
      pvalue.wilcox <- wilcox$p.value
		}
		else{
			side <- "left"
			pvalue.wilcox <- wilcox$p.value
		}

	}

	# Calculate the values of R(t) for points t in (0:Ni)/Ni

		grid <- seq(0,1,1/Ni)
		gamma <- seq(0,1,0.001)
		sol <- function(side,controls,cases){
			switch(side,
				right = ecdf(1-ecdf(controls)(cases))(grid),
				left = ecdf(ecdf(controls)(cases))(grid),
				both = apply(matrix(1-ecdf(1-ecdf(controls)(cases))(1-gamma%*%t(grid)) + ecdf(1-ecdf(controls)(cases))((1-gamma)%*%t(grid)),length(gamma),length(grid)), 2, max) )
			}
		ROC.t <- sol(side=side, controls=controls, cases=cases)

	# Calculate specificities and sensitivities in each cut-off point considered. In case of side="both", in each pair of cut-off points (xl,xu) such that xl<xu

		if(side=="right"){
			specificities <- sapply(points,function(c){sum(controls<=c)/length(controls)})
			sensitivities <- sapply(points,function(c){sum(cases>c)/length(cases)})
		}else if(side=="left"){
			specificities <- sapply(points,function(c){sum(controls>c)/length(controls)})
			sensitivities <- sapply(points,function(c){sum(cases<=c)/length(cases)})
		}else if(side=="both"){
			pairpoints <- matrix(0,sum(1:(length(points)-1)),2)
			index.pairpoints <- 0
			for(i in 1:(length(points)-1)){
				for(j in (i+1):length(points)){
					index.pairpoints <- index.pairpoints + 1
					pairpoints[index.pairpoints,] <- c(points[i],points[j])
				}
			}
			specificities <- sapply(1:index.pairpoints, function(i){round(ecdf(controls)(pairpoints[i,2]) - ecdf(controls)(pairpoints[i,1]), 5)})
			sensitivities <- sapply(1:index.pairpoints, function(i){round(ecdf(cases)(pairpoints[i,1]) + 1 - ecdf(cases)(pairpoints[i,2]), 5)})
  	}

	# Coordinates of points where the estimation of ROC curve has a step

		FPR <- c(0,sort(1-specificities),1)
		TPR <- c(0,sensitivities[order(1-specificities)],1)
		sFPR <- unique(FPR[which(duplicated(TPR)==0)])
		sTPR <- sapply(sFPR, function(x){max(TPR[which(abs(FPR-x) < sqrt(.Machine$double.eps))])})
		coordinates <- cbind(sFPR,sTPR)
		colnames(coordinates) <- c("FPR","TPR")

		if(side!="both"){
			points.coordinates.points <- sapply(1:nrow(coordinates), function(i){
					index <- which((abs(1-specificities-sFPR[i]) < sqrt(.Machine$double.eps)) & (abs(sensitivities-sTPR[i]) < sqrt(.Machine$double.eps)))
					points[index]})
			points.coordinates <- cbind(points.coordinates.points, coordinates)
			colnames(points.coordinates) <- c("c","FPR","TPR")
		}else{
			index <- unlist(lapply(1:nrow(coordinates), function(i){
					index <- which((abs(1-specificities-sFPR[i]) < sqrt(.Machine$double.eps)) & (abs(sensitivities-sTPR[i]) < sqrt(.Machine$double.eps)))
					}))
			pairpoints.coordinates.points <- matrix(pairpoints[index,], nrow=length(index), ncol=2)
			pairpoints.coordinates.coor <- matrix(c(1-specificities[index],sensitivities[index]), ncol=2)
			pairpoints.coordinates <- cbind(pairpoints.coordinates.points, pairpoints.coordinates.coor)
      colnames(pairpoints.coordinates) <- c("c1","c2","FPR","TPR")
		}

	# Calculate the area under the ROC curve (AUC)

		area <- mean(ROC.t[-1] + ROC.t[-Ni-1])/2


	# Display results

		if(side!="both"){
			results <- list(levels=levels, controls=controls, cases=cases, side=side, pvalue.wilcox=pvalue.wilcox, specificities=specificities, sensitivities=sensitivities, points=points, area=area, coordinates=coordinates, points.coordinates=points.coordinates, Ni=Ni, ROC.t=ROC.t)
		}else{
			results <- list(levels=levels, controls=controls, cases=cases, side=side, pvalue.wilcox=pvalue.wilcox, specificities=specificities, sensitivities=sensitivities, pairpoints=pairpoints, area=area, coordinates=coordinates, index=index, pairpoints.coordinates=pairpoints.coordinates, Ni=Ni, ROC.t=ROC.t)
		}

	attr(results, 'class') <- 'groc'

	# Plot ROC curve estimation

  if(plot.roc){
    if(!plot.density){
      par(mfrow=c(1,1))
    }
    plot(results, cex.lab=1.2)
  }

	return(results)

}
