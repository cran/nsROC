compareROCdep <- function(X, D, ...) {
  UseMethod("compareROCdep")
}

compareROCdep.default <- function(X, D, method=c("general.bootstrap","permutation","auc"), statistic=c("KS","L1","L2","CR","VK","other"), FUN.dist=function(g){max(abs(g))}, side=c("right","left"), Ni=1000, B=500, perm=500, seed=123, h.fun=function(H,x){H*sd(x)*length(x)^{-1/3}}, H=1, plot.roc=TRUE, type='s', lwd=3, lwd.curves=rep(2,ncol(X)), lty=1, lty.curves=rep(1,ncol(X)), col='black', col.curves=rainbow(ncol(X)), cex.lab=1.2, legend=c(sapply(1:ncol(X),function(i){eval(bquote(expression(hat(R)[.(i)](t))))}), expression(hat(R)(t))), legend.position='bottomright', legend.inset=0.03, cex.legend=1, ...){

  # Check if input arguments are correct

  method <- match.arg(method)
  statistic <- match.arg(statistic)
  side <- match.arg(side)
  if(statistic=="VK" && side=="left"){stop("Venkatraman method is implemented only for right side ROC curves.\n")}
  if(statistic=="VK"){method <- "permutation"}

  if(!is.numeric(B) || length(B)!=1 || B%%1!=0 || B <= 0){
    stop("B should be a positive integer.")
  }

  if(!is.numeric(perm) || length(perm)!=1 || perm%%1!=0 || perm <= 0){
    stop("perm should be a positive integer.")
  }

  if(!is.numeric(H) || length(H)!=1 || H < 0){
    stop("H should be a non-negative number.")
  }

  marker.samples <- lapply(1:ncol(X), function(i) X[,i])

  k <- length(marker.samples)

  if(length(lwd.curves)!=k){
    stop("lwd.curves should have the same length as the number of samples.")
  }

  if(length(lty.curves)!=k){
    stop("lty.curves should have the same length as the number of samples.")
  }

  if(length(col.curves)!=k){
    stop("col.curves should have the same length as the number of samples.")
  }

  sample.size <- unique(c(sapply(1:k, function(i){length(marker.samples[[i]])}), length(D)))
  index.na <- c(unlist(lapply(1:k, function(i){which(is.na(marker.samples[[i]]))})), which(is.na(D)))

  if(length(sample.size)!=1){
    stop("Response and all marker vectors should have the same length.")
  }else{
    marker.samples.matrix <- matrix(unlist(marker.samples), sample.size, k)
    if(length(index.na)!=0L){
      warning("Those subjects with any missing marker or response value have been removed.")
      index.na <- unique(index.na)
      marker.samples.matrix <- marker.samples.matrix[-index.na,]
      D <- D[-index.na]
    }
  }

  samples.check <- lapply(1:k, function(i){
			checkROC(marker.samples.matrix[,i],D)
			})

  n.controls <- samples.check[[1]]$n0
  n.cases <- samples.check[[1]]$n1

  if(!is.numeric(Ni) || length(Ni)!=1 || Ni%%1!=0 || Ni <= 0){
    stop("Ni should be a positive integer.")
  }else if(Ni < n.controls){
    warning("It is advisable to consider Ni higher than the number of controls.")
  }

  # ROC curves estimates

  controls.k <- sapply(1:k, function(i){samples.check[[i]]$controls})
  cases.k <- sapply(1:k, function(i){samples.check[[i]]$cases})

  cat("In the considered database there are", n.controls, "controls and", n.cases, "cases.\n")
  cat('\n')

  grid <- seq(0,1,1/Ni)
  roc.ord <- function(side,controls,cases){
		switch(side,
			right = ecdf(1-ecdf(controls)(cases))(grid),
			left = ecdf(ecdf(controls)(cases))(grid)
		  )
    }

  ROC.t <- sapply(1:k, function(i){roc.ord(side=side, controls=controls.k[,i], cases=cases.k[,i])})
  ROC <- apply(ROC.t, 1, mean)

  if(plot.roc){
    par(mfrow=c(1,1))
    plot(c(0,grid,1), c(0,ROC,1), type=type, lwd=lwd, lty=lty, xlab="False-Positive Rate", ylab="True-Positive Rate", cex.lab=cex.lab, main="ROC curves", col=col, ...)
    axis(1, at=seq(0,1,0.01), labels=F, tck=-0.01)
    axis(1, at=seq(0.1,0.9,0.1), labels=F, tck=-0.02)
    axis(2, at=seq(0,1,0.01), labels=F, tck=-0.01)
    axis(2, at=seq(0.1,0.9,0.1), labels=F, tck=-0.02)
    abline(0, 1, lty=1, col="gray")

    for(i in 1:k){lines(c(0,grid,1), c(0,ROC.t[,i],1), type=type, lwd=lwd.curves[i], lty=lty.curves[i], col=col.curves[i])}
    lines(c(0,grid,1), c(0,ROC,1), type=type, lwd=lwd, lty=lty, col=col)

    legend(legend.position, inset=legend.inset, legend=legend, lwd=c(lwd.curves,lwd), lty=c(lty.curves,lty), col=c(col.curves, col=col), cex=cex.legend)
  }


  if(method!="auc"){	# method = "general.bootstrap","permutation"

    if(statistic!="other"){		# statistic = "KS","L1","L2","CR","VK"

      if(statistic!="CR" && statistic!="VK"){		# statistic = "KS","L1","L2"

        stat.name <- function(statistic){
						switch(statistic,
							KS = cat("Statistic considered: Kolmogorov-Smirnov\n"),
							L1 = cat("Statistic considered: L1 measure\n"),
							L2 = cat("Statistic considered: L2 measure\n"))
          }
        stat.name(statistic)

        FUN.dist.stat <- function(statistic,g){
						switch(statistic,
							KS = max(abs(g)),
							L1 = mean(abs(g)),
							L2 = mean(g^2))
          }
        FUN.dist <- function(g){FUN.dist.stat(statistic=statistic, g)}

        stat <- sum(sapply(1:k,function(i){FUN.dist(g=sqrt(n.cases)*(ROC.t[,i] - ROC))}))

      }else{

        if(statistic=="CR"){		# statistic = "CR"

          cat("Statistic considered: Cramer-Von Mises\n")

          FUN.dist.CR <- function(g,h){sum(g[-length(g)]^2*(h[-1]-h[-length(h)]))}
          stat <- sum(sapply(1:k,function(i){FUN.dist.CR(g=sqrt(n.cases)*(ROC.t[,i] - ROC), h=ROC)}))

        }else{		# statistic = "VK"

          cat("Statistic considered: Venkatraman\n")

          X.k <- rbind(controls.k, cases.k)
          X.rank <- apply(X.k, 2, rank, ties.method='first')
          controls.rank <- X.rank[1:n.controls,]
          cases.rank <- X.rank[(n.controls+1):(n.controls+n.cases),]

          stat.i.j <- lapply(1:(k-1), function(i){sapply((i+1):k,function(j){
            total.errors.x <- c(sapply(1:(n.controls+n.cases-1), function(l){sum(cases.rank[,i]<=l) + sum(controls.rank[,i]>l)}), n.cases)
            total.errors.y <- c(sapply(1:(n.controls+n.cases-1), function(l){sum(cases.rank[,j]<=l) + sum(controls.rank[,j]>l)}), n.cases)
            E0 <- sum(abs(total.errors.y - total.errors.x))
            })})
          stat <- sum(unlist(stat.i.j))
			  }
		  }

  	}else{	# statistic = "other"

  		stat <- sum(sapply(1:k,function(i){FUN.dist(sqrt(n.cases)*(ROC.t[,i] - ROC))}))

  	}

  	if(method=="general.bootstrap"){	# method = "general.bootstrap"

  		cat("Method considered: General bootstrap (Martinez-Camblor and Corral, 2012)\n")

  		boot.cons <- function(H,X){
  		  sapply(1:k, function(i){
  		    h <- h.fun(H,X[,i])
  		    rnorm(nrow(X),0,h)
  		    })
  		  }

  		cat('\nProgress bar: Estimation of ROC curve in bootstrap iterations (B = ', B,')\n',sep='')
  		bar1 <- txtProgressBar(min = 0, max = B, style = 3)
  		set.seed(seed)
  		ROC.i.b <- lapply(1:B, function(b){
  		  setTxtProgressBar(bar1, b)
  		  controls.b <- controls.k[sample(1:n.controls, replace=TRUE),] + boot.cons(H=H, X=controls.k)
  		  cases.b <- cases.k[sample(1:n.cases, replace=TRUE),] + boot.cons(H=H, X=cases.k)
  		  sapply(1:k, function(i){roc.ord(side=side, controls=controls.b[,i], cases=cases.b[,i])})
  		  })
  		close(bar1)

  		ROC.i <- sapply(1:k, function(i){
  		  ROC.i.boot <- sapply(1:B, function(b){ROC.i.b[[b]][,i]})
  		  apply(ROC.i.boot, 1, mean)
  		  })

  		cat('Progress bar: Estimation of statistic value in bootstrap iterations (B = ',B,')\n',sep='')
  		bar2 <- txtProgressBar(min = 0, max = B, style = 3)
  		stat.boot <- sapply(1:B,function(b){
  		    setTxtProgressBar(bar2, b)
  		    if(statistic!="CR"){	# statistic = "KS","L1","L2","other"

  		      sum(sapply(1:k,function(i){
  		        g.i.j <- sapply(1:k, function(j){((i==j)-1/k)*sqrt(n.cases)*(ROC.i.b[[b]][,j] - ROC.i[,j])})
  		        FUN.dist(apply(g.i.j,1,sum))}))

  		    }else{	# statistic = "CR"

  		      ROC.b <- apply(ROC.i.b[[b]], 1, mean)
  	  	    sum(sapply(1:k,function(i){FUN.dist.CR(g=sqrt(n.cases)*(ROC.i.b[[b]][,i] - ROC.t[,i]), h=ROC.b)})) - k*(FUN.dist.CR(g=sqrt(n.cases)*(ROC.b - ROC), h=ROC.b))

    		  }
  		  })
  		close(bar2)

  		p.value <- mean(stat.boot >= stat)

  		cat("\nTest output:\n")
  		cat("Null hypothesis: The ", k, "considered ROC curves (paired) are equal.\n")
  		cat("Statistic value = ", stat,'\t\t',"p-value = ", p.value,'\n')

  		list(n.controls=n.controls, n.cases=n.cases, controls.k=controls.k, cases.k=cases.k, statistic=stat, stat.boot=stat.boot, p.value=p.value)

  	}else{	# method = "permutation"

  		cat("Method considered: Permutation (Venkatraman and Begg, 1996)\n")

  		cat('\nProgress bar: Estimation of statistic value in permutation iterations (perm = ', perm,')\n',sep='')
  		bar1 <- txtProgressBar(min = 0, max = perm, style = 3)
  		set.seed(seed)
  		stat.perm <- sapply(1:perm, function(b){

  		    setTxtProgressBar(bar1, b)
  				X.k <- rbind(controls.k, cases.k)
  				X.rank <- apply(X.k, 2, rank, ties.method='first')
  				X.perm <- apply(t(apply(X.rank, 1, sample)), 2, rank, ties.method='random')
  				controls.b <- X.perm[1:n.controls,]
  				cases.b <- X.perm[(n.controls+1):(n.controls+n.cases),]
  				ROC.i.b <- sapply(1:k, function(i){roc.ord(side=side, controls=controls.b[,i], cases=cases.b[,i])})
  				ROC.b <- apply(ROC.i.b, 1, mean)

  				if(statistic!="CR" && statistic!="VK"){	# statistic = "KS","L1","L2","other"

  				  sum(sapply(1:k,function(i){FUN.dist(sqrt(n.cases)*(ROC.i.b[,i] - ROC.b))}))

  				}else{

  				  if(statistic=="CR"){	# statistic = "CR"

  				    sum(sapply(1:k,function(i){FUN.dist.CR(g=sqrt(n.cases)*(ROC.i.b[,i] - ROC.b), h=ROC.b)}))

  				  }else{	# statistic = "VK"

  						stat.i.j.b <- lapply(1:(k-1), function(i){sapply((i+1):k,function(j){
  						    total.errors.x <- c(sapply(1:(n.controls+n.cases-1), function(l){sum(cases.b[,i]<=l) + sum(controls.b[,i]>l)}), n.cases)
  						    total.errors.y <- c(sapply(1:(n.controls+n.cases-1), function(l){sum(cases.b[,j]<=l) + sum(controls.b[,j]>l)}), n.cases)
  						    E0 <- sum(abs(total.errors.y - total.errors.x))
  						  })})

  						sum(unlist(stat.i.j.b))

  					}
  				}
  			})
  		cat('\n')

  		p.value <- mean(stat.perm >= stat)

  		cat("\nTest output:\n")
  		cat("Null hypothesis: The ", k, "considered ROC curves (paired) are equal.\n")
  		cat("Statistic value = ", stat,'\t\t',"p-value = ", p.value,'\n')

  		list(n.controls=n.controls, n.cases=n.cases, controls.k=controls.k, cases.k=cases.k, statistic=stat, stat.perm=stat.perm, p.value=p.value)

  	}

  }else{	# method = "auc"

    cat("Method considered: AUC comparison (DeLong, DeLong and Clarke-Pearson, 1988)\n")

  	cat('\nProgress bar: Estimation of statistic value in each variable (k = ', k,')\n',sep='')
  	bar1 <- txtProgressBar(min = 0, max = k, style = 3)
  	kernel.auc <- function(v){
    	  if(side=='right'){
    	    ifelse(isTRUE(all.equal(v[1],v[2])), 1/2, ifelse(v[2]<v[1], 1, 0))
    	  }else{
    	    ifelse(isTRUE(all.equal(v[1],v[2])), 1/2, ifelse(v[2]>v[1], 1, 0))
    	  }
    	}
  	kernel.pairs <- sapply(1:k,function(i){
  	  setTxtProgressBar(bar1, i)
  	  apply(expand.grid(cases.k[,i],controls.k[,i]),1,kernel.auc)
  	  })
  	stat <- 1/(n.cases*n.controls)*apply(kernel.pairs, 2, sum)
  	cat('\n')

  	X.component <- 1/n.controls*sapply(1:k,function(c){
  	  sapply(1:n.cases,function(i){
  	    sum(kernel.pairs[seq(i,n.controls*n.cases,n.cases),c])
  	    })
  	  })
  	Y.component <- 1/n.cases*sapply(1:k,function(c){
  	  sapply(1:n.controls,function(j){
  	    sum(kernel.pairs[seq(n.cases*(j-1)+1,n.cases*j,1),c])
  	    })
  	  })
  	dif.X.comp <- t(X.component)-matrix(rep(stat,n.cases),k,n.cases)
  	dif.Y.comp <- t(Y.component)-matrix(rep(stat,n.controls),k,n.controls)
  	S10 <- 1/(n.cases-1)*dif.X.comp%*%t(dif.X.comp)
  	S01 <- 1/(n.controls-1)*dif.Y.comp%*%t(dif.Y.comp)
  	S <- 1/n.cases*S10 + 1/n.controls*S01
  	L <- matrix((diag(k) + cbind(rep(0,k), -diag(k))[,1:k])[1:(k-1),], k-1, k)
  	M <- L%*%S%*%t(L)

  	stat.test <- stat%*%t(L)%*%solve(M)%*%L%*%as.matrix(stat)
  	p.value <- round(1-pchisq(stat.test, qr(M)$rank), 4)

  	cat("\nTest output:\n")
  	cat("Null hypothesis: The ", k, "considered ROC curves (paired) are equal.\n")
  	cat("Statistic value = ", stat.test,'\t\t',"p-value = ", p.value,'\n')

  	list(n.controls=n.controls, n.cases=n.cases, controls.k=controls.k, cases.k=cases.k, statistic=stat, test.statistic=stat.test, p.value=p.value)

  }

}
