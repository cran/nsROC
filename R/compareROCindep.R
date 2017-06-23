compareROCindep <- function(X, G, D, ...) {
  UseMethod("compareROCindep")
}

compareROCindep.default <- function(X, G, D, statistic=c("L1","L2","CR","other","VK","AUC"), FUN.stat.int=function(roc.i, roc){mean(abs(roc.i - roc))}, FUN.stat.cons=function(n.cases, n.controls){sqrt(n.cases)}, side=c("right","left"), Ni=1000, raw=FALSE, perm=500, seed=123, plot.roc=TRUE, type='s', lwd=3, lwd.curves=rep(2,length(table(G))), lty=1, lty.curves=rep(1,length(table(G))), col='black', col.curves=rainbow(length(table(G))), cex.lab=1.2, legend=c(sapply(1:length(table(G)),function(i){eval(bquote(expression(hat(R)[.(i)](t))))}), expression(hat(R)(t))), legend.position='bottomright', legend.inset=0.03, cex.legend=1, ...){

  # Check if input arguments are correct

  G <- as.factor(G)
  D <- as.factor(D)

  levels <- base::levels(D)

  if( missing(X) | missing(D) | missing(G) | is.null(X) | is.null(D) | is.null(G) | !is.numeric(X) | sum(!is.na(X))==0 | sum(!is.na(D))==0 | sum(!is.na(G))==0){
    if( missing(X) | is.null(X) | sum(!is.na(X))==0 ){stop("Marker values should be provided.")}
    else if( missing(D) | is.null(D) | sum(!is.na(D))==0 ){stop("Response values should be provided.")}
    else if( missing(G) | is.null(G) | sum(!is.na(G))==0 ){stop("Group values should be provided.")}
    else if( !is.numeric(X) ){stop("Marker values are not numeric.")}
  }
  else{
    if( sum(is.na(X))!=0 ){
      D <- D[!is.na(X)]
      G <- G[!is.na(X)]
      X <- X[!is.na(X)]
      warning("Some of marker values are NA. They have been removed.")
    }
    if( sum(is.na(D))!=0 ){
      X <- X[!is.na(D)]
      G <- G[!is.na(D)]
      D <- D[!is.na(D)]
      warning("Some of response values are NA. They have been removed.")
    }
    if( sum(is.na(G))!=0 ){
      X <- X[!is.na(G)]
      D <- D[!is.na(G)]
      G <- G[!is.na(G)]
      warning("Some of group values are NA. They have been removed.")
    }
    if( length(X)!=length(D)){
      stop("Marker and response vectors should have the same length.")
    }
    if( length(X)!=length(G)){
      stop("Marker and group vectors should have the same length.")
    }

    if( length(levels)<2 ){
      stop("Response vector must contain at least two different values.")
    }
    if( length(levels)>2 ){
      warning("There are more than two levels in response vector. Only two are considered.")
    }
  }

  statistic <- match.arg(statistic)
  side <- match.arg(side)
  if(statistic=="VK" && side=="left"){stop("Venkatraman method is implemented only for right side ROC curves.\n")}

  if(!is.numeric(perm) || length(perm)!=1 || perm%%1!=0 || perm <= 0){
    stop("perm should be a positive integer.")
  }


  # Calculate sizes and final data considered

  k <- length(base::levels(G))
  base::levels(G) <- 1:k
  samples.check <- lapply(1:k, function(i){
      tempX <- split(X,G)[[i]]
      tempD <- split(D,G)[[i]]
      controls <- split(tempX,tempD)[[levels[1]]]
      cases <- split(tempX,tempD)[[levels[2]]]
      n0 <- length(controls); n1 <- length(cases)
      list(controls=controls, cases=cases, n0=n0, n1=n1)
    })

  n.controls <- sapply(1:k, function(i){samples.check[[i]]$n0})
  n.cases <- sapply(1:k, function(i){samples.check[[i]]$n1})
  controls.k <- as.numeric(unlist(sapply(1:k, function(i){samples.check[[i]]$controls})))
  cases.k <- as.numeric(unlist(sapply(1:k, function(i){samples.check[[i]]$cases})))


  # Check the rest of input arguments

  if(length(lwd.curves)!=k){
    stop("lwd.curves should have the same length as the number of samples to compare.")
  }

  if(length(lty.curves)!=k){
    stop("lty.curves should have the same length as the number of samples to compare.")
  }

  if(length(col.curves)!=k){
    stop("col.curves should have the same length as the number of samples to compare.")
  }

  if(!is.numeric(Ni) || length(Ni)!=1 || Ni%%1!=0 || Ni <= 0){
    stop("Ni should be a positive integer.")
  }else if(Ni < max(n.controls)){
    warning("It is advisable to consider Ni higher than maximum number of controls.")
  }

  cat("The ", k, " ROC curves considered to compare come from these data: \n", sep="")
  sapply(1:k, function(i){
    cat("In the sample", i, "there are", n.controls[i], "controls and", n.cases[i], "cases.\n")
  })
  cat('\n')

  # ROC curves estimates

  if(raw){
  	controls.vector <- controls.k
  	cases.vector <- cases.k
  }else{
  	x.vector <- lapply(1:k, function(i){rank(c(samples.check[[i]]$controls, samples.check[[i]]$cases),ties.method='first')})
  	controls.vector <- as.numeric(unlist(sapply(1:k, function(i){x.vector[[i]][1:n.controls[i]]})))
  	cases.vector <- as.numeric(unlist(sapply(1:k, function(i){x.vector[[i]][(n.controls[i]+1):(n.controls[i]+n.cases[i])]})))
  }

  grid <- seq(0,1,1/Ni)
  if(statistic=="CR"){grid <- seq(0,1,1/(2*Ni))}
  roc.ord <- function(side,controls,cases){
      switch(side,
  			right = ecdf(1-ecdf(controls)(cases))(grid),
  			left = ecdf(ecdf(controls)(cases))(grid)
  		)
    }

  controls.i <- split(controls.vector,rep(1:k, n.controls))
  cases.i <- split(cases.vector,rep(1:k, n.cases))

  ROC.t <- sapply(1:k, function(i){roc.ord(side=side, controls=controls.i[[i]], cases=cases.i[[i]])})
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



  if(statistic!="VK" && statistic!="AUC"){   # statistic = "L1","L2","CR","other"

    if(statistic!="other"){   # statistic = "L1","L2","CR"

      stat.name <- function(statistic){
          switch(statistic,
              L1 = cat("Statistic considered: L1 measure\n"),
              L2 = cat("Statistic considered: L2 measure\n"),
              CR = cat("Statistic considered: Cramer-von Mises\n")
          )
        }
      stat.name(statistic)

      statistic.i <- function(statistic, roc.i, roc){
          switch(statistic,
    				  L1 = mean(abs(roc.i - roc)),
    				  L2 = mean((roc.i - roc)^2),
    				  CR = mean((roc.i[seq(2,2*Ni+1,2)] - roc[seq(2,2*Ni+1,2)])^2 * (roc[seq(3,2*Ni+1,2)] - roc[seq(1,2*Ni-1,2)]))
    			)
    		}
      statistic.integ <- sapply(1:k, function(i){statistic.i(statistic, ROC.t[,i], ROC)})
    	statistic.cons <- n.cases
    	if(statistic=="L1"){statistic.cons <- sqrt(n.cases)}

    }else{   # statistic = "other"

      statistic.integ <- sapply(1:k, function(i){FUN.stat.int(ROC.t[,i], ROC)})
      statistic.cons <- FUN.stat.cons(n.cases, n.controls)

    }

    stat <- sum(statistic.cons*statistic.integ)

    cat('\nProgress bar: Estimation of statistic value in permutation iterations (perm = ', perm,')\n',sep='')
    bar <- txtProgressBar(min = 0, max = perm, style = 3)
    set.seed(seed)
    stat.perm <- sapply(1:perm, function(b){
      setTxtProgressBar(bar, b)
      controls.perm <- split(controls.vector, sample(rep(1:k, n.controls)))
      cases.perm <- split(cases.vector, sample(rep(1:k, n.cases)))

      X.perm.breakties <- lapply(1:k, function(i){
        X.perm <- c(controls.perm[[i]], cases.perm[[i]])
        if(raw){
          X.perm
        }else{
          rank(X.perm, ties.method='random')
        }
      })

      controls.perm.bt <- lapply(1:k, function(i){X.perm.breakties[[i]][1:n.controls[i]]})
      cases.perm.bt <- lapply(1:k, function(i){X.perm.breakties[[i]][(n.controls[i]+1):(n.controls[i]+n.cases[i])]})

      ROC.t.perm <- sapply(1:k, function(i){roc.ord(side=side, controls=controls.perm.bt[[i]], cases=cases.perm.bt[[i]])})
      ROC.perm <- apply(ROC.t.perm, 1, mean)
      if(statistic!="other"){
        sum(statistic.cons * sapply(1:k, function(i){
          statistic.i(statistic, ROC.t.perm[,i], ROC.perm)
        }))
      }else{
        sum(statistic.cons * sapply(1:k, function(i){
          FUN.stat.int(ROC.t.perm[,i], ROC.perm)
        }))
      }
    })
    close(bar)

  }else{   # statistic = "VK","AUC"

      if(statistic=="VK"){   # statistic = "VK"

        cat("Statistic considered: Venkatraman (Venkatraman, 2000)\n")

        stat.i.j <- lapply(1:(k-1), function(i){sapply((i+1):k,function(j){
          xi <- sort(c(controls.i[[i]], cases.i[[i]])); xj <- sort(c(controls.i[[j]], cases.i[[j]]))
          f.xi <- ecdf(cases.i[[i]])(xi); f.xj <- ecdf(cases.i[[j]])(xj)
          g.xi <- ecdf(controls.i[[i]])(xi); g.xj <- ecdf(controls.i[[j]])(xj)
          kappa <- (n.cases[i]+n.cases[j])/(n.cases[i]+n.controls[i]+n.cases[j]+n.controls[j])
          p.xi.i <- kappa*f.xi + (1-kappa)*g.xi; p.xj.i <- kappa*f.xj + (1-kappa)*g.xj
          e.xi.i <- kappa*f.xi + (1-kappa)*(1-g.xi); e.xj.i <- kappa*f.xj + (1-kappa)*(1-g.xj)
          e.xi <- approxfun(c(0,p.xi.i,1),c(1-kappa,e.xi.i,kappa)); e.xj <- approxfun(c(0,p.xj.i,1),c(1-kappa,e.xj.i,kappa))
          integrate(function(p){abs(e.xi(p)-e.xj(p))}, 0, 1, subdivisions = 1000)$value
        })})

        stat <- sum(unlist(stat.i.j))	# If k>2 the statistic value is the sum of statistic values of each pair

        cat('\nProgress bar: Estimation of statistic value in permutation iterations (perm = ', perm,')\n',sep='')
        bar <- txtProgressBar(min = 0, max = perm, style = 3)
        set.seed(seed)
        stat.perm <- sapply(1:perm, function(b){

          setTxtProgressBar(bar, b)

          controls.perm <- split(controls.vector, sample(rep(1:k, n.controls)))
          cases.perm <- split(cases.vector, sample(rep(1:k, n.cases)))

          X.perm.breakties <- lapply(1:k, function(i){
            X.perm <- c(controls.perm[[i]], cases.perm[[i]])
            if(raw){
              X.perm
            }else{
              rank(X.perm, ties.method='random')
            }
          })

          controls.perm.bt <- lapply(1:k, function(i){X.perm.breakties[[i]][1:n.controls[i]]})
          cases.perm.bt <- lapply(1:k, function(i){X.perm.breakties[[i]][(n.controls[i]+1):(n.controls[i]+n.cases[i])]})

          stat.i.j.b <- lapply(1:(k-1), function(i){sapply((i+1):k,function(j){
            xi <- sort(X.perm.breakties[[i]]); xj <- sort(X.perm.breakties[[j]])
            f.xi <- ecdf(cases.perm.bt[[i]])(xi); f.xj <- ecdf(cases.perm.bt[[j]])(xj)
            g.xi <- ecdf(controls.perm.bt[[i]])(xi); g.xj <- ecdf(controls.perm.bt[[j]])(xj)
            kappa <- (n.cases[i]+n.cases[j])/(n.cases[i]+n.controls[i]+n.cases[j]+n.controls[j])
            p.xi.i <- kappa*f.xi + (1-kappa)*g.xi; p.xj.i <- kappa*f.xj + (1-kappa)*g.xj
            e.xi.i <- kappa*f.xi + (1-kappa)*(1-g.xi); e.xj.i <- kappa*f.xj + (1-kappa)*(1-g.xj)
            e.xi <- approxfun(c(0,p.xi.i,1),c(1-kappa,e.xi.i,kappa)); e.xj <- approxfun(c(0,p.xj.i,1),c(1-kappa,e.xj.i,kappa))
            integrate(function(p){abs(e.xi(p)-e.xj(p))}, 0, 1, subdivisions = 1000)$value
          })})

          sum(unlist(stat.i.j.b))		# If k>2 the statistic value is the sum of statistic values of each pair

        })
        close(bar)

      }else{   # statistic = "AUC"

        cat("Statistic considered: AUC\n")

        area.i <- sapply(1:k, function(i){
            X.i <- c(controls.i[[i]], cases.i[[i]])
            D.i <- c(rep(levels[1], n.controls[i]), rep(levels[2], n.cases[i]))
            gROC(X.i, D.i, side=side, Ni=Ni)$area
          })
        stat <- mean(abs(area.i - mean(area.i)))

        cat('\nProgress bar: Estimation of statistic value in permutation iterations (perm = ', perm,')\n',sep='')
        bar <- txtProgressBar(min = 0, max = perm, style = 3)
        set.seed(seed)
        stat.perm <- sapply(1:perm, function(b){
          setTxtProgressBar(bar, b)
          controls.perm <- split(controls.vector, sample(rep(1:k, n.controls)))
          cases.perm <- split(cases.vector, sample(rep(1:k, n.cases)))

          X.perm.breakties <- lapply(1:k, function(i){
            X.perm <- c(controls.perm[[i]], cases.perm[[i]])
            if(raw){
              X.perm
            }else{
              rank(X.perm, ties.method='random')
            }
          })

          controls.perm.bt <- lapply(1:k, function(i){X.perm.breakties[[i]][1:n.controls[i]]})
          cases.perm.bt <- lapply(1:k, function(i){X.perm.breakties[[i]][(n.controls[i]+1):(n.controls[i]+n.cases[i])]})

          area.i.perm <- sapply(1:k, function(i){
            D.i <- c(rep(levels[1], n.controls[i]), rep(levels[2], n.cases[i]))
            gROC(X.perm.breakties[[i]], D.i, side=side, Ni=Ni)$area
          })
          mean(abs(area.i.perm - mean(area.i.perm)))
        })
        close(bar)

      }

  }

  p.value <- mean(stat.perm > stat)

  cat("\nTest output:\n")
  cat("Null hypothesis: The", k, "considered ROC curves (independent) are equal.\n")
  cat("Statistic value = ", stat,'\t\t',"p-value = ", p.value,'\n')

  list(n.controls=n.controls, n.cases=n.cases, controls.k=controls.k, cases.k=cases.k, statistic=stat, stat.perm=stat.perm, p.value=p.value)

}
