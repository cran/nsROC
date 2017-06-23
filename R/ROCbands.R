# The only required package is "sde" to simulate Brownian Bridges (BBridge function) at Jensen confidence bands estimation method.

# install.packages("sde")
# library(sde)

ROCbands <- function(groc, ...) {
  UseMethod("ROCbands")
}

ROCbands.default <- function(groc, method=c("PSN","JMS","DEK"), conf.level=0.95, B=500, bootstrap.bar=TRUE, alpha1=NULL, s=1, a.J=1/1000, b.J=1-1/1000, plot.bands=FALSE, plot.var=FALSE, ...){

  # Check if input arguments are correct

  method <- match.arg(method)

  if(!is.numeric(conf.level) || length(conf.level)!=1){
    stop("conf.level should be a number.")
  }else if(conf.level > 1 && conf.level <= 100){
    conf.level <- conf.level/100
  }
  if(conf.level < 0 || conf.level > 1){
    stop("conf.level should be in the unit interval.")
  }

  alpha <- 1-conf.level

  if( !is.null(alpha1) ){
    if(!is.numeric(alpha1) || length(alpha1)!=1){
      stop("alpha1 should be a number.")
    }else if(alpha1<0 | alpha1>1){
      stop("alpha1 should be in the unit interval.")
    }else if(alpha1>alpha){
      stop("alpha1 should be lower than 1-conf.level.")
    }
  }

  if(!is.numeric(B) || length(B)!=1 || B%%1!=0 || B <= 0){
    stop("B should be a positive integer.")
  }

  if(!is.numeric(a.J) || length(a.J)!=1 || a.J < 0 || a.J > 1){
    stop("a.J should be a number in the unit interval.")
  }

  if(!is.numeric(b.J) || length(b.J)!=1 || b.J < 0 || b.J > 1){
    stop("b.J should be a number in the unit interval.")
  }

  if(!is.numeric(s) || length(s)!=1 || s < 0){
    stop("B should be a non-negative number.")
  }

  # ROC curve estimate

  side <- groc$side; controls <- groc$controls; cases <- groc$cases
  ROC.t <- groc$ROC.t; Ni <- groc$Ni

  grid <- seq(0,1,1/Ni)
  gamma <- seq(0,1,0.001)
  sol <- function(side,controls,cases){
		switch(side,
			right = ecdf(1-ecdf(controls)(cases))(grid),
			left = ecdf(ecdf(controls)(cases))(grid),
			both = apply(matrix(1-ecdf(1-ecdf(controls)(cases))(1-gamma%*%t(grid)) + ecdf(1-ecdf(controls)(cases))((1-gamma)%*%t(grid)),length(gamma),length(grid)), 2, max) )
		}


		# PSN confidence bands

		if(method=="PSN"){
			n <- length(cases); m <- length(controls)
			hn <- s*min(n,m)^{-1/5}*sd(cases); hm <- s*min(n,m)^{-1/5}*sd(controls)
			cases.b <- sapply(1:B, function(b){sample(cases, replace=TRUE) + rnorm(n,0,hn)})
			controls.b <- sapply(1:B, function(b){sample(controls, replace=TRUE) + rnorm(m,0,hm)})
			if(side=="both" && bootstrap.bar==TRUE){
				cat('Progress bar of bootstrap iterations (B = ', B,")\n",sep="")
				bar <- txtProgressBar(min = 0, max = B, style = 3)
			}
			ROC.B <- sapply(1:B, function(b){
							if(side=="both" && bootstrap.bar==TRUE){
								setTxtProgressBar(bar, b)
							}
							sol(side=side, controls=controls.b[,b], cases=cases.b[,b])
							})
			if(side=="both" && bootstrap.bar==TRUE){close(bar)}
			X.B <- ROC.B - ROC.t
			sigma.B <- apply(sqrt(n)*X.B,1,sd)
			sigma.B[sigma.B==0] <- .Machine$double.eps
			X.B[sigma.B==0,] <- .Machine$double.eps
			U.B <- apply(sqrt(n)*X.B/sigma.B,2,max,na.rm=TRUE)
			L.B <- apply(sqrt(n)*X.B/sigma.B,2,min,na.rm=TRUE)

			if(is.null(alpha1)){
				c1.v <- sapply(seq(0,alpha,0.005),function(x){quantile(U.B, na.rm=TRUE, 1-x)})
				c2.v <- sapply(seq(alpha,0,-0.005),function(x){quantile(L.B, na.rm=TRUE, x)})
				pos.min <- which.min(c1.v-c2.v)
				alpha1 <- seq(0,alpha,0.005)[pos.min]
				alpha2 <- seq(alpha,0,-0.005)[pos.min]
				c1 <- as.numeric(c1.v[pos.min])
				c2 <- as.numeric(c2.v[pos.min])
        fixed.alpha1 <- FALSE
			}else{
				alpha1 <- alpha1
				alpha2 <- alpha-alpha1
				c1.v <- quantile(U.B, na.rm=TRUE, 1-alpha1)
				c2.v <- quantile(L.B, na.rm=TRUE, alpha2)
				c1 <- as.numeric(c1.v)
				c2 <- as.numeric(c2.v)
				fixed.alpha1 <- TRUE
			}

			L <- ROC.t - c1*sigma.B/sqrt(n)
			U <- ROC.t - c2*sigma.B/sqrt(n)
			L[which(L < 0)] <- 0
			U[which(U > 1)] <- 1
			L[which(L > 0.95)] <- 0.95	# Inferior band correction
			U[which(U < 0.05)] <- 0.05	# Superior band correction

			practical.area <- mean(U[-1]-L[-1]+U[-length(U)]-L[-length(L)])/2
			theoretical.area <- (c1-c2)/sqrt(n)*mean(sigma.B[-1]+sigma.B[-Ni-1],na.rm=TRUE)/2

			if(plot.var){
			  if(plot.bands){
			    par(mfrow=c(1,2))
			    }else{
			      par(mfrow=c(1,1))
			      }
			  plot(grid, sigma.B, type='l', xlab="t", ylab=expression(sigma(t)), main=expression(Standard~deviation~of~X(omega,t)~along~bootstrap~runs), cex.main=1.2-0.2*plot.var, font.main=2)
			  }else{
			    if(plot.bands){
			      par(mfrow=c(1,1))
			    }
			    }

			results <- list(method=method, B=B, conf.level=conf.level, ROC.t=ROC.t, s=s, alpha1=alpha1, alpha2=alpha2,  fixed.alpha1=fixed.alpha1, c1=c1, c2=c2, ROC.B=ROC.B, sd.PSN=sigma.B, practical.area=practical.area, theoretical.area=theoretical.area, L=L, U=U, Ni=Ni, ROC.t=ROC.t)


		# JMS confidence bands

		}else if(method=="JMS"){

			if(side!="right"){
				stop("Jensen et al. method may only be used to right-side ROC curves.")
			}

			if(!any(abs(grid - a.J) < .Machine$double.eps)){
				stop("a.J value must be stated in function of Ni value.")
			}else if(!any(abs(grid - b.J) < .Machine$double.eps)){
				stop("b.J value must be stated in function of Ni value.")
			}

			n <- length(cases); m <- length(controls)
			p <- grid[(which(abs(grid - a.J) < .Machine$double.eps) + 1) : (which(abs(grid - b.J) < .Machine$double.eps) - 1)]
			naiveROC <- ecdf(1-ecdf(controls)(cases))(p)
			Sp.step <- which(!duplicated(naiveROC))
			l.Sp.step <- length(Sp.step)
			smoothedROC1 <- (naiveROC[Sp.step[-l.Sp.step]] + naiveROC[Sp.step[-1]])/2
			smoothedROC2 <- naiveROC[Sp.step[-1]]
			smoothedROC <- c(naiveROC[Sp.step[1]]/2,as.vector(rbind(smoothedROC1,smoothedROC2)))
			Sp.smoothedROC <- c(Sp.step[1],as.vector(rbind(Sp.step[-l.Sp.step], (Sp.step[-l.Sp.step] + Sp.step[-1] )/2 )))
			smoothROC <- approxfun(c(0,p[Sp.smoothedROC],1), c(0,smoothedROC,1))

			f.controls <- approxfun(c(-1/.Machine$double.eps,density(controls)$x,1/.Machine$double.eps),c(0,density(controls)$y,0))
			f.cases <- approxfun(c(-1/.Machine$double.eps,density(cases)$x,1/.Machine$double.eps),c(0,density(cases)$y,0))

			F.controls <- ecdf(controls)
			F.cases <- ecdf(cases)

			lambda <- m/(n+m)
			C1 <- quantile(controls,1-p)
			C2 <- f.cases(C1); C3 <- f.controls(C1); C4 <- C2/C3
			C5 <- F.cases(C1)

			Var <- 1/lambda*C4^2*p*(1-p) + 1/(1-lambda)*C5*(1-C5)

			Stand.Psi <- sapply(1:B,function(b){
						B1.b <- as.vector(sde::BBridge(x=0, y=0, t0=0, T=1, N=Ni))
						B2.b <- as.vector(sde::BBridge(x=0, y=0, t0=0, T=1, N=Ni))
						B1 <- approxfun(grid,B1.b)
						B2 <- approxfun(grid,B2.b)
						Psi <- sqrt(1/lambda)*C4*B1(1-p) + sqrt(1/(1-lambda))*B2(C5)
						Stand.Psi <- abs(Psi)/sqrt(Var)
						})

			Extremum <- apply(Stand.Psi,2,max,na.rm=TRUE)
			K.alpha <- quantile(Extremum,1-alpha/2)

			smoothROC.p <- smoothROC(p)
			U <- smoothROC.p + K.alpha*sqrt(Var/(n+m))
			L <- smoothROC.p - K.alpha*sqrt(Var/(n+m))
			L[which(L < 0)] <- 0
			U[which(U > 1)] <- 1

			practical.area <- mean(U[-1]-L[-1]+U[-length(U)]-L[-length(L)])/2

			if(plot.var){
				if(plot.bands){
					par(mfrow=c(1,2))
				}else{
				  par(mfrow=c(1,1))
				}
					plot(p, Var, type='l', xlab="t", ylab=expression(sigma^{2}~(t)), main=expression(Variance~of~R(omega,t)), cex.main=1.2-0.2*(plot.var), font.main=2)
			}else{
        if(plot.bands){
          par(mfrow=c(1,1))
        }
			}

			results <- list(method=method, B=B, conf.level=conf.level, ROC.t=ROC.t, a.J=a.J, b.J=b.J, p=p, smoothROC.p=smoothROC.p, K.alpha=K.alpha, var.JMS=Var, practical.area=practical.area, L=L, U=U, Ni=Ni, ROC.t=ROC.t)


		# DEK confidence bands

		}else if(method=="DEK"){

			if(side!="right"){
				stop("Demidenko method may only be used to right-side ROC curves.")
			}

			superior.band.index <- function(xy) {
								i <- order(xy[, 1], xy[, 2], decreasing=FALSE)
								y <- xy[i, 2]
								frontier <- which(cummax(y) <= y)
								y.0 <- y[frontier]
								frontier <- frontier[c(TRUE, y.0[-1] != y.0[-length(y.0)])]
								return(i[frontier])
							}

			inferior.band.index <- function(xy) {
								i <- order(xy[, 1], xy[, 2], decreasing=TRUE)
								y <- xy[i, 2]
								frontier <- which(cummin(y) >= y)
								y.0 <- y[frontier]
 								frontier <- frontier[c(TRUE, y.0[-1] != y.0[-length(y.0)])]
								return(i[frontier])
							}

			n <- length(cases); m <- length(controls)
			m0 <- mean(controls); m1 <- mean(cases)
			s0 <- sd(controls); s1 <- sd(cases)
			lambda <- 1-alpha
			Gamma <- function(c,mu,s){(mu-c)/s}
			points <- seq(min((density(c(controls,cases)))$x), max((density(c(controls,cases)))$x), length.out=5000)
			l <- length(points)

			x <- matrix(rep(0,2*l),l,2); y <- matrix(rep(0,2*l),l,2)
			U <- matrix(rep(0,2*l),l,2); V <- matrix(rep(0,2*l),l,2)
			G0 <- matrix(rep(0,2*l),l,2); G1 <- matrix(rep(0,2*l),l,2)

			for(i in 2:(l-1)){
				c <- points[i]
				g0 <- Gamma(c,m0,s0); g1 <- Gamma(c,m1,s1)
				G0[i,] <- g0; G1[i,] <- g1
				A0 <- 1/(1/m + g0^2/(2*(m-1))); A1 <- 1/(1/n + g1^2/(2*(n-1)))
				B0 <- A0/s0; B1 <- A1/s1
				D0 <- g0*A0^2/(s0*(m-1)); D1 <- g1*A1^2/(s1*(n-1))
				q <- qchisq(lambda,2)

				p0 <- -B0^2*A0*q + D0^2*q^2
				p1 <- -2*A0*D0*B1*q
				p2 <- A0*A1*B0^2 + 2*D0*A0*D1*q + B1^2*A0^2 - 2*D0^2*A1*q
				p3 <- 2*D0*A0*A1*B1 - 2*A0^2*D1*B1
				p4 <- A1^2*D0^2 - 2*D0*D1*A0*A1 + A0^2*D1^2
				nu <- polyroot(c(p0,p1,p2,p3,p4))
				v <- suppressWarnings(as.numeric(nu[abs(Im(nu)) < sqrt(.Machine$double.eps)]))
				u <- ((A0*D1 - D0*A1)*v^2 - A0*B1*v - D0*q)/(A0*B0)
				V[i,] <- v; U[i,] <- u
				x[i,] <- sort(u+g0); y[i,] <- sort(v+g1)
			}

			U <- U[-c(1,l),]; V <- V[-c(1,l),]
			x <- x[-c(1,l),]; y <- y[-c(1,l),]
			G0 <- G0[-c(1,l),]; G1 <- G1[-c(1,l),]

			extremes1.matrix <- t(rbind(pnorm(as.vector(t(G0))),pnorm(as.vector(t(y)))))
			extremes2.matrix <- t(rbind(pnorm(as.vector(t(x))),pnorm(as.vector(t(G1)))))
			extremes.matrix <- rbind(extremes1.matrix,extremes2.matrix)
			sup.index <- superior.band.index(extremes.matrix)
			superior.band.points <- rbind(c(0,0), extremes.matrix[sup.index, , drop=FALSE], c(1,1))
			inf.index <- inferior.band.index(extremes.matrix)
			inferior.band.points <- rbind(c(1,1), extremes.matrix[inf.index, , drop=FALSE], c(0,0))
			inferior.band <- approxfun(inferior.band.points); L <- inferior.band(grid)
			superior.band <- approxfun(superior.band.points); U <- superior.band(grid)
			practical.area <- mean(U[-1]-L[-1]+U[-length(U)]-L[-length(U)])/2

      if(plot.bands){
        par(mfrow=c(1,1))
      }

			results <- list(method=method, conf.level=conf.level, DEK.fpr=pnorm(G0), DEK.tpr=pnorm(G1), practical.area=practical.area, L=L, U=U, Ni=Ni, ROC.t=ROC.t, x=x, y=y, inferior.band.points=inferior.band.points, superior.band.points=superior.band.points)

		}

  attr(results, 'class') <- 'rocbands'


  if(plot.bands){
    plot(results, cex.lab=1.2)
  }

return(results)

}

