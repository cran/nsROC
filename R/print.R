print.groc <- function(x, ...){

  obj <- x
  cat("Data was encoded with", obj$levels[1], "(controls) and", obj$levels[2], "(cases).\n")
  if (!is.null(obj$pvalue.wilcox)) {
    if (obj$side == "right") {
      cat("Wilcoxon rank sum test:\n \t alternative hypothesis: median(controls) < median(cases); p-value =",  as.numeric(format(mean(obj$pvalue.wilcox), digits = 4)), "\n")
    }
    else {
      cat("Wilcoxon rank sum test:\n \t alternative hypothesis: median(cases) < median(controls); p-value =", as.numeric(format(mean(obj$pvalue.wilcox), digits = 4)), "\n")
    }
  }
  printside <- function(side) {
    switch(side, right = cat("It is assumed that larger values of the marker indicate larger confidence that a given subject is a case.\n"),
           left = cat("It is assumed that lower values of the marker indicate larger confidence that a given subject is a case.\n"),
           both = cat("It is assumed that both lower and larges values of the marker indicate larger confidence that a given subject is a case.\n"),
           both2 = cat("It is assumed that both lower and larges values of the marker indicate larger confidence that a given subject is a controls.\n"))
  }
  printside(obj$side)
  if (!is.null(obj$Ni)) {
    cat("Detailed ROC curve estimation was computed (time consuming) with", obj$Ni, "equidistant FPRs in (0,1).")
  }
  if (!is.null(obj$Ni) && obj$side == "both" && length(obj$index) != nrow(obj$coordinates)) {
    cat("There are pairs of points returning the same sensitivity and specificity (See pairpoints.coordinates)", "\n")
  }
  cat("There are", length(obj$controls), "controls and", length(obj$cases), "cases.\n")
  cat("The area under the ROC curve (AUC) is ", round(obj$auc, 3), ".\n", sep = "")
}




print.rocbands <- function(x, ...){

  obj <- x

  if(obj$method=="PSN"){

    cat("The method considered to build confidence bands is the one proposed in Martinez-Camblor et al. (2016).\n")
    cat("Confidence level (1-alpha): ", obj$conf.level, ".\n", sep="")
    cat("Bootstrap replications: ", obj$B, ".\n", sep="")
    cat("Scale parameter (bandwidth construction): ", obj$s, ".\n", sep="")
    if(obj$fixed.alpha1==FALSE){
      cat("The optimal confidence band is reached for alpha1 = ", obj$alpha1," and alpha2 = ", obj$alpha2,".\n", sep="")
    }else{
      cat("alpha1: ", obj$alpha1,".\n", sep="")
    }
    cat("The area between the confidence bands is ", round(obj$practical.area,4) ," (theoretically ", round(obj$theoretical.area,4), ").\n", sep="")

  }else{

    if(obj$method=="JMS"){

      cat("The method considered to build confidence bands is the one proposed in Jensen et al. (2000).\n")
      cat("Confidence level (1-alpha): ", obj$conf.level, ".\n", sep="")
      cat("Bootstrap replications: ", obj$B, ".\n", sep="")
      cat("Interval in which compute the regional confidence bands: (", obj$a.J, ",", obj$b.J, ").\n", sep="")
      cat("K.alpha: ", obj$K.alpha, ".\n", sep="")
      cat("The area between the confidence bands is ", round(obj$practical.area,4),".\n", sep="")

    }else{

      cat("The method considered to build confidence bands is the one proposed in Demidenko (2012).\n")
      cat("Confidence level (1-alpha): ", obj$conf.level, ".\n", sep="")
      cat("The area between the confidence bands is ", round(obj$practical.area,4),".\n", sep="")

    }
  }

}


print.cdroc <- function(x, ...){

  obj <- x

  cat('cdroc object: \n')
  cat('Number of cut points:', length(obj$cutPoints), '\n')
  cat('Method:', obj$method, '\n')
  cat('predict.time:', obj$predict.time, '\n')
  cat('AUC:', obj$auc, '\n')
  if(obj$method=='wKM'){
    cat('Kernel function:', obj$kernel, '\n')
    cat('Kernel bandwidth:', obj$h, '\n')
  }
  if (!is.null(obj$ci)){
    cat('Bootstrap number of replicates:', obj$boot.n, '\n')
    cat('Bootstrap AUC:', round(obj$meanAuc,3), '\n')
    cat('Bootstrap AUC Confidence Interval:', round(obj$ciAuc,3), '\n')
    cat('Bootstrap AUC Confidence Level:', obj$conf.level, '\n')
    cat('Bootstrap seed:', obj$seed, '\n')
  }

}
