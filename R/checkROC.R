
checkROC <- function(X, D){

  D <- as.factor(D)
  levels <- base::levels(D)

  if( missing(X) | missing(D) | is.null(X) | is.null(D) | !is.numeric(X) | sum(!is.na(X))==0 | sum(!is.na(D))==0 ){
    if( missing(X) | is.null(X) | sum(!is.na(X))==0 ){stop("Marker values should be provided.")}
    else if( missing(D) | is.null(D) | sum(!is.na(D))==0 ){stop("Response values should be provided.")}
    else if( !is.numeric(X) ){stop("Marker values are not numeric.")}
  }
  else{
    if( sum(is.na(X))!=0 ){
      D <- D[!is.na(X)]
      X <- X[!is.na(X)]
      warning("Some of marker values are NA. They have been removed.")
    }
    if( sum(is.na(D))!=0 ){
      X <- X[!is.na(D)]
      D <- D[!is.na(D)]
      warning("Some of response values are NA. They have been removed.")
    }
    if( length(X)!=length(D)){
      stop("Marker and response vectors should have the same length.")
    }

    if( length(levels)<2 ){
      stop("Response vector must contain at least two different values.")
    }
    if( length(levels)>2 ){
      warning("There are more than two levels in response vector. Only two are considered.")
    }
  }

  # Get controls and cases data and specify which level corresponds to each one

  controls <- split(X,D)[[levels[1]]]
  cases <- split(X,D)[[levels[2]]]

  # Calculate sizes and final data considered

  n0 <- length(controls); n1 <- length(cases)
  X <- c(controls,cases)
  D <- c(rep(levels[1],n0),rep(levels[2],n1))

  list(levels=levels, controls=controls, cases=cases, n0=n0, n1=n1, X=X, D=D)

}
