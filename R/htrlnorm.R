#' @name htrlnorm
#' @rdname htrlnorm
#' 
#' @aliases mle.htrlnorm,dhtrlnorm,phtrlnorm,qhtrlnorm
#'
#' @title Hurdle Truncated Log-Normal Model
#' 
#' @description mle parameter estimation, random generation, density and 
#' quantile function for the hurdle truncated log-normal distribution.
#'
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations.
#' @param meanlog mean in the log scale.
#' @param sd standard deviation in log scale.
#' @param warning.silent logical; if TRUE suppress all warning messages from 
#' \code{\link{fitdistr}} during the mle process.
NULL



#' @rdname htrlnorm
#' @importFrom MASS fitdistr
#' @export
mle.htrlnorm <- function(x, warning.silent=TRUE){
  
  if(any(x<0)) stop("find negative values in x")
  if(!is.vector(x)) stop("x must be a vector")
  
  x <- log(x+1)
  x1 <- x[x>0]
  
  n <- length(x)
  n1 <- length(x1)
  n0 <- n - n1
  
  phi <- n0/n
  
  param.trnorm.mle <- NULL
  try(param.trnorm.mle <- MASS::fitdistr(x=x1, densfun=truncnorm::dtruncnorm,
                                         start=list(mean=mean(x1),sd=sd(x1)),
                                         a=0,b=Inf),
      silent=warning.silent
  )
  
  if(!is.null(param.trnorm.mle) && phi!=0){
    loglik <- n0*log(phi) + n1*log(1-phi) + param.trnorm.mle$loglik
  }else if (!is.null(param.trnorm.mle) && phi==0){
    loglik <- param.trnorm.mle$loglik
  }else{
    loglik <- NA
  }
  
  
  if(!is.null(param.trnorm.mle)){
    return(c("phi"    =phi,
             "meanlog"=as.numeric(param.trnorm.mle$estimate[1]),
             "sdlog"  =as.numeric(param.trnorm.mle$estimate[2]),
             "loglik" =loglik,
             "success"= 1)
    )
  } else {
    return(c("phi"    =phi,
             "meanlog"=NA,
             "sdlog"  =NA,
             "loglik" =NA,
             "success"=0)
    )
  }
  
}



#'@rdname htrlnorm
#'@importFrom truncnorm dtruncnorm
#'@export
dhtrlnorm <- function(xi, phi=0, meanlog=0, sdlog=1){
  
  #CHECK ARGUMENTS
  if(!is.numeric(phi) | !is.numeric(meanlog) | !is.numeric(sdlog)) stop("phi, meanlog, sdlog must be all number")
  if(phi<0 || phi>1) stop("phi must be in range [0,1]")
  if(sd<=0) stop("sd must be greater than 0")
  
  ans <- (1-phi)*dtruncnorm(xi,a=0,b=Inf,mean=meanlog,sd=sdlog)
  ans[y==0] <- phi
  
  return(ans)
}

#'@rdname htrlnorm
#'@importFrom truncnorm qtruncnorm
#'@export
qhtrlnorm <- function(p, phi=0, meanlog=0, sdlog=1){
  
  #CHECK ARGUMENTS
  if(!is.numeric(phi) | !is.numeric(meanlog) | !is.numeric(sdlog)) stop("phi, meanlog, sdlog must be all number")
  if(phi<0 || phi>1) stop("phi must be in range [0,1]")
  if(sdlog<0) stop("sdlog must be greater than 0")
  
  ans <- rep(NA_real_,length(p))
  ans[p<=phi] <- 0
  
  idx <- is.na(ans)
  ans[idx]  <- qtruncnorm(p=(p[idx]-phi)/(1-phi), a=0, b=Inf,
                          mean=meanlog, sd=sdlog)
  
  ans <- exp(ans) - 1
  ans[intersect(which(ans>0),which(ans<1))] <- 1
  ans <- round(ans)
  
  return(ans)
}

#' @rdname htrlnorm
#' @importFrom truncnorm rtruncnorm
#' @export
rhtrlnorm <- function(N, phi=0, meanlog=0, sdlog=1){
  
  if(!is.numeric(N) | !is.numeric(phi) | !is.numeric(meanlog) | !is.numeric(sdlog)) stop("N, phi, meanlog, sdlog must be all number")
  if(round(N)!=N) stop("N must be an integer")
  if(N<1) stop("N must be a positive integer number")
  if(phi<0 || phi>1) stop("phi must be a number in [0,1]")
  if(sdlog<0) stop("sdlog must be a positive number")
  
  ans <- stats::rbinom(n=N, size=1, prob=1-phi)
  m <- length(ans[ans>0])
  ans[ans==1] <- truncnorm::rtruncnorm(n=m, a=0, b=Inf,
                                       mean=meanlog,sd=sdlog)
  
  ans <- exp(ans) - 1
  ans[intersect(which(ans>0),which(ans<1))] <- 1
  ans <- round(ans)
  
  return(ans)
}
