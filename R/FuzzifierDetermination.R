#' @title FuzzifierDetermination
#' @name FuzzifierDetermination
#' @description
#' To determine the optimized fuzzifier value for FCM running.
#'
#' @usage FuzzifierDetermination(x, wv, max.m=10, stand=FALSE, dmetric='sqeuclidean')
#'
#' @param x Data.frame. the input Rrs data
#' @param wv Wavelength of X
#' @param max.m Set max.m as for determination of m.mub. Default as 10
#' @param stand Whether the data set is standardized (T needs not standardization)
#' @param dmetric Distance method. Default as 'sqeuclidean'
#'
#' @return \code{FD} list contains several result by \code{FuzzifierDetermination}:
#'   \itemize{
#'     \item \strong{x} The raw input Rrs dataframe with unit sr^-1
#'     \item \strong{x.stand} The standardized Rrs dataframe, if \code{stand=F}
#'     \item \strong{wv} Wavelength with unit nm
#'     \item \strong{max.m} The maximum fuzzifier of FCM as a restriction
#'     \item \strong{stand} A logic value for choosing whether to use standardization
#'     \item \strong{dmetric} A string value for choosing which distance metric to be used
#'     \item \strong{Area} A numeric vector for trapezoidal integral values
#'     \item \strong{m.ub} The upper boundary of fuzzifier(m) value
#'     \item \strong{m.used} The desired value of fuzzifier(m) value
#'   }
#'
#' @export
#'
#' @examples 
#' \dontrun{
#' library(FCMm)
#' library(tidyverse)
#' data("Nechad2015")
#' w <- Nechad2015 %>% names %>%
#'     str_extract(.,pattern="\\d") %>%
#'     is.na %>% {!.}
#' wv <- w %>% names(Nechad2015)[.] %>%
#'     gsub('X','',.) %>% as.numeric
#' x <- w %>% Nechad2015[,.]
#' names(x) <- wv
#' rm(w)
#' FD <- FuzzifierDetermination(x, wv, stand=F)
#' }
#'
#' @references 
#' \itemize{
#'   \item Bi S, Li Y, Xu J, et al. Optical classification of inland waters based on
#'     an improved Fuzzy C-Means method[J]. Optics Express, 2019, 27(24): 34838-34856.
#'   \item Dembele D. Multi-objective optimization for clustering 3-way gene
#'     expression data[J]. Advances in Data Analysis and Classification, 2008, 2(3):
#'     211-225.
#' }
#'

FuzzifierDetermination <- function(x, wv, max.m=10, stand=FALSE, dmetric="sqeuclidean"){

  # about x
  if(missing(x))
    stop("Missing input data set!")
  if(is.null(x))
    stop("Data set is null")
  if((is.data.frame(x))||(is.vector(x)))
    x <- as.matrix(x)
  if(!is.matrix(x))
    stop("Data set must be a vector, data frame or matrix")
  if(any(is.na(x)))
    stop("Data set should not contain NA values. Please remove NAs and try again")
  if(!is.numeric(x))
    stop("Data set must be a numeric vector, data frame or matrix")

  # about wv
  if(missing(wv))
    stop("Missing wavelength!")
  if(is.null(wv))
    stop("Wavelength is null")
  if(!is.vector(wv)){
    if(dim(wv)[2] != dim(x)[2]){
      stop("the dimension of wavelength and input data set is different! --A")}
    wv <- t(wv)
  }
  if(dim(t(wv))[2] != dim(x)[2]){
    stop("the dimension of wavelength and input data set is different! --B")}

  if(!is.numeric(wv))
    stop("Wavelength must be a numeric vector")

  # about max.m
  if(!is.matrix(x))
    stop("max.m must be a numeric")
  if(max.m > 50)
    stop("Are you kidding me? The max.m is too large!")

  # stand
  if(!is.logical(stand))
    stop("Stand must be TRUE or FALSE")

  # demteric
  if(dmetric != "sqeuclidean")
    stop("Sorry, this function could only support sqeuclidean distance metric!")

  Area <- NULL
  x.stand <- x
  if(stand==FALSE){
    Area <- .trapz(wv,x)
    x.stand <- x
    for(i in 1:ncol(x)){
      x.stand[,i] = x[,i] / Area
    }
  }

  m.ub <- .commub(x.stand,max.m)

  if(m.ub >= 10){
    warning("The upon boundary fuzzifier m is larger than 10!")
    stop("Maybe this Rrs data set is not good for clustering, Please check outliers and run again")
  }

  m.used <- 1+m.ub/10

  m <- m.used

  result = list()
  result$x <- x
  result$x.stand <- x.stand # x.stand is raw matrix if stand.flag is TRUE
  result$wv <- wv
  result$max.m <- max.m
  result$stand <- stand
  result$dmetric <- dmetric
  result$Area <- Area
  result$m.ub <- m.ub
  result$m.used <- m.used

  return(result)

}

#' @export
.commub <- function(x,max.m){
  dismt <- .compdismt(x)
  ii <- 1
  Ym <- NULL
  m.arr <- seq(from=1.1,to=max.m,by=0.1)
  for(m in m.arr){
    Ym[ii] <- abs(sd(dismt^(1/(m-1)))/mean(dismt^(1/(m-1))) - 0.03*dim(x)[2])/(0.03*dim(x)[2])
    ii <- ii + 1
  }
  if(which.min(Ym)==length(Ym))
    warning("The calculated m.ub is equal to max.m. Please set a larger max.m")
  m.ub <- m.arr[which.min(Ym)]
  return(m.ub)
}

#' @export
.compdismt <-function(x){
  x <- as.matrix(x)
  n <- dim(x)[1]
  arr <- array(dim=c(n^2-n,1))
  ind <- 1
  for(i in 1:n){
    for(j in 1:n){
      if(i != j){
        a <- as.vector(x[i,])
        b <- as.vector(x[j,])
        arr[ind] <- sum(t(a-b) * (a-b))
        # arr[ind] <- .compdist(a,b)
        ind <- ind+1
      }
    }
  }
  return(arr)
}
