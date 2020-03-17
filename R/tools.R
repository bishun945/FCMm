#' @title Assess the performance of algorithms
#' @name cal.metrics
#' @usage cal.metrics(x,y,name="all",log10=FALSE)
#' @param x True value OR Actual value
#' @param y Estimated value OR Predict value
#' @param name The name of metrics
#' @param log10 Logical. Whether the input x and y should be log10-transformed
#' @export
#' @importFrom stats cor cor.test median na.omit sd
#' @importFrom stats confint cov lm qt var
#' @family Utils
#' @note (2020-02-09) All functions used log10 transformation was assgined to the 
#'   Key parameter `log10`
#' 
cal.metrics <- function(x,y,name="all",log10=FALSE){

  name <- match.arg(name, cal.metrics.names())
  
  if(log10){
    x <- log10(x)
    y <- log10(y)
  }
  
  # RMSE Series (Root mean squared error)
  .cal.rmse <- function(x,y){
    return(sqrt(mean((x-y)^2, na.rm=T)))
  }
  
  .cal.crmse<-function(x,y){
    return(sqrt(mean((((x-mean(x, na.rm=T))-(y-mean(y, na.rm=T))))^2, na.rm=T)))
  }
  
  # MAE Series (Mean absolute error)
  .cal.mae<-function(x,y){
    return(mean(abs(y-x), na.rm=T))
  }
  
  .cal.mdae<-function(x,y){
    return(median(abs(y-x), na.rm=T))
  }

  .cal.nmae_sd <- function(x,y){
    return(mean(abs(y-x), na.rm=T)/sd(y, na.rm=T))
  }
  
  .cal.nmdae_sd <- function(x,y){
    return(median(abs(y-x), na.rm=T)/sd(y, na.rm=T))
  }
  
  # MAPE Series (Mean absolute percent error)
  .cal.mape<-function(x,y){
    return(mean(abs((y-x)/x), na.rm=T)*100)
  }
  
  .cal.smape <- function(x,y){
    return(mean(       abs(2*(y-x))/(abs(x)+abs(mean(x, na.rm=T)))       , na.rm=T)*100)
  }
  
  .cal.smdape <- function(x,y){
    return(median(       abs(2*(y-x))/(abs(x)+abs(mean(x, na.rm=T)))       , na.rm=T)*100)
  }
  
  .cal.smrpe <- function(x,y){
    return(mean(       (2*(y-x))/((x)+(mean(x, na.rm=T)))       , na.rm=T)*100)
  }
  
  .cal.smdrpe <- function(x,y){
    return(median(       (2*(y-x))/((x)+(mean(x, na.rm=T)))       , na.rm=T)*100)
  }
  
  .cal.mdape<-function(x,y){
    return(median(abs((y-x)/x), na.rm=T)*100)
  }
  
  .cal.mrpe<-function(x,y){
    return(mean((y-x)/x, na.rm=T)*100)
  }
  
  .cal.mdrpe<-function(x,y){
    return(median((y-x)/x, na.rm=T)*100)
  }
  
  .cal.umrpe <- function(x,y){
    return(mean((y-x)/(0.5*x+0.5*y)*100, na.rm=T))
  }
  
  .cal.umdrpe <- function(x,y){
    return(median((y-x)/(0.5*x+0.5*y)*100, na.rm=T))
  }
  
  # Bias Series 
  .cal.bias<-function(x,y){
    return(mean((y-x), na.rm=T))
  }
  
  .cal.md_bias<-function(x,y){
    return(median((y-x), na.rm=T))
  }
  
  .cal.ratio <- function(x,y){
    return(mean(y/x, na.rm=T))
  }
  
  .cal.md_ratio <- function(x,y){
    return(median(y/x, na.rm=T))
  }
  
  .cal.cv<-function(x,y){
    return(sd(y, na.rm=T) / mean(y, na.rm=T) * 100)
  }
  
  .cal.R<-function(x,y){
    result <- cor.test(x,y,method = "pearson")
    result <- as.numeric(result$estimate)
    return(result)
  }
  
  .cal.R2<-function(x,y){
    result <- cor.test(x,y,method = "pearson")
    result <- as.numeric(result$estimate) ^ 2
    return(result)
  }
  
  # Linear regression relationship Series
  .cal.slope<-function(x,y){
    tmp <- na.omit(data.frame(x=x,y=y))
    if(dim(tmp)[1] <= 2){
     result <- NA 
    }else{
      # ft <- lm(x~y,data=tmp)
      ft <- lmodel2_(x~y,data=tmp)
      if(ft$r < 0.5){
        result <- NA
      }else{
        result <- as.numeric(ft$regression.results[which(ft$regression.results$Method == "SMA"),"Slope"])
      }
    }
    return(result)
  }
  
  .cal.intercept<-function(x,y){
    tmp <- na.omit(data.frame(x=x,y=y))
    if(dim(tmp)[1] <= 2){
      result <- NA 
    }else{
      # ft <- lm(x~y,data=tmp)
      ft <- lmodel2_(x~y,data=tmp)
      if(ft$r < 0.5){
        result <- NA
      }else{
        result <- as.numeric(ft$regression.results[which(ft$regression.results$Method == "SMA"),"Intercept"])
      }
    }
    return(result)
  }
  
  .cal.R2_SMA<-function(x,y){
    tmp <- na.omit(data.frame(x=x,y=y))
    if(dim(tmp)[1] <= 2){
      result <- NA 
    }else{
      # ft <- lm(x~y,data=tmp)
      ft <- lmodel2_(x~y,data=tmp)
      result <- ft$rsquare
    }
    return(result)
  }
  
  .cal.eta<-function(x,y){
    return(length(which(is.na(y) == FALSE)) / length(y))
  }
  
  .cal.chi2<-function(x,y){
    return(sum((x-y)^2/y, na.rm=T))
  }

  if(name == "RMSE"){
    result <- .cal.rmse(x,y)
  }else if(name == "CRMSE"){
    result <- .cal.crmse(x,y)
  }else if(name == "MAE"){
    result <- .cal.mae(x,y)
  }else if(name == "MDAE"){
    result <- .cal.mdae(x,y)
  }else if(name == "NMAE_SD"){
    result <- .cal.nmae_sd(x,y)
  }else if(name == "NMDAE_SD"){
    result <- .cal.nmdae_sd(x,y)
  }else if(name == "MAPE"){
    result <- .cal.mape(x,y)
  }else if(name == "SMAPE"){
    result <- .cal.smape(x,y)
  }else if(name == "SMDAPE"){
    result <- .cal.smdape(x,y)
  }else if(name == "SMRPE"){
    result <- .cal.smrpe(x,y)
  }else if(name == "SMDRPE"){
    result <- .cal.smdrpe(x,y)
  }else if(name == "MDAPE"){
    result <- .cal.mdape(x,y)
  }else if(name == "MRPE"){
    result <- .cal.mrpe(x,y)
  }else if(name == "MDRPE"){
    result <- .cal.mdrpe(x,y)
  }else if(name == "UMRPE"){
    result <- .cal.umrpe(x,y)
  }else if(name == "UMDRPE"){
    result <- .cal.umdrpe(x,y)
  }else if(name == "BIAS"){
    result <- .cal.bias(x,y)
  }else if(name == "MD_BIAS"){
    result <- .cal.md_bias(x,y)
  }else if(name == "RATIO"){
    result <- .cal.ratio(x,y)
  }else if(name == "MD_RATIO"){
    result <- .cal.md_ratio(x,y)
  }else if(name == "CV"){
    result <- .cal.cv(x,y)
  }else if(name == "R"){
    result <- .cal.R(x,y)
  }else if(name == "R2"){
    result <- .cal.R2(x,y)
  }else if(name == "SLOPE"){
    result <- .cal.slope(x,y)
  }else if(name == "INTERCEPT"){
    result <- .cal.intercept(x,y)
  }else if(name == "R2_SMA"){
    result <- .cal.R2_SMA(x,y)
  }else if(name == "eta"){
    result <- .cal.eta(x,y)
  }else if(name == "chi2"){
    result <- .cal.chi2(x,y)
  }
  else if(name == "all"){
    result <- list(RMSE      = .cal.rmse(x,y),
                   CRMSE     = .cal.crmse(x,y),
                   MAE       = .cal.mae(x,y),
                   MDAE      = .cal.mdae(x,y),
                   NMAE_SD   = .cal.nmae_sd(x,y),
                   NMDAE_SD  = .cal.nmdae_sd(x,y),
                   MAPE      = .cal.mape(x,y),
                   SMAPE     = .cal.smape(x,y),
                   SMDAPE    = .cal.smdape(x,y),
                   SMRPE     = .cal.smrpe(x,y),
                   SMDRPE    = .cal.smdrpe(x,y),
                   MDAPE     = .cal.mdape(x,y),
                   MRPE      = .cal.mrpe(x,y),
                   MDRPE     = .cal.mdrpe(x,y),
                   UMRPE     = .cal.umrpe(x,y),
                   UMDRPE    = .cal.umdrpe(x,y),
                   BIAS      = .cal.bias(x,y),
                   MD_BIAS   = .cal.md_bias(x,y),
                   RATIO     = .cal.ratio(x,y),
                   MD_RATIO  = .cal.md_ratio(x,y),
                   CV        = .cal.cv(x,y),
                   R         = .cal.R(x,y),
                   R2        = .cal.R2(x,y),
                   SLOPE     = .cal.slope(x,y),
                   INTERCEPT = .cal.intercept(x,y),
                   R2_SMA    = .cal.R2_SMA(x,y),
                   eta       = .cal.eta(x,y),
                   chi2      = .cal.chi2(x,y)
    )
  }
  return(result)
}

#' @title List all metrics in function cal.metrics
#' @name cal.metrics.names
#' @export
#' @return strings
#' @family Utils
cal.metrics.names <- function(){
  c('RMSE','CRMSE',
    'MAE','MDAE','NMAE_SD','NMDAE_SD',
    'MAPE', 'SMAPE', 'SMRPE', 'MDAPE',"MRPE",'MDRPE',"UMRPE",'UMDRPE',
    'SMDAPE','SMDRPE',
    'BIAS', 'MD_BIAS', 'RATIO', 'MD_RATIO',
    'CV','R','R2',
    'SLOPE', 'INTERCEPT', 'R2_SMA',
    'eta','chi2','all')
}

#' @title Assess the performance of algorithms as vector
#' @name cal.metrics.vector
#' @usage cal.metrics.vector(x,y,name="all",log10=FALSE)
#' @param x True value OR Actual value
#' @param y Estimated value OR Predict value
#' @param name The name of metrics
#' @param log10 Logical. Whether the input x and y should be log10-transformed
#' @export
#' @importFrom stats cor cor.test median na.omit sd
#' @family Utils
#' 
cal.metrics.vector <- function(x,y,name="all",log10=FALSE){

  name <- match.arg(name, cal.metrics.vector.names())
  
  if(log10){
    x <- log10(x)
    y <- log10(y)
  }
  
  .cal.mae<-function(x,y){
    return(abs(y-x))
  }
  
  .cal.nmae<-function(x,y){
    return(abs(y-x) / sd(y, na.rm=T))
  }
  
  .cal.mape<-function(x,y){
    return(abs((y-x)/x)*100)
  }
  
  .cal.smape <- function(x,y){
    return(abs(2*(y-x))/(abs(x)+abs(mean(x, na.rm=T)))*100)
  }
  
  .cal.smrpe <- function(x,y){
    return((2*(y-x))/((x)+(mean(x, na.rm=T)))*100)
  }
  
  .cal.nmape<-function(x,y){
    return(abs((y-x)/x) / sd(y, na.rm=T) *100)
  }
  
  .cal.mrpe<-function(x,y){
    return((y-x)/x*100)
  }
  
  .cal.nmrpe<-function(x,y){
    return((y-x)/x / sd(y, na.rm=T) *100)
  }
  
  .cal.umrpe <- function(x,y){
    return((y-x)/(0.5*x+0.5*y)*100)
  }
  
  .cal.bias<-function(x,y){
    return(y-x)
  }
  
  .cal.ratio<-function(x,y){
    return(y/x)
  }

  if(name == "MAE"){
    result <- .cal.mae(x,y)
  }else if(name == "NMAE"){
    result <- .cal.nmae(x,y)
  }else if(name == "MAPE"){
    result <- .cal.mape(x,y)
  }else if(name == "SMAPE"){
    result <- .cal.smape(x,y)
  }else if(name == "SMRPE"){
    result <- .cal.smrpe(x,y)
  }else if(name == "NMAPE"){
    result <- .cal.nmape(x,y)
  }else if(name == "MRPE"){
    result <- .cal.mrpe(x,y)
  }else if(name == "NMRPE"){
    result <- .cal.nmrpe(x,y)
  }else if(name == "UMRPE"){
    result <- .cal.umrpe(x,y)
  }else if(name == "BIAS"){
    result <- .cal.bias(x,y)
  }else if(name == "RATIO"){
    result <- .cal.ratio(x,y)
  }
  else if(name == "all"){
    result <- list(MAE   = .cal.mae(x,y),
                   NMAE  = .cal.nmae(x,y),
                   MAPE  = .cal.mape(x,y),
                   SMAPE = .cal.smape(x,y),
                   SMRPE = .cal.smrpe(x,y),
                   NMAPE = .cal.nmape(x,y),
                   MRPE  = .cal.mrpe(x,y),
                   NMRPE = .cal.nmrpe(x,y),
                   UMRPE = .cal.umrpe(x,y),
                   BIAS  = .cal.bias(x,y),
                   RATIO = .cal.ratio(x,y)
    )
  }
  return(result)
}

#' @title List all metrics in function cal.metrics
#' @name cal.metrics.vector.names
#' @export
#' @return strings
#' @family Utils
cal.metrics.vector.names <- function(){
  c('MAE', 'NMAE',
    'MAPE', 'SMAPE', 'NMAPE',
    'MRPE', 'SMRPE', 'NMRPE', 'UMRPE',
    'BIAS', 'RATIO', 'all')
}

#' @title Standard deviation by trimmed values
#' @name trim_sd
#' @param x input numeric vector
#' @param trim percent of length(x) to be trimmed
#' @param na.rm default as TRUE to remove NA values
#' @export
#' @return sd value
#' @family Utils
trim_sd <- function(x, trim=0.05, na.rm=T){
  
  stopifnot(is.vector(x) & is.numeric(x) & is.numeric(trim))
  
  x <- as.numeric(na.omit(x))
  
  n_trim <- ceiling(length(x) * trim)
  
  head = n_trim + 1
  tail = length(x) - n_trim
  
  x_sort <- sort(x)
  
  x_new <- x_sort[head:tail]
  
  if(length(x_new) <= 1){
    return(sd(x))
  }else{
    return(sd(x_new))
  }

}


#' @title Convert dataframe with factor to character
#' @name .level_to_variable
#' @param dt dataframe
#' @param warn warnning option
#' @export
#' @family Utils
.level_to_variable <- function(dt, warn=F){
  if(!is.data.frame(dt))
    stop('Input must be a data.frame')
  w <- NULL
  for(i in names(dt)){w<-c(w,is.factor(dt[,i]))}
  if(sum(w) != 0){
    for(i in names(dt)[w]){
      tmp <- levels(dt[,i])[dt[,i]]
      dt[,i] <- tmp
      }
    if(warn) message('Factors turned to variables are: ')
    if(warn) message(paste(names(dt)[w],collapse=' '))
    return(dt)
  }else{
    if(warn) message('No factor in this data.frame.')
    return(dt)
  }
}

.compdist <- function(a, b, dmetric="sqeuclidean", pw=2){
  if(missing(a))
    stop("Missing data arguments to compute the distances")
  if(missing(b))
    b <- a
  if(!is.numeric(a) || !is.numeric(b))
    stop("Input data arguments must be numeric to compute the distances")
  dmetrics <- c("euclidean", "sqeuclidean", "manhattan", "minkowski", "chebyshev",
                "pearsonchi", "neymanchi", "sqchi", "divergence", "addsymchi", "prosymchi", "clark",
                "canberra", "sorensen", "lorentzian", "sqchord", "cosine", "correlation")
  dmetric <- match.arg(dmetric, dmetrics)
  if(dmetric=="euclidean")
    distance <- sqrt(sum(t(a-b) * (a-b)))
  else if (dmetric=="sqeuclidean")
    distance <- sum(t(a-b) * (a-b))
  else if (dmetric=="manhattan")
    distance <- sum(abs(a-b))
  else if (dmetric=="minkowski")
    distance <- sum(abs(a-b)^pw)^(1/pw)
  else if (dmetric=="chebyshev")
    distance <- max(abs(a-b))
  else if (dmetric=="pearsonchi")
    distance <- sum((t(a-b) * (a-b) / b))
  else if (dmetric=="neymanchi")
    distance <- sum(t(a-b) * (a-b) / a)
  else if (dmetric=="sqchi")
    distance <- sum(t(a-b) * (a-b) / (a+b))
  else if (dmetric=="addsymchi")
    distance <- sum((t(a-b) * (a-b)) * ((a+b) / (a*b)))
  else if (dmetric=="prosymchi")
    distance <- 2 * sum((t(a-b) * (a-b)) / (a+b))
  else if (dmetric=="divergence")
    distance <- 2 * sum((t(a-b) * (a-b)) / (t(a+b) * (a+b)))
  else if (dmetric=="clark")
    distance <- sqrt(sum(abs((a-b) / (a+b))^2))
  else if (dmetric=="sqchord")
    distance <- sum(sqrt(a)-sqrt(b))^2
  else if (dmetric=="canberra")
    distance <- sum(abs(a-b)/(a+b))
  else if (dmetric=="sorensen")
    distance <- sum(abs(a-b))/sum(a+b)
  else if (dmetric=="lorentzian")
    distance <- sum(log(1+abs(a-b), exp(1)))
  else if (dmetric=="cosine")
    distance <- sum(a%*%b)/(sqrt(sum(a*a)) * sqrt(sum(b*b)))
  else if (dmetric=="correlation")
    distance <- (1-cor(a,b))/2
  return(distance)
}

.trapz <- function(wv,x){
  Area <-  as.matrix(x[,1])
  for(i in 1:dim(x)[1]){
    y <- x[i,]
    len <- length(y)
    Area[i] <- 0.5 * sum(diff(wv,1)*(y[1:(len-1)]+y[2:len]))
  }
  return(as.matrix(Area))
}


#' @title Run SMA linear regression
#' @name lmodel2_
#' @export
#' @return A lmodel2 function list
#' @param formula See package \code{lmodel2} for more details
#' @param data See package \code{lmodel2} for more details
#' @param range.y See package \code{lmodel2} for more details
#' @param range.x See package \code{lmodel2} for more details
#' @param nperm See package \code{lmodel2} for more details
#' @note This function was modified by the package \code{lmodel2} (Version 1.7-3) since 
#'   the message as output is unnecessary in most of scene. I just remove this messages 
#'   in by comment.
#' @import lmodel2
#' @family Utils
lmodel2_ <- function(formula, data = NULL, range.y = NULL, range.x = NULL, nperm = 0)
{
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action", 
               "offset"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  y <- mf[, 1]
  x <- mf[, 2]
  n <- nrow(mf)
  RMA <- FALSE
  if ((length(range.y) != 0) && (length(range.x) != 0)) {
    RMA <- TRUE
  }
  else {
    if (length(range.y) != 0) 
      stop("Give a value (relative or interval) to parameter 'range.x' for ranging")
    if (length(range.x) != 0) 
      stop("Give a value (relative or interval) to parameter 'range.y' for ranging")
  }
  if (!RMA) 
    # message("RMA was not requested: it will not be computed.", 
    #         "\n")
  if (nperm < 0) 
    stop("'nperm' cannot be negative")
  if (length(x) != n) 
    stop("Vectors y and x are not of the same length")
  if (var(y) == 0)
    stop("Variance(y) = 0. The regression coefficients cannot be computed")
  if (var(x) == 0)
    stop("Variance(x) = 0. The regression coefficients cannot be computed")
  ybar <- mean(y)
  xbar <- mean(x)
  yx <- cbind(y, x)
  cov.mat <- cov(yx)
  vary <- cov.mat[1, 1]
  varx <- cov.mat[2, 2]
  r <- cor(y, x)
  rsquare <- r^2
  info.slope <- 0
  info.CI <- 0
  t <- qt(0.025, (n - 2), lower.tail = FALSE)
  epsilon <- .Machine$double.eps
  met <- "OLS"
  temp <- lm(y ~ x)
  b.ols <- summary(temp)$coefficients[2, 1]
  angle <- atan(b.ols) * 180/pi
  P.param <- summary(temp)$coefficients[2, 4]
  temp2 <- confint(temp)
  res1 <- summary(temp)$coefficients[1, 1]
  res2 <- b.ols
  res3 <- angle
  res4 <- temp2[1, 1]
  res5 <- temp2[1, 2]
  res6 <- temp2[2, 1]
  res7 <- temp2[2, 2]
  res8 <- NA
  met <- c(met, "MA")
  MA.res <- MA.reg(b.ols, r, rsquare, ybar, xbar, epsilon, 
                   1)
  b.ma <- MA.res[2]
  if (rsquare <= epsilon) {
    info.slope <- 1
    warning("MA:  R-square = 0", "\n")
  }
  CL <- CLma(b.ma, n, r, rsquare, xbar, ybar, cov.mat, t, 1, 
             epsilon)
  lambda1 <- CL[5]
  lambda2 <- CL[6]
  H <- CL[7]
  A <- CL[8]
  if (CL[9] == 1) 
    info.CI <- 1
  res1 <- c(res1, MA.res[1])
  res2 <- c(res2, MA.res[2])
  res3 <- c(res3, MA.res[3])
  res4 <- c(res4, CL[1])
  res5 <- c(res5, CL[2])
  res6 <- c(res6, CL[3])
  res7 <- c(res7, CL[4])
  res8 <- c(res8, NA)
  met <- c(met, "SMA")
  SMA.res <- SMA.reg(b.ols, r, rsquare, ybar, xbar, epsilon)
  b.sma <- SMA.res[2]
  if (rsquare <= epsilon) {
    info.slope <- 1
    warning("SMA: R-square = 0", "\n")
  }
  CL <- CLsma(b.sma, n, r, rsquare, xbar, ybar, t, epsilon)
  res1 <- c(res1, SMA.res[1])
  res2 <- c(res2, SMA.res[2])
  res3 <- c(res3, SMA.res[3])
  res4 <- c(res4, CL[1])
  res5 <- c(res5, CL[2])
  res6 <- c(res6, CL[3])
  res7 <- c(res7, CL[4])
  res8 <- c(res8, NA)
  yx.2 = yx
  if (RMA) {
    met = c(met, "RMA")
    range.yx = apply(yx, 2, range)
    if (range.y == "relative") {
      if (range.yx[1, 1] < 0) 
        stop("Negative values in 'y'. Use 'interval' for ranging")
      range.yx[1, 1] <- 0
    }
    if (range.x == "relative") {
      if (range.yx[1, 2] < 0) 
        stop("Negative values in 'x'. Use 'interval' for ranging")
      range.yx[1, 2] <- 0
    }
    yx.1 <- sweep(yx, 2, range.yx[1, ], "-")
    yx.2 <- sweep(yx.1, 2, (range.yx[2, ] - range.yx[1, ]), 
                  "/")
    ran.y <- (range.yx[2, 1] - range.yx[1, 1])
    ran.x <- (range.yx[2, 2] - range.yx[1, 2])
    ratio <- ran.y/ran.x
    temp.ranged <- lm(yx.2[, 1] ~ yx.2[, 2])
    b.ols.ranged <- summary(temp.ranged)$coefficients[2, 
                                                      1]
    RMA.res <- MA.reg(b.ols.ranged, r, rsquare, ybar, xbar, 
                      epsilon, ratio)
    b.rma <- RMA.res[2]
    if (rsquare <= epsilon) {
      info.slope <- 1
      warning("RMA: R-square = 0")
    }
    cov.rma <- cov(yx.2)
    CL <- CLma(b.rma, n, r, rsquare, xbar, ybar, cov.rma, 
               t, ratio, epsilon)
    if (CL[9] == 1) 
      info.CI <- 1
    res1 <- c(res1, RMA.res[1])
    res2 <- c(res2, RMA.res[2])
    res3 <- c(res3, RMA.res[3])
    res4 <- c(res4, CL[1])
    res5 <- c(res5, CL[2])
    res6 <- c(res6, CL[3])
    res7 <- c(res7, CL[4])
    res8 <- c(res8, NA)
  }
  sdy <- sqrt(vary)
  sdx <- sqrt(varx)
  theta <- 90 - sign(r) * (atan(r * sdx/sdy) * 180/pi + atan(r * 
                                                               sdy/sdx) * 180/pi)
  if (nperm == 0) 
    # message("No permutation test will be performed")
  if ((nperm > 0) & (rsquare > epsilon)) {
    res8 <- permutest.lmodel2(yx, yx.2, b.ols, b.ma, b.rma/ratio, 
                              RMA, ratio, nperm, epsilon)
  }
  if (H == 0) 
    H <- NA
  reg.res <- data.frame(met, res1, res2, res3, res8)
  CI.res <- data.frame(met, res4, res5, res6, res7)
  colnames(reg.res) <- c("Method", "Intercept", "Slope", "Angle (degrees)", 
                         "P-perm (1-tailed)")
  colnames(CI.res) <- c("Method", "2.5%-Intercept", "97.5%-Intercept", 
                        "2.5%-Slope", "97.5%-Slope")
  out <- list(y = y, x = x, regression.results = reg.res, confidence.intervals = CI.res, 
              eigenvalues = c(lambda1, lambda2), H = H, n = n, r = r, 
              rsquare = rsquare, P.param = P.param, theta = theta, 
              nperm = nperm, epsilon = epsilon, info.slope = info.slope, 
              info.CI = info.CI, call = match.call())
  class(out) <- "lmodel2"
  out
}

# Import from lmodel2
MA.reg <-
  function(b.ols, r, rsquare, ybar, xbar, epsilon, ratio)
  {
    if(rsquare > epsilon) {
      d <- (b.ols^2 - rsquare) /(rsquare * b.ols)
      b.ma <- 0.5*(d + sign(r)*sqrt(d^2+4))
      b.ma <- b.ma*ratio
      b0 <- ybar - b.ma*xbar
      angle <- atan(b.ma)*180/pi
    } else { 
      b0 <- NA
      b.ma <- NA
      angle <- NA
    }
    MA.res <- c(b0, b.ma, angle)
  }

# Import from lmodel2
SMA.reg <-
  function(b.ols, r, rsquare, ybar, xbar, epsilon)
  {
    if(rsquare > epsilon) { 
      b.sma <- sign(r) * b.ols/r 
      b0 <- ybar - b.sma*xbar
      angle <- atan(b.sma)*180/pi
    } else { 
      b0 <- NA
      b.sma <- NA
      angle <- NA
    }
    SMA.res <- c(b0, b.sma, angle)
  }

# Import from lmodel2
CLma <-
  function(b.ma, n, r, rsquare, xbar, ybar, cov.mat, t, ratio, epsilon)
  {
    H <- 0
    A <- NA
    info.CI <- 0
    ## Compute eigenvalues by eigen(covariance matrix), as in PCA
    mat.eig <- eigen(cov.mat)
    lambda1 <- mat.eig$values[1]
    lambda2 <- mat.eig$values[2]
    #
    if(rsquare > epsilon) {
      b.ma <- b.ma/ratio
      if(lambda2 <= epsilon) {               # If eigenvalue2 = 0
        b1inf <-  b.ma
        b1sup <- b.ma
      } else {
        if( (lambda1-lambda2) < epsilon) {  # If equal eigenvalues
          tmp <- atan(b.ma)*180/pi
          b1inf <- tan((tmp-90)*pi/180)
          b1sup <- tan((tmp+90)*pi/180)
        } else {
          H <- t^2 / (((lambda1/lambda2)+(lambda2/lambda1)-2)*(n-2))
          if(H >= 1) {                     # If H >= 1
            b1inf <- NA
            b1sup <- NA
          } else {
            A <- sqrt(H/(1-H))
            if((b.ma*A) == -1) {          # angle = -90 deg., b1inf = -infinity
              b1inf <- -Inf
              info.CI <- 1
            } else {
              b1inf <- ratio * (b.ma-A) / (1+b.ma*A)
            }
            if((b.ma*A) == 1) {           # angle = +90 deg., b1sup = +infinity
              b1sup <- Inf
              info.CI <- 1
            } else {
              b1sup <- ratio * (b.ma+A) / (1-b.ma*A)
            }
          }
        }
      }
      if((H == 0) || (H >= 1)) {             # H >= 1
        b0inf <- NA
        b0sup <- NA
      } else {
        if(xbar >= 0) {
          b0inf <- ybar - b1sup*xbar
          b0sup <- ybar - b1inf*xbar
        } else {
          b0inf <- ybar - b1inf*xbar
          b0sup <- ybar - b1sup*xbar
        }
      }
    } else {
      b1inf <- NA
      b1sup <- NA
      b0inf <- NA
      b0sup <- NA
    }
    CL <- c(b0inf, b0sup, b1inf, b1sup, lambda1, lambda2, H, A, info.CI)
  }

# Import lmodel2
CLsma <-
  function(b.sma, n, r, rsquare, xbar, ybar, t, epsilon)
  {
    if(rsquare > epsilon) {
      B <- t^2 * (1-rsquare) / (n-2)
      if(b.sma > 0) {
        b1inf <- b.sma * (sqrt(B+1) - sqrt(B))
        b1sup <- b.sma * (sqrt(B+1) + sqrt(B))
      } else {
        b1inf <- b.sma * (sqrt(B+1) + sqrt(B))
        b1sup <- b.sma * (sqrt(B+1) - sqrt(B))
      }
      if(xbar >= 0) {
        b0inf <- ybar - b1sup*xbar
        b0sup <- ybar - b1inf*xbar
      } else {
        b0inf <- ybar - b1inf*xbar
        b0sup <- ybar - b1sup*xbar
      }      
    } else {
      b1inf <- NA
      b1sup <- NA
      b0inf <- NA
      b0sup <- NA
    }
    CL <- c(b0inf, b0sup, b1inf, b1sup)
  }