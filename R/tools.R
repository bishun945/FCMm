#' @title Assess the performance of algorithms
#' @name cal.metrics
#' @usage cal.metrics(x,y,name="all",log10=FALSE)
#' @param x True value
#' @param y Estimated value
#' @param name The name of metrics
#' @param log10 Logical. Whether the input x and y should be log10-transformed
#' @export
#' @importFrom stats cor cor.test median na.omit sd
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
  
  .cal.rmse <- function(x,y){
    return(sqrt(sum((x-y)^2,na.rm=T)/length(x)))
  }
  .cal.crmse<-function(x,y){
    return(sqrt(sum((((x-mean(x,na.rm=T))-(y-mean(y,na.rm=T))))^2,na.rm = T)/length(x)))
  }
  .cal.mae<-function(x,y){
    return(result<-sum(abs(y-x),na.rm = T)/length(x))
  }
  .cal.mape<-function(x,y){
    return(result<-sum(abs((y-x)/x),na.rm = T)/length(x)*100)
  }
  .cal.mrpe<-function(x,y){
    return(sum((y-x)/x,na.rm=T)/length(x)*100)
  }
  .cal.umrpe <- function(x,y){
    return(mean((y-x)/(0.5*x+0.5*y)*100))
  }
  .cal.bias<-function(x,y){
    return(sum((y-x),na.rm = T)/length(x))
  }
  .cal.ratio <- function(x,y){
    return(mean(y/x))
  }
  .cal.cv<-function(x,y){
    return(sd(y,na.rm=T) / mean(y, na.rm=T) * 100)
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
  .cal.eta<-function(x,y){
    return(length(which(is.na(y) == FALSE)) / length(y))
  }
  .cal.chi2<-function(x,y){
    return(sum((x-y)^2/y,na.rm = T))
  }

  if(name == "RMSE"){
    result <- .cal.rmse(x,y)
  }else if(name == "CRMSE"){
    result <- .cal.crmse(x,y)
  }else if(name == "MAE"){
    result <- .cal.mae(x,y)
  }else if(name == "MAPE"){
    result <- .cal.mape(x,y)
  }else if(name == "MRPE"){
    result <- .cal.mrpe(x,y)
  }else if(name == "UMRPE"){
    result <- .cal.umrpe(x,y)
  }else if(name == "BIAS"){
    result <- .cal.bias(x,y)
  }else if(name == "RATIO"){
    result <- .cal.ratio(x,y)
  }else if(name == "CV"){
    result <- .cal.cv(x,y)
  }else if(name == "R"){
    result <- .cal.R(x,y)
  }else if(name == "R2"){
    result <- .cal.R2(x,y)
  }else if(name == "eta"){
    result <- .cal.eta(x,y)
  }else if(name == "chi2"){
    result <- .cal.chi2(x,y)
  }
  else if(name == "all"){
    result <- list(RMSE  = .cal.rmse(x,y),
                   CRMSE = .cal.crmse(x,y),
                   MAE   = .cal.mae(x,y),
                   MAPE  = .cal.mape(x,y),
                   UMRPE = .cal.umrpe(x,y),
                   BIAS  = .cal.bias(x,y),
                   RATIO = .cal.ratio(x,y),
                   CV    = .cal.cv(x,y),
                   R     = .cal.R(x,y),
                   R2    = .cal.R2(x,y),
                   eta   = .cal.eta(x,y),
                   chi2  = .cal.chi2(x,y)
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
  c('RMSE','CRMSE','MAE','MAPE',"MRPE","UMRPE",
    'BIAS', 'RATIO', 'CV','R','R2','eta','chi2','all')
}

#' @title Assess the performance of algorithms as vector
#' @name cal.metrics.vector
#' @usage cal.metrics.vector(x,y,name="all",log10=FALSE)
#' @param x True value
#' @param y Estimated value
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
  .cal.mape<-function(x,y){
    return(abs((y-x)/x)*100)
  }
  .cal.mrpe<-function(x,y){
    return((y-x)/x*100)
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
  }else if(name == "MAPE"){
    result <- .cal.mape(x,y)
  }else if(name == "MRPE"){
    result <- .cal.mrpe(x,y)
  }else if(name == "UMRPE"){
    result <- .cal.umrpe(x,y)
  }else if(name == "BIAS"){
    result <- .cal.bias(x,y)
  }else if(name == "RATIO"){
    result <- .cal.ratio(x,y)
  }
  else if(name == "all"){
    result <- list(MAE   = .cal.mae(x,y),
                   MAPE  = .cal.mape(x,y),
                   MRPE  = .cal.mrpe(x,y),
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
  c('MAE','MAPE','MRPE', 'UMRPE',
    'BIAS', 'RATIO', 'all')
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
