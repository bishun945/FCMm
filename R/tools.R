#' @title cal.metrics
#' @name cal.metrics
#' @usage cal.metrics(x,y,name="all")
#' @param x True value
#' @param y Estimated value
#' @param name The name of metrics
#' @export
#' @importFrom stats cor cor.test median na.omit sd
cal.metrics <- function(x,y,name="all"){

  .cal.rmse <- function(x,y){
    result <- sqrt(sum((x-y)^2,na.rm=T)/length(x))
    return(result)
  }
  .cal.crmse<-function(x,y){
    result<-sqrt(sum((((x-mean(x,na.rm=T))-(y-mean(y,na.rm=T))))^2,na.rm = T)/length(x))
    return(result)
  }
  .cal.mae1<-function(x,y){
    result<-sum(abs(x-y),na.rm = T)/length(x)
    return(result)
  }
  .cal.mae2<-function(x,y){
    result<-10^(sum(abs((log10(x))-log10(y)),na.rm = T)/length(x))
    return(result)
  }
  .cal.mape<-function(x,y){
    result<-sum(abs((x-y)/y),na.rm = T)/length(x)*100
    return(result)
  }
  .cal.mdape<-function(x,y){
    result<-sum(abs((y-median(x,na.rm = T))/x),na.rm = T)/length(x)*100
    return(result)
  }
  .cal.bias1<-function(x,y){
    result<-sum((y-x),na.rm = T)/length(x)
    return(result)
  }
  .cal.bias2<-function(x,y){
    result<-10^(sum((log10(y)-log10(x)),na.rm = T)/length(x))
    return(result)
  }
  .cal.cv<-function(x,y){
    # result<-sqrt(sum(y-mean(y,na.rm = T),na.rm = T)^2/length(x)) / abs(mean(y,na.rm = T))
    result <- sd(y,na.rm=T) / mean(y, na.rm=T) * 100
    return(result)
  }
  .cal.R<-function(x,y){
    result <- cor.test(x,y,method = "pearson")
    result <- as.numeric(result$estimate)
    return(result)
  }
  .cal.R2<-function(x,y){
    # result<-1-(sum((y-mean(x,na.rm = T))^2,na.rm = T)/sum((x-mean(x,na.rm = T))^2,na.rm = T))
    result <- cor.test(x,y,method = "pearson")
    result <- as.numeric(result$estimate) ^ 2
    return(result)
  }
  .cal.eta<-function(x,y){
    result<-length(which(is.na(y) == FALSE)) / length(y)
    return(result)
  }
  .cal.chi2<-function(x,y){
    result<-sum((x-y)^2/y,na.rm = T)
    return(result)
  }

  if(name == "RMSE"){
    result <- .cal.rmse(x,y)
  }else if(name == "CRMSE"){
    result <- .cal.crmse(x,y)
  }else if(name == "MAE1"){
    result <- .cal.mae1(x,y)
  }else if(name == "MAE2"){
    result <- .cal.mae2(x,y)
  }else if(name == "MAPE"){
    result <- .cal.mape(x,y)
  }else if(name == "MDAPE"){
    result <- .cal.mdape(x,y)
  }else if(name == "BIAS1"){
    result <- .cal.bias1(x,y)
  }else if(name == "BIAS2"){
    result <- .cal.bias2(x,y)
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
    result <- list(RMSE=.cal.rmse(x,y),
                   CRMSE=.cal.crmse(x,y),
                   MAE1=.cal.mae1(x,y),
                   MAE2=.cal.mae2(x,y),
                   MAPE=.cal.mape(x,y),
                   MDAPE=.cal.mdape(x,y),
                   BIAS1=.cal.bias1(x,y),
                   BIAS2=.cal.bias2(x,y),
                   CV=.cal.cv(x,y),
                   R=.cal.R(x,y),
                   R2=.cal.R2(x,y),
                   eta=.cal.eta(x,y),
                   chi2=.cal.chi2(x,y)
    )
  }
  return(result)
}

#' @title .level_to_variable
#' @name .level_to_variable
#' @param dt dataframe
#' @param warn warnning option
#' @export
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
