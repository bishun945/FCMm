#' @name cal.metrics
#' @title Assess the performance of algorithms
#' @param x True value OR Actual value
#' @param y Estimated value OR Predict value
#' @param name The name of metrics
#' @param log10 Logical. Whether the input x and y should be log10-transformed
#' @param c.value Compensated value for CMAPE and CMRPE. Default is the mean of x.
#' @export
#' @return Result of \code{cal.metrics} is a list of selected metric values.
#' @family Utils
#' @note (2020-02-09) All functions used log10 transformation was assigned to the 
#'   Key parameter `log10`.
#' @examples 
#' x = runif(20)
#' y = runif(20)
#' result = cal.metrics(x, y)
#' str(result)
#' @importFrom stats cor cor.test median na.omit sd
#' @importFrom stats confint cov lm qt var
#' @importFrom magrittr %>% %<>%
cal.metrics <- function(x, y, 
                        name = "all", 
                        log10 = FALSE, 
                        c.value = mean(x, na.rm=TRUE)){
  
  name <- match.arg(name, cal.metrics.names())
  
  if(log10){
    x <- log10(x)
    y <- log10(y)
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
    # Symmetric series
  }else if(name == "SMAPE"){
    result <- .cal.smape(x,y)
  }else if(name == "SMDAPE"){
    result <- .cal.smdape(x,y)
  }else if(name == "SMRPE"){
    result <- .cal.smrpe(x,y)
  }else if(name == "SMDRPE"){
    result <- .cal.smdrpe(x,y)
    # Compensated series
  }else if(name == "CMAPE"){
    result <- .cal.cmape(x, y, c.value)
  }else if(name == "CMDAPE"){
    result <- .cal.cmdape(x, y, c.value)
  }else if(name == "CMRPE"){
    result <- .cal.cmrpe(x, y, c.value)
  }else if(name == "CMDRPE"){
    result <- .cal.cmdrpe(x, y, c.value)
    # 
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
                   # Symmetric seires
                   SMAPE     = .cal.smape(x,y),
                   SMDAPE    = .cal.smdape(x,y),
                   SMRPE     = .cal.smrpe(x,y),
                   SMDRPE    = .cal.smdrpe(x,y),
                   # Compensated seires
                   CMAPE     = .cal.cmape(x, y, c.value),
                   CMDAPE    = .cal.cmdape(x, y, c.value),
                   CMRPE     = .cal.cmrpe(x, y, c.value),
                   CMDRPE    = .cal.cmdrpe(x, y, c.value),
                   # 
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

# RMSE Series (Root mean squared error)
#' @noRd
.cal.rmse <- function(x,y){
  return(sqrt(mean((x-y)^2, na.rm=TRUE)))
}

#' @noRd
.cal.crmse<-function(x,y){
  return(sqrt(mean((((x-mean(x, na.rm=TRUE))-(y-mean(y, na.rm=TRUE))))^2, na.rm=TRUE)))
}

# MAE Series (Mean absolute error)
#' @noRd
.cal.mae<-function(x,y){
  return(mean(abs(x-y), na.rm=TRUE))
}

#' @noRd
.cal.mdae<-function(x,y){
  return(median(abs(x-y), na.rm=TRUE))
}

#' @noRd
.cal.nmae_sd <- function(x,y){
  return(mean(abs(x-y), na.rm=TRUE)/sd(y, na.rm=TRUE))
}

#' @noRd
.cal.nmdae_sd <- function(x,y){
  return(median(abs(x-y), na.rm=TRUE)/sd(y, na.rm=TRUE))
}

# MAPE Series (Mean absolute percent error)
#' @noRd
.cal.mape<-function(x, y){
  return(mean(abs((x-y)/x), na.rm=TRUE)*100)
}

## Compensated series
#' @noRd
.cal.cmape <- function(x, y, c.value){
  return(mean(       abs(2*(x-y))/(abs(x) + abs(c.value))       , na.rm=TRUE)*100)
}

#' @noRd
.cal.cmdape <- function(x, y, c.value){
  return(median(       abs(2*(x-y))/(abs(x) + abs(c.value))       , na.rm=TRUE)*100)
}

#' @noRd
.cal.cmrpe <- function(x, y, c.value){
  return(mean(       2*(x-y)/(x+c.value)       , na.rm=TRUE)*100)
}

#' @noRd
.cal.cmdrpe <- function(x, y, c.value){
  return(median(       2*(x-y)/(x+c.value)      , na.rm=TRUE)*100)
}

## Symmetric series
#' @noRd
.cal.smape <- function(x, y){
  return(mean(       abs(2*(x-y))/(abs(x)+abs(y))       , na.rm=TRUE)*100)
}

#' @noRd
.cal.smdape <- function(x,y){
  return(median(       abs(2*(x-y))/(abs(x)+abs(y))       , na.rm=TRUE)*100)
}

#' @noRd
.cal.smrpe <- function(x,y){
  return(mean(       (2*(x-y))/(x+y)       , na.rm=TRUE)*100)
}

#' @noRd
.cal.smdrpe <- function(x,y){
  return(median(       (2*(x-y))/(x+y)       , na.rm=TRUE)*100)
}

## 
#' @noRd
.cal.mdape<-function(x,y){
  return(median(abs((x-y)/x), na.rm=TRUE)*100)
}

#' @noRd
.cal.mrpe<-function(x,y){
  return(mean((x-y)/x, na.rm=TRUE)*100)
}

#' @noRd
.cal.mdrpe<-function(x,y){
  return(median((x-y)/x, na.rm=TRUE)*100)
}

#' @noRd
.cal.umrpe <- function(x,y){
  return(mean((x-y)/(0.5*x+0.5*y)*100, na.rm=TRUE))
}

#' @noRd
.cal.umdrpe <- function(x,y){
  return(median((x-y)/(0.5*x+0.5*y)*100, na.rm=TRUE))
}

# Bias Series
#' @noRd
.cal.bias<-function(x,y){
  return(mean((x-y), na.rm=TRUE))
}

#' @noRd
.cal.md_bias<-function(x,y){
  return(median((x-y), na.rm=TRUE))
}

#' @noRd
.cal.ratio <- function(x,y){
  return(mean(x/y, na.rm=TRUE))
}

#' @noRd
.cal.md_ratio <- function(x,y){
  return(median(x/y, na.rm=TRUE))
}

#' @noRd
.cal.cv<-function(x,y){
  return(sd(y, na.rm=TRUE) / mean(y, na.rm=TRUE) * 100)
}

#' @noRd
.cal.R<-function(x,y){
  result <- cor.test(x,y,method = "pearson")
  result <- as.numeric(result$estimate)
  return(result)
}

#' @noRd
.cal.R2<-function(x,y){
  result <- cor.test(x,y,method = "pearson")
  result <- as.numeric(result$estimate) ^ 2
  return(result)
}

# Linear regression relationship Series
#' @noRd
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

#' @noRd
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

#' @noRd
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

#' @noRd
.cal.eta<-function(x,y){
  return(length(which(is.na(y) == FALSE)) / length(y))
}

#' @noRd
.cal.chi2<-function(x,y){
  return(sum((x-y)^2/y, na.rm=TRUE))
}




#' @export
#' @return The call of \code{cal.metrics.names} returns the metrics names in \link{cal.metrics}
#' @rdname cal.metrics
#' 
cal.metrics.names <- function(){
  c('RMSE','CRMSE',
    'MAE','MDAE','NMAE_SD','NMDAE_SD',
    'MAPE', 'SMAPE', 'SMRPE', 'MDAPE',"MRPE",'MDRPE',"UMRPE",'UMDRPE',
    'SMDAPE','SMDRPE',
    'CMAPE', 'CMRPE','CMDAPE','CMDRPE',
    'BIAS', 'MD_BIAS', 'RATIO', 'MD_RATIO',
    'CV','R','R2',
    'SLOPE', 'INTERCEPT', 'R2_SMA',
    'eta','chi2','all')
}