#' @title Assess the performance of algorithms as vector
#' @name cal.metrics.vector
#' @param x True value OR Actual value
#' @param y Estimated value OR Predict value
#' @param name The name of metrics
#' @param log10 Logical. Whether the input x and y should be log10-transformed
#' @param c.value Compensated value for CMAPE and CMRPE. Default is the mean of x.
#' @export
#' @return Result of \code{cal.metrics.vector} is a list of selected metric values.
#' @importFrom stats cor cor.test median na.omit sd
#' @family Utils
#' @examples 
#' x = runif(100)
#' y = runif(100)
#' result = cal.metrics.vector(x,y)
#' 
cal.metrics.vector <- function(x, y, 
                               name = "all", log10 = FALSE,
                               c.value = mean(x, na.rm=TRUE)){
  
  name <- match.arg(name, cal.metrics.vector.names())
  
  if(log10){
    x <- log10(x)
    y <- log10(y)
  }

  if(name == "MAE"){
    result <- .cal.mae.v(x,y)
  }else if(name == "NMAE"){
    result <- .cal.nmae.v(x,y)
  }else if(name == "MAPE"){
    result <- .cal.mape.v(x,y)
    ## Symmetric
  }else if(name == "SMAPE"){
    result <- .cal.smape.v(x,y)
  }else if(name == "SMRPE"){
    result <- .cal.smrpe.v(x,y)
    ## Compensated
  }else if(name == "CMAPE"){
    result <- .cal.cmape.v(x, y, c.value)
  }else if(name == "CMRPE"){
    result <- .cal.cmrpe.v(x, y, c.value)
    ##
  }else if(name == "NMAPE"){
    result <- .cal.nmape.v(x,y)
  }else if(name == "MRPE"){
    result <- .cal.mrpe.v(x,y)
  }else if(name == "NMRPE"){
    result <- .cal.nmrpe.v(x,y)
  }else if(name == "UMRPE"){
    result <- .cal.umrpe.v(x,y)
  }else if(name == "BIAS"){
    result <- .cal.bias.v(x,y)
  }else if(name == "RATIO"){
    result <- .cal.ratio.v(x,y)
  }
  else if(name == "all"){
    result <- list(MAE   = .cal.mae.v(x,y),
                   NMAE  = .cal.nmae.v(x,y),
                   MAPE  = .cal.mape.v(x,y),
                   # symmetric
                   SMAPE = .cal.smape.v(x,y),
                   SMRPE = .cal.smrpe.v(x,y),
                   # compensated
                   CMAPE = .cal.cmape.v(x, y, c.value),
                   CMRPE = .cal.cmrpe.v(x, y, c.value),
                   # 
                   NMAPE = .cal.nmape.v(x,y),
                   MRPE  = .cal.mrpe.v(x,y),
                   NMRPE = .cal.nmrpe.v(x,y),
                   UMRPE = .cal.umrpe.v(x,y),
                   BIAS  = .cal.bias.v(x,y),
                   RATIO = .cal.ratio.v(x,y)
    )
  }
  
  return(result)
  
}

#' @noRd
.cal.mae.v <- function(x,y){
  return(abs(y-x))
}

#' @noRd
.cal.nmae.v <- function(x,y){
  return(abs(y-x) / sd(y, na.rm=TRUE))
}

#' @noRd
.cal.mape.v <- function(x,y){
  return(abs((y-x)/x)*100)
}

## Symmetric sereis
#' @noRd
.cal.smape.v <- function(x,y){
  return(abs(2*(y-x))/(abs(x)+abs(y))*100)
}

#' @noRd
.cal.smrpe.v <- function(x,y){
  return((2*(y-x))/((x)+(y))*100)
}

## Compensated sereis
#' @noRd
.cal.cmape.v <- function(x, y, c.value){
  return(abs(2*(y-x))/(abs(x)+abs(c.value))*100)
}

#' @noRd
.cal.cmrpe.v <- function(x, y, c.value){
  return((2*(y-x))/((x)+(c.value))*100)
}

#' @noRd
.cal.nmape.v <- function(x,y){
  return(abs((y-x)/x) / sd(y, na.rm=TRUE) *100)
}

#' @noRd
.cal.mrpe.v <- function(x,y){
  return((y-x)/x*100)
}

#' @noRd
.cal.nmrpe.v <- function(x,y){
  return((y-x)/x / sd(y, na.rm=TRUE) *100)
}

#' @noRd
.cal.umrpe.v <- function(x,y){
  return((y-x)/(0.5*x+0.5*y)*100)
}

#' @noRd
.cal.bias.v <- function(x,y){
  return(y-x)
}

#' @noRd
.cal.ratio.v <- function(x,y){
  return(y/x)
}


#' @export
#' @return The call of \code{cal.metrics.vector.names} retunrs the metrics 
#'   names in \link{cal.metrics.vector}
#' @rdname cal.metrics.vector
#' 
cal.metrics.vector.names <- function(){
  c('MAE', 'NMAE',
    'MAPE', 'SMAPE', 'NMAPE', 'CMAPE','CMRPE',
    'MRPE', 'SMRPE', 'NMRPE', 'UMRPE',
    'BIAS', 'RATIO', 'all')
}