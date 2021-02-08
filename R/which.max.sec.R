#' @noRd
which.max.sec <- function(x) {
  
  if(length(x) == 1) {
    
    return(which.max(x))
    
  } else {
    
    y = sort.int(x, decreasing = TRUE, index.return = TRUE)
    return(y$ix[2])
    
  }

}

#' @noRd
max.sec <- function(x) {
  
  if(length(x) == 1) {
    
    return(max(x))
    
  } else {
    
    y = sort(x, decreasing = TRUE)
    return(y[2])
    
  }
  
}