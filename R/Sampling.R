#' @name Sampling_via_cluster
#' @title Stratified sampling by clusters or types
#' @param x A vector represents the cluster or type, only support numeric value now.
#' @param num The number of sampling. Default as \code{length(x)}
#' @param replace The way of sampling. See \code{help(sample)} for more details.
#' @param order.value The vector used for ordering. Default is \code{NULL}.
#' @param ... Parameters of function \code{Sampling_by_sort}.
#' @return Result of \code{Sampling_via_cluster} is the index of sampled 
#'   items from the input \code{x}.
#' @family Algorithm assessment
#' @examples 
#' 
#' library(FCMm)
#' x = round(runif(100, 1, 10))
#' table(x)
#' w_sampled = Sampling_via_cluster(x, 20)
#' table(x[w_sampled])
#' 
#' @export
#' 
Sampling_via_cluster <- function(x, num=length(x), 
                                 replace=FALSE, 
                                 order.value= NULL,
                                 ...){
  
  if(num > length(x))
    stop("Sample size is larger than input size.")
  
  if(anyNA(x))
    stop("NA values in input variable.")
  
  x <- as.vector(x)
  types <- sort(unique(x))
  
  types_num_final <- Sampling_each_cluster_num(x, num)
  
  # Start to sample
  ind <- seq(length(x))
  
  result <- NULL
  
  for(i in 1:length(types_num_final)){
    
    x_type <- ind[which(x == types[i])]
    v_type <- order.value[x_type]
    
    if(is.null(order.value)) {
      
      tmp_sample <- sample(x_type, types_num_final[i], replace=replace)
      result <- c(result, tmp_sample)
      
    } else {
      
      tmp_sample <- Sampling_by_sort(v_type, 
                                     types_num_final[i], 
                                     replace=replace,
                                     ...)
      tmp_sample <- x_type[tmp_sample]
      result <- c(result, tmp_sample)
      
    }
    
    
  }
  
  return(result)
  
}

#' @noRd
Sampling_each_cluster_num <- function(x, num) {
  
  # Define the number of each cluster
  types <- sort(unique(x))
  raw_num <- as.numeric(table(x))
  types_num <- as.numeric(table(x)) / length(x) * num
  types_num <- floor(types_num)
  resi <- num - sum(types_num)
  
  if(resi < 0){
    
    stop('Stratifed sampling total number is greater than the sample number.')
    
  }else if(resi == 0){
    
    types_num_final <- types_num
    
  }else if(resi > 0){
    
    resi_ <- resi
    types_num_final <- types_num
    condition1 = !all(types_num_final <= raw_num) | sum(types_num_final) < num
    while(condition1){
      w_candidate <- which(types_num_final < raw_num)
      rand_select <- sample(w_candidate, 1, replace=TRUE)
      if(types_num_final[rand_select] < raw_num[rand_select]){
        types_num_final[rand_select] <- types_num_final[rand_select] + 1
      }
      condition1 = !all(types_num_final <= raw_num) | sum(types_num_final) < num
    }
  }
  
  return(types_num_final)
  
}


#' @rdname Sampling_via_cluster
#' @param log10 Whether to log10-trans the \code{order.value}. Default is \code{FALSE}.
#' @param n_group Number of group to be devided. Default is 3.
#' @export
#' @examples 
#' 
#' set.seed(1)
#' N <- 100
#' x <- stats::rlnorm(N)
#' ind1 <- Sampling_by_sort(x, N/2)
#' ind2 <- sample.int(length(x), N/2)
#' 
#' library(ggplot2)
#' ggplot() + 
#' geom_density(aes(x = x, color = "Total")) + 
#' geom_density(aes(x = x[ind1], color = "ind1"), alpha = 0.8) + 
#' geom_density(aes(x = x[ind2], color = "ind2"), alpha = 0.6) +
#' scale_x_log10()
#' 
#' 
#' 
Sampling_by_sort <- function(x, 
                             num     = length(x), 
                             log10   = FALSE,
                             n_group = 3, 
                             replace = TRUE) {
  
  xlen   <- length(x)
  
  if(log10) {
    x <- log10(x)
  } 
  
  if(num > xlen) {
    
    stop("The sample select number is greater than the x length!")
    
  }
  
  if(xlen <= n_group) {
    
    warning(sprintf("The length of `x` = %s is less than that of `n_group` = %s",
                    xlen, n_group))
    
    return(sample.int(xlen, num, replace = replace))
    
  } else {
    
    z_sort <- sort.int(x, index.return = TRUE)
    group  <- ggplot2::cut_number(z_sort$x, n_group)
    
    while( all(as.numeric(table(group)) / xlen * num < 1) ) {
      if(n_group == 1) break
      n_group <- n_group - 1
    }
      
    group <- ggplot2::cut_number(z_sort$x, n_group)
    
    
    w_res  <- stats::aggregate(z_sort$ix, list(group), function(x) {
      sample.int(length(x), floor(num*length(x)/xlen), replace = replace)
    }, simplify = FALSE) %>% .[, 2] #%>% unlist() %>% as.vector() %>% sort()
    
    result <- NULL
    for(i in 1:length(levels(group))) {
      result <- c(result, (1:xlen)[group %in% levels(group)[i]][w_res[[i]]])
    }
    
    return(result)
    
  }
  
}
