#' @name Assessment_via_cluster
#' @title Assessment each algorithm for every cluster
#' @param pred prediction of Chla
#' @param meas in-situ measurement of Chla
#' @param memb membership value matrix
#' @param metrics metrics need to be calculated
#' @param log10 Should pred and meas be log10-transformed? (default as \code{FALSE})
#' @param total Whether to calculate summarized metrics (default as \code{TRUE})
#' @param hard.mode If \code{FALSE}, the membership values are used to calculate validation metrics
#' @param cal.precision Whether to calculate precision (only support for vectorized metrics), default as \code{FALSE}
#' @param na.process na.process and choose to statistic na value percent
#' @param plot.col option to plot col result for chosed metrics (default as \code{FALSE})
#' 
#' @note If the \code{cal.precision} is \code{TRUE}, the \code{hard.mode == TRUE} is used. In that case,
#'   mean and sd calculation is conducted for hard mode based on result from \link{cal.metrics.vector}.
#'   
#' @export
#' @return Results of \code{Assessment_via_cluster} are returned as a list including:
#' \item{Metrcs}{A list of the selected metric values for all algorithms. \code{Valid_percent} 
#'   would be included if \code{na.process} are set as \code{TRUE}}
#' \item{res_plot}{Bar plots by using ggplot function for metrics value at every cluster.}
#' \item{res_plot_dt}{Dataframe for plotting \code{res_plot}. I just keep it in case of plotting other types}
#' \item{res_plot_facet}{\code{res_plot} added on \code{facet_wrap}.}
#' \item{input}{input parameters of \link{Assessment_via_cluster}}
#' 
#' @family Algorithm assessment
#' 
#' @import ggplot2  
#' @importFrom stats runif quantile
#' @importFrom reshape2 melt
#' 
#' @examples 
#' library(FCMm) 
#' library(ggplot2) 
#' library(magrittr)
#' library(stringr)
#' data("Nechad2015")
#' x <- Nechad2015[,3:11]
#' wv <- gsub("X","",names(x)) %>% as.numeric
#' set.seed(1234)
#' w <- sample(1:nrow(x), 100) 
#' x <- x[w, ]
#' names(x) <- wv
#' nb = 4 # Obtained from the vignette "Cluster a new dataset by FCMm"
#' set.seed(1234)
#' FD <- FuzzifierDetermination(x, wv, stand=FALSE)
#' result <- FCM.new(FD, nb, fast.mode = TRUE)
#' p.spec <- plot_spec(result, show.stand=TRUE)
#' print(p.spec$p.cluster.spec)
#' Chla <- Nechad2015$X.Chl_a..ug.L.[w]
#' Chla[Chla >= 999] <- NA
#' dt_Chla <- run_all_Chla_algorithms(x) %>% as.data.frame
#' dt_Chla <- data.frame(Chla_true = Chla, 
#' BR_Gil10 = dt_Chla$BR_Gil10, 
#' OC4_OLCI = dt_Chla$OC4_OLCI, 
#' OCI_Hu12 = dt_Chla$OCI_Hu12, 
#' NDCI_Mi12= dt_Chla$NDCI_Mi12) %>% round(3)
#' w = which(!is.na(dt_Chla$Chla_true))
#' dt_Chla = dt_Chla[w,]
#' memb = result$res.FCM$u[w,] %>% round(4)
#' Asses_soft <- Assessment_via_cluster(pred = dt_Chla[,-1], 
#' meas = dt_Chla[,1], memb = memb, log10 = TRUE, hard.mode = FALSE, 
#' na.process = TRUE, plot.col = TRUE)
#' Asses_soft$res_plot_facet
#' knitr::kable(Asses_soft$MAE %>% round(3))
#' knitr::kable(Asses_soft$MAPE %>% round(2))
#' 
#' 
Assessment_via_cluster <- function(pred, meas, memb,
                                   metrics = c('MAE','MAPE'),
                                   log10 = FALSE, 
                                   total = TRUE,
                                   hard.mode= TRUE,
                                   cal.precision = FALSE,
                                   na.process = FALSE,
                                   plot.col = FALSE){
  
  pred <- as.data.frame(pred)
  meas <- as.data.frame(meas)
  
  if(nrow(pred) != nrow(meas) | nrow(pred) != nrow(memb))
    stop('Rows of input are different!')
  
  if(!na.process){
    if(anyNA(pred) | anyNA(meas) | anyNA(memb)){
      stop('Not choose to process NA values. However, predicted, measured or membership values including NA values!')
    }
  }else{
    if(anyNA(meas)){
      stop('Choose to process NA values. But measured values including NA values!')
    }
    if(length(which(meas == 0)) > 0){
      stop('Choose to process NA values. But measured values including ZERO values!')
    }
  }
  
  if(hard.mode == TRUE){
    for(i in 1:length(metrics))
      metrics[i] <- match.arg(metrics[i], cal.metrics.names())
  }else{
    for(i in 1:length(metrics))
      metrics[i] <- match.arg(metrics[i], cal.metrics.vector.names())
  }
  
  if(cal.precision == TRUE){ # precision calculation is only for vectorized metrics
    for(i in 1:length(metrics)){
      metrics[i] <- match.arg(metrics[i], cal.metrics.vector.names())
    }
    hard.mode = TRUE # default to hard.mode for mean and sd calculation in each cluster
  }
  
  # generate the output dataframe  
  model_names <- colnames(pred)
  cluster_names <- colnames(memb)
  cluster_crisp <- apply(memb,1,which.max)
  
  validation <- matrix(data=0,
                       nrow=length(cluster_names), 
                       ncol=length(model_names)) %>% as.data.frame()
  colnames(validation) <- model_names
  rownames(validation) <- cluster_names
  
  # output is a list, create it!
  result <- list()
  if(cal.precision == TRUE){
    
    for(i in 1:length(metrics)){
      result[[metrics[i]]] <- validation
      precision_name <- paste0(metrics[i], '_p')
      result[[precision_name]] <- validation
    }
    
  }else{ # no need for precision
    
    for(i in 1:length(metrics))
      result[[metrics[i]]] <- validation
  
  }
  
  
  if(na.process){
    result[["Valid_percent"]] <- validation
  }
  
  # strat loop via model and cluster
  for(model in 1:ncol(pred)){
    for(cluster in 1:ncol(memb)){
      
      w <- which(cluster_crisp == cluster)
      x <- meas[w,] # true
      y <- pred[w, model] # pred
      num_raw <- length(x)
      
      if(na.process){
        w_finite <- which(is.finite(y) & y > 0)
        x <- x[w_finite]
        y <- y[w_finite]
        num_new <- length(w_finite)
        if(num_new > num_raw)
          stop("Error! The subseted sample number is smaller the raw.")
      }
      
      # calculate precision
      if(cal.precision == TRUE){
      
        for(metric in metrics){
          precision_name <- paste0(metric, '_p')
          metric_value <- cal.metrics.vector(x,y,metric,log10=log10)
          if(log10){
            x_ = log10(x)
          }else{
            x_ = x
          }
          # quant <- quantile(metric_value, c(0.001, 0.999))
          # metric_value_ = metric_value[metric_value > quant[1] & metric_value < quant[2]]
          metric_value_ = metric_value # for cluster specific metrics, their point numbers are too small to limit
          result[[metric]][cluster, model] <- mean(metric_value_, na.rm=TRUE)
          if(length(metric_value) == 0){
            result[[precision_name]][cluster, model] <- NA
          }else{
            if(median(metric_value, na.rm=TRUE) / median(x_, na.rm=TRUE) > 10){
              # if the error are extremely biased, for instance, 10 fold than measurement,
              #   the precision (standard deviation) will reduce since they are uncorrectly high.
              #   Given that, the precision for this condition is assigned to NA values.
              #   Same as total calculation.
              result[[precision_name]][cluster, model] <- NA
            }else{
              result[[precision_name]][cluster, model] <- sd(metric_value_, na.rm=TRUE)
            }
          }
        }
           
      }else{ # no need for precision
        
        for(metric in metrics){
          result[[metric]][cluster, model] <- cal.metrics(x,y,metric,log10=log10)
        }
        
      }
      
      if(na.process){
        result[["Valid_percent"]][cluster, model] <- num_new / num_raw * 100
      }
    }
    
    if(total == TRUE){
      
      x <- meas[, 1]
      y <- pred[, model]
      num_raw <- length(x)
      
      if(na.process){
        w_finite <- which(is.finite(y) & y > 0)
        x <- x[w_finite]
        y <- y[w_finite]
        num_new <- length(w_finite)
      }
      
      for(metric in metrics){
        
        if(rownames(result[[metric]])[nrow(result[[metric]])] != 'SUM'){
          result[[metric]] %<>% rbind(.,NA)
          rownames(result[[metric]])[nrow(result[[metric]])] <- 'SUM'
        }
        
        # calculate precision for total
        if(cal.precision == TRUE){
          
          precision_name <- paste0(metric, '_p')
          if(rownames(result[[precision_name]])[nrow(result[[precision_name]])] != 'SUM'){
            result[[precision_name]] %<>% rbind(.,NA)
            rownames(result[[precision_name]])[nrow(result[[precision_name]])] <- 'SUM'
          }
          
          for(metric in metrics){
            metric_value <- cal.metrics.vector(x,y,metric,log10=log10)
            # Fot total calculation, it is unfair when outliers are included.
            # I add a quantile calculation to limit the statistic range which is betwene 0.1% and 99.9%
            quant <- quantile(metric_value, c(0.001, 0.999))
            metric_value_ = metric_value[metric_value > quant[1] & metric_value < quant[2]]
            result[[metric]]['SUM', model] <- mean(metric_value_, na.rm=TRUE)
            precision_name <- paste0(metric, '_p')
            if(log10){
              x_ = log10(x)
            }else{
              x_ = x
            }
            if(median(metric_value, na.rm=TRUE) / median(x_, na.rm=TRUE) > 10){
              result[[precision_name]]['SUM', model] <- NA
            }else{
              result[[precision_name]]['SUM', model] <- sd(metric_value_, na.rm=TRUE)
            }
            
          }
          
        }else{ # no need for precision  total
          
          result[[metric]]['SUM', model] <- cal.metrics(x,y,metric,log10=log10)  
          
        }
        
      }
      
      if(na.process){
        result[["Valid_percent"]]['SUM', model] <- num_new / num_raw * 100
      }
    }
  }
  
  
  
  
  # Fuzzy mode
  if(hard.mode == FALSE){
    
    result_fz <- result
    
    for(i in 1:ncol(memb)){
      for(j in 1:ncol(pred)){
        for(metric in metrics){
          
          x <- meas[, 1]
          y <- pred[, j]
          
          w_cluster <- which(cluster_crisp == i)
          
          if(all(is.na(y[w_cluster]))){ # all prediction from j model in cluster i are NA values. Return metrics as NA
            
            result_fz[[metric]][i,j] <- NA
            
          }else{
            
            if(na.process){
              w <- which(is.finite(y))
              x <- x[w]
              y <- y[w]
              memb_ <- memb[w,]
            }
            
            if(dim(memb_)[1] != length(x))
              stop("The fuzzy metrics are calculated with different rows from memb and pred")
            
            Er <- cal.metrics.vector(x,y,metric,log10)
            result_fz[[metric]][i,j] <- sum(memb_[,i] * Er, na.rm=TRUE) / sum(memb_[,i], na.rm=TRUE)
            
          } # ENDIF
    
        }
      }
    }
    
    result <- result_fz
    
  }
  
  # Plot work
  if(plot.col == TRUE){
    
    res_plot <- list()
    res_plot_dt <- list()
    
    num.model <- ncol(pred)
    set.seed(1234)
    ind = runif(num.model) %>% sort.int(., index.return=TRUE) %>% .$ix
    cp = Spectral(num.model)[ind]

    for(metric in names(result)){
      
      tmp <- result[[metric]]
      
      tmp <- cbind(x=rownames(tmp), tmp)
      
      tmp <- melt(tmp, id="x", variable.name='Models')
      
      p <- ggplot(tmp) + 
        geom_col(aes(x=x,y=value,group=Models,fill=Models),
                 position="dodge") + 
        scale_fill_manual(values=cp) + 
        labs(y=metric) + 
        theme_bw()
      
      res_plot[[metric]] <- p
      
      names(tmp)[3] <- metric
      res_plot_dt[[metric]] <- tmp
      
    }
    
    result$res_plot <- res_plot
    result$res_plot_dt <- res_plot_dt
    
    tmp <- data.frame()
    for(i in 1:length(res_plot_dt)){
      if(i == 1){
        tmp = res_plot_dt[[1]]
      }else{
        tmp <- data.frame(tmp, res_plot_dt[[i]][,3])
        names(tmp)[ncol(tmp)] <- names(res_plot_dt[[i]])[3]
      }
    }
    tmpp <- melt(tmp, id=c("x","Models"), variable.name="Metric")
    
    result$res_plot_facet <- 
      ggplot(tmpp) + 
      geom_col(aes(x=x,y=value,group=Models,fill=Models), position='dodge') + 
      scale_fill_manual(values=cp) + 
      labs(x=NULL, y=NULL) + 
      theme_bw() + 
      facet_wrap(~Metric, scales="free_y")
    
  }else{
    
    result$res_plot <- NULL
    result$res_plot_dt <- NULL
    result$res_plot_facet <- NULL
    
  }
  
  result$input <- list(
    
    pred           = pred,
    meas           = meas,
    memb           = memb,
    metrics        = metrics,
    log10          = log10,
    total          = total,
    hard.mode      = hard.mode,
    cal.precision  = cal.precision,
    na.process     = na.process
    
  )
  
  return(result)
  
}



#' @name Score_algorithms_interval
#' @title Score a vector of error metrics for algorithms by the interval
#' @param x Input vector of error metrics (NA values are allowed)
#' @param trim whether to run a trim process to calculate mean and standard deviation of 
#'   input vector x (Default as \code{FALSE})
#' @param reward.punishment Whether to conduct the reward and punishment mechanism in scoring
#'   system (Default as \code{TRUE})
#' @param decreasing the order of the good metric to be evaluated. For instance, MAE should use
#'   \code{decreasing = TRUE} (Default) since the algorithm performs better when MAE becomes smaller. 
#'   However, when comes to \code{Rsquare} from linear regression (maximum is 1), it should be
#'   \code{FALSE}
#' @param hundred.percent A variable constrain for metrics that the maximun of input \code{x} 
#'   should not be greater than 100.
#' @export
#' 
#' @return Results of \code{Score_algorithms_interval()} are returned as a list including:
#' \item{p}{A ggplot list of the scoring result.}
#' \item{score}{The final score from by the interval score.}
#' \item{u}{Trimmed mean of input x with \code{NA} values removed.}
#' \item{bds}{Up and low boundaries for determining scores.}
#' \item{x}{The input x.}
#' 
#' @family Algorithm assessment
#' 
#' @importFrom stats qt sd
#' @import ggplot2
#' 
#' @examples 
#' set.seed(1234)
#' x = runif(10)
#' result = Score_algorithms_interval(x)
#' 

Score_algorithms_interval <- function(x, trim=FALSE, reward.punishment=TRUE, 
                             decreasing=TRUE, hundred.percent=FALSE
                             ){
  
  x <- as.numeric(x)
  
  if(length(x) <= 4){
    stop("Candidate number are less equal than 4, no need use this assessment system.")
  }else{
    if(anyNA(x)){
      n <- length(which(is.na(x) == FALSE))
      n_trim <- ceiling(n/2*0.1)
    }else{
      n <- length(x)
      n_trim <- ceiling(n/2*0.1)
    }
  }
  
  if(trim==TRUE){
    x_ = sort(x)[(n_trim+1):(n-n_trim)]
    u <- mean(x_, na.rm=TRUE)
    sig <- sd(x_, na.rm=TRUE)
  }else{
    u <- mean(x, na.rm=TRUE)
    sig <- sd(x, na.rm=TRUE)
  }
  
  tCL <- qt(0.975,length(x)-1) # 95% Student distribution
  
  if(hundred.percent == FALSE){
    UB = u + tCL*sig/sqrt(n)
    LB = u - tCL*sig/sqrt(n)
    bds = c(UB,LB)
  }else{
    bds = c(50,60,75,95)
    decreasing = FALSE
    reward.punishment = FALSE
  }
  
  tmp = data.frame(x=seq(1,length(x)), y=x)
  
  p <- ggplot() + 
    geom_point(data=tmp, aes(x,y)) + 
    geom_hline(yintercept=u,color='blue') + 
    geom_hline(yintercept=bds,color='blue',linetype=2) + 
    scale_x_continuous(breaks=seq(length(x))) +
    labs(x='Model ID', y='Error metric') + 
    theme_bw()
  
  score <- x * 0
  
  if(hundred.percent == TRUE){
    if(max(x) > 100 | min(x) < 0){
      stop("The input x has a number larger than 100 in Valid_percent metric!")
    }
    score[which(x >= 95)]          = 0
    score[which(x >= 75 & x < 95)] = -0.5 * 2
    score[which(x >= 60 & x < 75)] = -1   * 2   
    score[which(x >= 50 & x < 60)] = -1.5 * 2
    score[which(x <  50)]          = -2   * 2
    
  }else{
    
    if(decreasing == TRUE){ # The smaller value, the better performance
      score[which(x >= UB)] <- 0
      score[which(x < UB & x >= LB)] <- 1
      score[which(x < LB)] <- 2
      # Reward and punishment mechanism for worst and best model
      if(reward.punishment == TRUE){
        score[which.max(x)] <- score[which.max(x)] - 0.5
        score[which.min(x)] <- score[which.min(x)] + 0.5
      }    
    }else{
      score[which(x >= UB)] <- 2
      score[which(x < UB & x >= LB)] <- 1
      score[which(x < LB)] <- 0
      # Reward and punishment mechanism for worst and best model
      if(reward.punishment == TRUE){
        score[which.max(x)] <- score[which.max(x)] + 0.5
        score[which.min(x)] <- score[which.min(x)] - 0.5
      }  
    }
    
    # punish the NA results with -1 score
    if(hundred.percent == TRUE){
      score[which(is.na(x) == TRUE)] <- -2
    }else{
      score[which(is.na(x) == TRUE)] <- -1      
    }
    
    
  }
  
  result <- list(p     = p,
                 score = score,
                 u     = u,
                 bds   = bds,
                 x     = x)
  
  return(result)
}


#' @name Score_algorithms_sort
#' @title Score a vector of error metrics for algorithms by the sort
#' 
#' @param x Input vector of error metrics (NA values are allowed)
#' @param decreasing the order of the good metric to be evaluated. For instance, MAE should use
#'   \code{decreasing = TRUE} (Default) since the algorithm performs better when MAE becomes smaller. 
#'   However, when comes to \code{Rsquare} from linear regression (maximum is 1), it should be
#'   \code{FALSE}
#'   
#' @export
#' 
#' @return Results of \code{Score_algorithms_sort()} are returned as a vector presenting score values.
#' 
#' @family Algorithm assessment
#' 
#' @examples 
#' set.seed(1234)
#' x = runif(10)
#' result = Score_algorithms_sort(x)
Score_algorithms_sort <- function(x, decreasing = TRUE){
  
  x <- as.numeric(x)
  
  # decreasing == TRUE for smaller metrics that are better
  if(decreasing == TRUE){
    na.last = FALSE
  }else{
    na.last = TRUE
  }
  
  r <- sort.int(x, decreasing=decreasing, index.return = TRUE, na.last = na.last)
  
  score <- NULL

  for(i in 1:length(x)){
    if(is.na(x[i])){
      score[i] = 0
    }else{
      score[i] = which(x[i] == r$x)
    }
  }
  
  return(score)

}



#' @name Sampling_via_cluster
#' @title Stratified sampling by clusters or types
#' @param x A vector represents the cluster or type, only support numeric value now.
#' @param num The number of sampling.
#' @param replace The way of sampling. See \code{help(sample)} for more details.
#' @export
#' @return Result of \code{Sampling_via_cluster} is the index of sampled 
#'   items from the input \code{x}.
#' @family Algorithm assessment
#' @examples 
#' set.seed(1234)
#' x = round(runif(100,1,10))
#' table(x)
#' w_sampled = Sampling_via_cluster(x, 20)
#' table(x[w_sampled])
#' 
Sampling_via_cluster <- function(x, num, replace=FALSE){
  
  if(num > length(x))
    stop("Sample size is larger than input size.")
  
  if(anyNA(x))
    stop("NA values in input variable.")
  
  x <- as.numeric(x) 
  
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
  
  # Start to sample
  ind <- seq(length(x))
  
  result <- NULL
  
  for(i in 1:length(types_num_final)){
    
    x_type <- ind[which(x == types[i])]
    tmp_sample <- sample(x_type, types_num_final[i], replace=replace)
    result <- c(result, tmp_sample)
    
  }
  
  return(result)
  
}


#' @name Getting_Asses_results
#' @title Get the result of function Assessment_via_cluster
#' @description This function mainly use function \link{Assessment_via_cluster} to get
#'   assessments both from fuzzy and hard mode. Specifically, it will return the accuracy and precision of 
#'   \code{MAE},\code{CMAPE},\code{BIAS}, and \code{CMRPE} which would be seemed as the input value of 
#'   function \link{Scoring_system}.
#' @param sample.size Sample size. This supports a bootstrap way to run the function 
#'   \link{Assessment_via_cluster}. The number should not be larger than the row number 
#'   of pred or so.
#' @param replace Logical, replace, default as \code{FALSE}
#' @param pred Prediction matrix or data.frame
#' @param meas Measured (actual) matrix or data.frame
#' @param memb Membership matrix
#' @param cluster Cluster vector. Could be calculated by the parameter \code{memb}. Will be deprecated later.
#' @param seed Seed number for fixing the random process. See \code{help(set.seed)} for more details.
#' @note The row number of \code{pred}, \code{meas}, \code{memb}, and \code{cluster} should be the same. 
#'   This function is designed for bootstrapping process to get Chla algorithms assessment. Therefore, 
#'   parameters of \link{Assessment_via_cluster} is set as fixed such as \code{log10 = TRUE}, 
#'   \code{na.process = TRUE}. Given that, I will not export this function in latter to avoid confuses.
#' @export
#' @return A list containing fuzzy and hard results from \link{Assessment_via_cluster}
#' @family Algorithm assessment
#' 
#' @examples 
#' library(FCMm) 
#' library(ggplot2) 
#' library(magrittr)
#' library(stringr)
#' data("Nechad2015")
#' x <- Nechad2015[,3:11]
#' wv <- gsub("X","",names(x)) %>% as.numeric
#' set.seed(1234)
#' w <- sample(1:nrow(x), 100)
#' x <- x[w, ]
#' names(x) <- wv
#' nb = 4 # Obtained from the vignette "Cluster a new dataset by FCMm"
#' set.seed(1234)
#' FD <- FuzzifierDetermination(x, wv, stand=FALSE)
#' result <- FCM.new(FD, nb, fast.mode = TRUE)
#' p.spec <- plot_spec(result, show.stand=TRUE)
#' print(p.spec$p.cluster.spec)
#' Chla <- Nechad2015$X.Chl_a..ug.L.[w]
#' Chla[Chla >= 999] <- NA
#' dt_Chla <- run_all_Chla_algorithms(x) %>% as.data.frame
#' dt_Chla <- data.frame(Chla_true = Chla,
#' BR_Gil10 = dt_Chla$BR_Gil10, 
#' OC4_OLCI = dt_Chla$OC4_OLCI, 
#' OCI_Hu12 = dt_Chla$OCI_Hu12, 
#' NDCI_Mi12= dt_Chla$NDCI_Mi12) %>% round(3)
#' w = which(!is.na(dt_Chla$Chla_true))
#' dt_Chla = dt_Chla[w,]
#' memb = result$res.FCM$u[w,] %>% round(4)
#' cluster =  result$res.FCM$cluster[w]
#' Asses_results <- Getting_Asses_results(sample.size=20, 
#' pred = dt_Chla[,-1], meas = data.frame(dt_Chla[,1]), memb = memb, 
#' cluster = cluster)
#' 
Getting_Asses_results <- function(sample.size, replace = FALSE,
                                  pred, meas, memb,
                                  cluster = apply(memb, 1, which.max), 
                                  seed = NULL){
  
  meas <- data.frame(meas)
  
  if(sample.size > nrow(pred))
    stop("Enter a smaller sample size to subset the data set.")
  
  if(var(c(nrow(pred), nrow(meas), nrow(memb), nrow(cluster))) != 0)
    stop("Row numbers of pred, meas, and memb are different.")
  
  # Stratified sampling by cluster
  if(is.null(seed)){
    w <- Sampling_via_cluster(cluster, num=sample.size, replace=replace)
  }else{
    set.seed(as.numeric(seed))
    w <- Sampling_via_cluster(cluster, num=sample.size, replace=replace)
  }
  
  w <- sort(w)
  
  pred_ <- pred[w,]
  meas_ <- meas[w,]
  memb_ <- memb[w,]
  cluster_ <- cluster[w]
  
  Asses_fz <-  Assessment_via_cluster(pred=pred_,
                                      meas=meas_,
                                      memb=memb_,
                                      metrics = c("MAE","CMAPE","BIAS",'CMRPE'),
                                      log10 = TRUE,
                                      hard.mode = FALSE,
                                      na.process = TRUE,
                                      plot.col = FALSE)
  
  Asses_p <-   Assessment_via_cluster(pred=pred_,
                                      meas=meas_,
                                      memb=memb_,
                                      metrics = c("MAE","CMAPE","BIAS",'CMRPE'),
                                      log10 = TRUE,
                                      hard.mode = FALSE,
                                      cal.precision = TRUE,
                                      na.process = TRUE,
                                      plot.col = FALSE)
  
  result <- list(
                  # Asses_hd=Asses_hd, 
                  Asses_fz = Asses_fz,
                  Asses_p  = Asses_p 
                  )
  
  return(result)
}


#' @name Scoring_system
#' 
#' @title The main function for algorithms scoring based on accuracy, precision, and effectiveness.
#' 
#' @param Inputs The list returned form function \link{Getting_Asses_results}
#' @param method The method selected to score algorithms: 
#'   \itemize{
#'     \item \code{sort-based} (default) which is scored by the sort of accuracy and precision 
#'       metrics (see more in \link{Score_algorithms_sort}). 
#'     \item \code{interval-based} which is relatively scored by the interval of accuracy and 
#'       precision (used by Brewin et al. (2015) and Neil et al. (2019)). 
#'       See more in \link{Score_algorithms_interval}).
#'   }
#' @param param_sort The parameters of function \link{Score_algorithms_sort}
#' @param param_interval The parameters of function \link{Score_algorithms_interval}
#' @param remove.negative Option to replace the negative score as zero (default as \code{FALSE})
#' @param Times Parameter of \code{Scoring_system_bootstrap}. The bootstrap time for running 
#'   \code{Scoring_system} (default as \code{1000})
#' @param replace Parameter of \code{Scoring_system_bootstrap}. The sample method for bootstrap running.
#'   Default as \code{TRUE}. See more details in \link{sample}.
#' 
#' @export
#' 
#' @return The result of \code{Scoring_system} are including:
#' \itemize{
#'   \item \strong{Total_score} Data.frame of final score result with algorithm as column and cluster as row.
#'   \item \strong{Accuracy} Data.frame of \code{Accuracy} score with algorithm as column and cluster as row.
#'   \item \strong{Precision} Data.frame of \code{Precision} score with algorithm as column and cluster as row.
#'   \item \strong{Effectiveness} Data.frame of \code{Effectiveness} score with algorithm as column and cluster as row.
#'   \item \strong{Accuracy.list} List including data.frames of used \code{Accuracy} metrics.
#'   \item \strong{Precision.list} List including data.frames of used \code{Precision} metrics.
#'   \item \strong{Total_score.melt} Melted data.frame of \strong{Total_score} for plotting.
#'   \item \strong{opt_algorithm} The optimal algorithm names for each cluster.
#'   \item \strong{Inputs} Inputs of this function.
#' } 
#' 
#' @note 
#'   \code{Scoring_system_bootstrap} is the bootstrap mode of \code{Scoring_system} which is useful when
#'     the outcome is unstable for large number of samples. The default boostrap time in \code{Scoring_system_bootstrap}
#'     is set as \code{1000} and the result of it is the list of several aggregated data.frames and standard deviations.
#'   
#' @details  
#' The \code{Accuracy} and \code{Precision} is newly defined in \code{FCMm} package (referred by 
#'   Hooker et al. (2005)):
#' \itemize{
#'   \item \code{Accuracy} is the estimation of how close the result of the experiment is to the 
#'     true value.
#'   \item \code{Precision} is the estimation of how excatly the result is determined independently 
#'     of any true value.
#' }
#' In other words, \code{Accuracy} is telling a story truthfully and precision is how similarly
#'   the story is represented over and over again.
#' Here we use AE, a vector for each sample, for instance:
#' \itemize{
#'   \item \code{Accuracy} is the aggregation (no matter mean or median, in fuzzy calculation process), 
#'     we use mean to some extent. 
#'   \item \code{Precision} is actually the \strong{stability} of AE (reproducebility) which means the error
#'     produced by the algorithm is under certain control.
#' }
#' Finally, the function will multiply the total score (\code{Accuracy} + \code{Precision}) by the 
#'   effectiveness (i.e., Valid_percent returned by \link{Assessment_via_cluster}).
#' 
#' @family Algorithm assessment
#' 
#' @importFrom reshape2 melt
#' @importFrom magrittr %>%
#' 
#' @references
#' \itemize{
#'   \item Hooker S B. Second SeaWiFS HPLC analysis round-robin experiment (SeaHARRE-2)[M]. 
#'     National Aeronautics and Space Administration, Goddard Space Flight Center, 2005. 
#'   \item Seegers B N, Stumpf R P, Schaeffer B A, et al. Performance metrics for the assessment 
#'     of satellite data products: an ocean color case study[J]. Optics express, 2018, 26(6): 7404-7422.
#'   \item Neil C, Spyrakos E, Hunter P D, et al. A global approach for chlorophyll-a retrieval 
#'     across optically complex inland waters based on optical water types[J]. Remote Sensing of 
#'     Environment, 2019, 229: 159-178.
#'   \item Brewin R J W, Sathyendranath S, Müller D, et al. The Ocean Colour Climate Change 
#'     Initiative: III. A round-robin comparison on in-water bio-optical algorithms[J]. Remote Sensing of 
#'     Environment, 2015, 162: 271-294.
#'   \item Moore T S, Dowell M D, Bradt S, et al. An optical water type framework for selecting and 
#'     blending retrievals from bio-optical algorithms in lakes and coastal waters[J]. Remote sensing of 
#'     environment, 2014, 143: 97-111.
#' } 
#' 
#' @examples 
#' library(FCMm) 
#' library(ggplot2) 
#' library(magrittr)
#' library(stringr)
#' data("Nechad2015")
#' x <- Nechad2015[,3:11]
#' wv <- gsub("X","",names(x)) %>% as.numeric
#' set.seed(1234)
#' w <- sample(1:nrow(x), 100)
#' x <- x[w, ]
#' names(x) <- wv
#' nb = 4 # Obtained from the vignette "Cluster a new dataset by FCMm"
#' set.seed(1234)
#' FD <- FuzzifierDetermination(x, wv, stand = FALSE)
#' result <- FCM.new(FD, nb, fast.mode = TRUE)
#' p.spec <- plot_spec(result, show.stand = TRUE, show.ribbon = TRUE)
#' print(p.spec$p.cluster.spec)
#' Chla <- Nechad2015$X.Chl_a..ug.L.[w]
#' Chla[Chla >= 999] <- NA
#' dt_Chla <- run_all_Chla_algorithms(x) %>% as.data.frame
#' dt_Chla <- data.frame(Chla_true = Chla, 
#' BR_Gil10 = dt_Chla$BR_Gil10, 
#' OC4_OLCI = dt_Chla$OC4_OLCI, 
#' OCI_Hu12 = dt_Chla$OCI_Hu12, 
#' NDCI_Mi12= dt_Chla$NDCI_Mi12) %>% round(3)
#' w = which(!is.na(dt_Chla$Chla_true))
#' dt_Chla = dt_Chla[w,]
#' memb = result$res.FCM$u[w,] %>% round(4)
#' cluster =  result$res.FCM$cluster[w]
#' Asses_results <- Getting_Asses_results(sample.size=20, pred = dt_Chla[,-1], 
#' meas = data.frame(dt_Chla[,1]), memb = memb)
#' Score = Scoring_system(Asses_results)
#' 
Scoring_system <- function(Inputs, 
                           method = 'sort-based',
                           param_sort = list(decreasing = TRUE),
                           param_interval = list(trim=FALSE, reward.punishment=TRUE,
                                             decreasing=TRUE, hundred.percent=FALSE),
                           remove.negative = FALSE){
  
  Asses_fz <- Inputs$Asses_fz
  Asses_p  <- Inputs$Asses_p # If on precision, the mode is hard
  
  method = match.arg(method, c('sort-based', 'interval-based'))
  if(method == 'sort-based'){
    Score_algorithms <- function(x){
      return(Score_algorithms_sort(x, decreasing = param_sort$decreasing))
    }
  }else if(method == 'interval-based'){
    Score_algorithms <- function(x){
      a = param_interval
      r = Score_algorithms_interval(x, trim=a$trim, reward.punishment = a$reward.punishment,
                                decreasing=a$decreasing, hundred.percent = a$hundred.percent)
      return(r$score)
    }
  }else{
    stop('Method selection error.')
  }

  # Note 2020-03-06:

  #
  Accuracy_CMRPE <- Accuracy_BIAS <- Accuracy_MAE <- Accuracy_CMAPE <- Asses_fz$MAE * NA
  for(i in 1:nrow(Accuracy_MAE)){
    Accuracy_MAE[i,]   <- Score_algorithms(Asses_fz$MAE[i,])
    Accuracy_BIAS[i,]  <- Score_algorithms(abs(Asses_fz$BIAS[i,]))
    Accuracy_CMAPE[i,] <- Score_algorithms(Asses_fz$CMAPE[i,])
    Accuracy_CMRPE[i,] <- Score_algorithms(abs(Asses_fz$CMRPE[i,]))
  }
  Accuracy <- Accuracy_MAE + Accuracy_BIAS + Accuracy_CMRPE + Accuracy_CMAPE
  
  
  # Precision part
  Precision_CMRPE <- Precision_BIAS <- Precision_MAE <- Precision_CMAPE <- Asses_p$MAE_p * NA
  for(i in 1:nrow(Precision_MAE)){
    Precision_MAE[i,]   <- Score_algorithms(Asses_p$MAE_p[i,])
    Precision_BIAS[i,]  <- Score_algorithms(abs(Asses_p$BIAS_p[i,]))
    Precision_CMAPE[i,] <- Score_algorithms(Asses_p$CMAPE_p[i,])
    Precision_CMRPE[i,] <- Score_algorithms(abs(Asses_p$CMRPE_p[i,]))
  }
  Precision <- Precision_MAE + Precision_BIAS + Precision_CMRPE + Precision_CMAPE
  
  Total_score <- (Accuracy + Precision) * Asses_fz$Valid_percent / 100
  
  if(remove.negative == TRUE){
    w = which(Total_score < 0, arr.ind=TRUE)
    for(i in 1:nrow(w))
      Total_score[w[i,1],w[i,2]] <- 0
  }
  
  w = which(is.na(Total_score), arr.ind=TRUE)
  if(nrow(w) != 0){
    for(i in 1:nrow(w))
      Total_score[w[i,1],w[i,2]] <- 0
  }
  
  # output optimal algorithm names
  opt_algorithm <- Total_score %>% apply(., 1, which.max) %>% 
    names(Total_score)[.] %>% 
    setNames(., rownames(Total_score)) %>%
    .[-length(.)]
  
  Total_score.melt <- melt(cbind(x=rownames(Total_score), Total_score), id='x')
  
  result <- list(
                 Total_score           = Total_score,
                 Accuarcy              = Accuracy,
                 Precision             = Precision,
                 Effectiveness         = Asses_fz$Valid_percent,
                 Accuracy.list         = list(Accuracy_MAE = Accuracy_MAE,
                                              Accuracy_CMAPE = Accuracy_CMAPE,
                                              Accuracy_BIAS = Accuracy_BIAS,
                                              Accuracy_CMRPE = Accuracy_CMRPE),
                 Precision.list        = list(Precision_MAE = Precision_MAE,
                                              Precision_CMAPE = Precision_CMAPE,
                                              Precision_BIAS = Precision_BIAS,
                                              Precision_CMRPE = Precision_CMRPE),
                 Total_score.melt      = Total_score.melt,
                 opt_algorithm         = opt_algorithm,
                 Inputs                = Inputs
                )
  
  return(result)
  
}



#' @export
#' @rdname Scoring_system
#' @return The result of \code{Scoring_system_bootstrap} are including:
#' \itemize{
#'   \item \strong{Times} The times of bootstrap running.
#'   \item \strong{Score_all_clusters} The total score for algorithms across all clusters.
#'   \item \strong{Score_list} All times of bootstrapping results are recorded in it.
#'   \item \strong{Score_list_melt} Melted \code{Score_list}.
#'   \item \strong{Opt_algorithm_list} The optimal algorithm for every runing.
#'   \item \strong{Opt_algorithm} The optimal algorithm defined by mode of \code{Opt_algorithm_list}
#'     for each cluster.
#'   \item \strong{Remove_algorithm} The algorithms to be removed when blending.
#'   \item \strong{plot_col} The col plot of \code{Score_list_melt}.
#'   \item \strong{plot_scatter} The scatter plot of measured and predicted Chla concentration
#'     colored by clusters.
#'   \item \strong{plot_scatter_opt} The scatter plot of measured and predicted Chla concentration
#'     colored by clusters for optimized algorithms.
#'   \item \strong{Blend_result} The results from the inherent function \code{Chla_algorithms_blend}.
#'   \item \strong{dt_Chla} Data.frame with combination of candidate algortihms and blended results.
#'   \item \strong{Chla_blend} The blended Chla concentration by score results.
#'   \item \strong{Results_of_scoring_system} A list including all results of \link{Scoring_system} function.
#' } 
#' 
#' @importFrom stats aggregate
#' @importFrom stringr str_split
#' 
Scoring_system_bootstrap <- function(Times = 1000, 
                                     Inputs, replace = TRUE, 
                                     method = 'sort-based',
                                     param_sort = list(decreasing = TRUE),
                                     param_interval = list(trim=FALSE, reward.punishment=TRUE,
                                                           decreasing=TRUE, hundred.percent=FALSE),
                                     remove.negative = FALSE){
  
  if(Times <= 1) stop("Times should be larger than 1")
  
  pred = round(Inputs$Asses_fz$input$pred, 4)
  meas = round(Inputs$Asses_fz$input$meas, 4)
  memb = round(Inputs$Asses_fz$input$memb, 4)
  cluster = apply(memb, 1, which.max)
  K = ncol(memb)
  
  Score_list <- NULL
  Opt_algorithm_list <- NULL
  Results_of_scoring_system <- list()
  
  for(i in 1:Times){
    
    cat(i, "/", Times, "\n")
    Asses_list <- Getting_Asses_results(sample.size= nrow(memb),
                                        pred, meas, memb, replace = replace)
    Score <- Scoring_system(Inputs = Asses_list,
                           method = method,
                           param_sort = param_sort,
                           param_interval = param_interval,
                           remove.negative = remove.negative)
    
    Results_of_scoring_system[[i]] <- Score
    
    if(i == 1){
      Score_list <- Score$Total_score.melt
      names(Score_list)[ncol(Score_list)] <- i
    }else{
      Score_list <- cbind(Score_list, Score$Total_score.melt[,3])
      names(Score_list)[ncol(Score_list)] <- i
    }
    
    Opt_algorithm_list <- rbind(Opt_algorithm_list, Score$opt_algorithm)
    
  }
  
  # Opt algorithms for each time
  Opt_algorithm_list <- cbind(seq(1, Times), Opt_algorithm_list)
  Opt_algorithm_list <- data.frame(Opt_algorithm_list, stringsAsFactors = FALSE)
  colnames(Opt_algorithm_list) <- c("Times", names(Score$opt_algorithm))
  

  # Calculate the score statistic information
  u_arr   = apply(Score_list[, 3:ncol(Score_list)], 1, mean) %>% round(2)
  sig_arr = apply(Score_list[, 3:ncol(Score_list)], 1, sd) %>% round(2)
  ymin    = u_arr - sig_arr
  ymax    = u_arr + sig_arr
  
  res_Score <- cbind(Score_list[,1:2],
                     value    = u_arr,
                     sig_arr  = sig_arr) %>% level_to_variable
  Score_all_clusters = res_Score[res_Score$x == "SUM", -1]
  res_Score = res_Score[!(res_Score$x %in% 'SUM'),]
  res_Score$x_f = factor(res_Score$x, levels = paste("Cluster", (1:ncol(memb))))
  
  # The optimal algorithm defined by maximum scores per cluster
  tmpdt = stats::aggregate(res_Score$value, list(res_Score$x), which.max)
  tmpdt = tmpdt[tmpdt$Group.1 != "SUM",]
  Opt_algorithm = rep("", K) %>% setNames(., paste("Cluster", 1:K))
  for(i in 1:K){
    tmptmp = res_Score[which(res_Score$x == names(Opt_algorithm)[i]),]
    Opt_algorithm[i] <- tmptmp$variable[which.max(tmptmp$value)]
  }
  
  ## To-do: 
  ## In some extents, the optimal could not return finite value especially for image rasters,
  ##   mostly because of inappropriate atmospheric correction. Given that, it is necessary to select 
  ##   several (two or more) substitutes as a reinforcement.
  
  # define the position of optimal and removed algorithms
  res_Score$pos_opt <- res_Score$pos_removed <- res_Score$sig_arr * NA
  for(i in 1:length(Opt_algorithm)){
    w = which(res_Score$x == names(Opt_algorithm)[i] & res_Score$variable == Opt_algorithm[i])
    res_Score$pos_opt[w] = res_Score$value[w]
  }
  
  # The removed algorithms defined by the lower half res_Score per cluster
  remove_ratio = 1/2
  if(remove_ratio > 1) stop("The remove ratio should be between zero and one!")
  for(i in levels(res_Score$x_f)){
    tmpdt = res_Score[res_Score$x == i, ]
    tmpdt = tmpdt[tmpdt$variable %in% Opt_algorithm, ]
    num_removed = floor(length(unique(tmpdt$variable)) * remove_ratio)
    w_removed = which(res_Score$x == i & res_Score$value <= sort(tmpdt$value)[num_removed])
    res_Score$pos_removed[w_removed] = res_Score$value[w_removed]
  }
  ws_removed <- which(!is.na(res_Score$pos_removed))
  w_cancel <- which(!(res_Score$variable[ws_removed] %in% Opt_algorithm))
  res_Score$pos_removed[ws_removed][w_cancel] <- NA
 
  # For these algorithms, they are not taken into account when blending
  Remove_algorithm <- list()
  for(i in 1:length(unique(res_Score$x))){
    Remove_algorithm[[i]] <- res_Score$variable[which(!is.na(res_Score$pos_removed) & res_Score$x == unique(res_Score$x)[i])]
  }
  
  # add box
  res_Score$box <- res_Score$sig_arr * NA
  w_box <- which(res_Score$variable %in% Opt_algorithm)
  res_Score$box[w_box] <- res_Score$value[w_box]
  
  
  # plot Score_list
  num.model <- Score_list$variable %>% unique %>% length
  set.seed(1234)
  ind = runif(num.model) %>% sort.int(., index.return=TRUE) %>% .$ix
  
  dodge <- position_dodge(width=0.9)
  
  # plot_col <- ggplot(data=res_Score) +
  #   geom_col(aes(x=x_f, y=value, group=variable, fill=variable), position=dodge) + 
  #   geom_errorbar(aes(x=x_f, ymin=value-sig_arr, ymax=value+sig_arr, group=variable),
  #                 position=dodge, width=0.6) + 
  #   scale_fill_manual(values=Spectral(num.model)[ind]) +
  #   labs(y='Score', fill='Algorithm', x=NULL)+
  #   facet_wrap(~x_f, scales='free_x') + 
  #   theme_bw() + 
  #   theme(axis.text.x.bottom = element_blank(),
  #         axis.ticks.x.bottom = element_blank(),
  #         strip.background = element_rect(fill='white', color='white'),
  #         strip.text = element_text(face='bold'),
  #         text = element_text(size=13),
  #         legend.position = "right")
  plot_col <- 
    ggplot() +
    geom_col(data=res_Score, 
             aes(x=x_f, y=value, group=variable, fill=variable),
             alpha=0.75, position=dodge) + 
    geom_col(data=res_Score, 
             aes(x=x_f, y=box, group=variable, fill=variable), color="gray20",
             position=dodge) +
    geom_errorbar(data=res_Score, 
                  aes(x=x_f, ymin=value-sig_arr, ymax=value+sig_arr, group=variable),
                  position=dodge, width=0.6, alpha=1) + 
    geom_point(data=res_Score, 
               aes(x=x_f, y=pos_opt, group=variable),
               position=dodge, color="green", fill=NA, shape="diamond", stroke = 0.5, size=2) + 
    geom_point(data=res_Score, 
               aes(x=x_f, y=pos_removed, group=variable),
               position=dodge, color="red", fill=NA, shape="cross", stroke = 1) + 
    scale_fill_manual(values=Spectral(num.model)[ind]) +
    scale_y_continuous(limits=c(0, NA)) + 
    labs(y='Score', fill='Algorithm', x=NULL)+
    facet_wrap(~x_f, scales='free_x', nrow=3) + 
    theme_bw() + 
    theme(axis.text.x.bottom = element_blank(),
          axis.ticks.x.bottom = element_blank(),
          panel.grid.major.x = element_blank(),
          strip.background = element_rect(fill='white', color='white'),
          strip.text = element_text(face='bold'),
          panel.spacing.y = unit(0.5, "lines"),
          text = element_text(size=13),
          # legend.position = c(1,0),
          legend.position = "right",
          # legend.justification = c(1,0),
          legend.justification = c(0.5, 0.5),
          legend.text = element_text(size=8),
          legend.title = element_text(size=9, face="bold")) +
    # guides(fill = guide_legend(ncol=3))
    guides(fill = guide_legend(ncol=1))
  
  
  #-----------------------#
  #                       #
  #     Blending work     #
  #                       #
  #-----------------------#
  
  dt_Chla_opt <- pred[, Opt_algorithm]
  Blend_result <- Chla_algorithms_blend(dt_Chla_opt, memb, Opt_algorithm, Remove_algorithm)
  Chla_blend <- Blend_result$Chla_blend
  
  # draw a scatter plot colored by clusters
  Chla_limits <- (range(meas) * c(0.9, 1.1)) %>% round(4)
  if(Chla_limits[1] <= 0) Chla_limits[1] = 0.1
  Chla_plots <- data.frame(Chla_true=meas, cluster = cluster, pred, Chla_blend)
  Chla_plots$cluster <- factor(as.character(Chla_plots$cluster), levels = 1:K, ordered = TRUE)
  names(Chla_plots)[1] <- "Chla_true"
  dt_Chla_ <- reshape2::melt(Chla_plots, id=c("Chla_true", "cluster"))
  dt_Chla_$value_opt <- NA
  for(i in 1:K){
    w = which(dt_Chla_$cluster == i & dt_Chla_$variable == Opt_algorithm[i])
    dt_Chla_$value_opt[w] <- dt_Chla_$value[w]
  }
  w = which(dt_Chla_$variable == "Chla_blend")
  dt_Chla_$value_opt[w] <- dt_Chla_$value[w]
  
  cp = RdYlBu(K)
  names(cp) <- levels(dt_Chla_$cluster)
  plot_scatter <- 
    ggplot() + 
    geom_abline(intercept=0, slope=1, linetype=2) + 
    geom_point(data=dt_Chla_, 
               aes(x=Chla_true, y=value, color=cluster), 
               # shape = "circle open", stroke = 1, 
               alpha=0.3) + 
    geom_point(data=dt_Chla_, 
               aes(x=Chla_true, y=value_opt, fill=cluster), 
               shape = "circle filled", color = "black",
               alpha=1) + 
    geom_rug(data=dt_Chla_, 
             aes(x=Chla_true, y=value, color=cluster), alpha=0.2, show.legend = F) + 
    scale_x_log10(limits=Chla_limits) + scale_y_log10(limits=Chla_limits) +
    labs(x='Measured Chla [ug/L]', y='Predicted Chla [ug/L]', color='Cluster', fill="Cluster") + 
    scale_color_manual(values=cp) +
    scale_fill_manual(values=cp) + 
    theme_bw() + 
    facet_wrap(~variable) +
    # facet_grid(variable~cluster) + 
    theme(legend.position = "right",
          strip.background = element_rect(fill='white', color='white'),
          strip.text = element_text(face='bold'),
          text = element_text(size=13)) + 
    guides(col = guide_legend(ncol=1))
  
  dt_Chla_sub <- dt_Chla_[dt_Chla_$variable %in% c(
    "Chla_blend", unique(Opt_algorithm)
  ),]
  plot_scatter_opt <- 
    ggplot() + 
    geom_abline(intercept=0, slope=1, linetype=2) + 
    geom_point(data=dt_Chla_sub, 
               aes(x=Chla_true, y=value, color=cluster), 
               # shape = "circle open", stroke = 1, 
               alpha=0.3) + 
    geom_point(data=dt_Chla_sub, 
               aes(x=Chla_true, y=value_opt, fill=cluster), 
               shape = "circle filled", color = "black",
               alpha=1) + 
    geom_rug(data=dt_Chla_sub, 
             aes(x=Chla_true, y=value, color=cluster), alpha=0.2, show.legend = F) + 
    scale_x_log10(limits=Chla_limits) + scale_y_log10(limits=Chla_limits) +
    labs(x='Measured Chla [ug/L]', y='Predicted Chla [ug/L]', color='Cluster', fill="Cluster") + 
    scale_color_manual(values=cp) +
    scale_fill_manual(values=cp) + 
    theme_bw() + 
    # facet_wrap(~variable) +
    facet_grid(variable~cluster) +
    theme(legend.position = "right",
          strip.background = element_rect(fill='white', color='white'),
          strip.text = element_text(face='bold'),
          text = element_text(size=13)) + 
    guides(col = guide_legend(ncol=1))
  
  # return outputs
  result = list(
    Times              = Times,
    Score_all_clusters = Score_all_clusters,
    Score_list         = Score_list, 
    Score_list_melt    = res_Score,
    Opt_algorithm_list = Opt_algorithm_list,
    Opt_algorithm      = Opt_algorithm,
    Remove_algorithm   = Remove_algorithm,
    plot_col           = plot_col,
    plot_scatter       = plot_scatter,
    plot_scatter_opt   = plot_scatter_opt,
    Blend_result       = Blend_result, 
    dt_Chla            = Chla_plots,
    Chla_blend         = Chla_blend,
    Results_of_scoring_system = Results_of_scoring_system
  )
  
  return(result)
  
}
  

#' @name Chla_algorithms_blend
#' 
#' @title Algorithms blending via membership values based on optimal candidates
#' 
#' @param dt_Chla_opt The data.frame includes estimation by optimal algorithms per cluster
#' @param memb Membership values
#' @param Opt_algorithm Optimal algorithm name
#' @param Remove_algorithm Removed algorithm name
#' 
#' @noRd
#' @return The result of \code{Chla_algorithm_blend} is a list including:
#'   \itemize{
#'     \item \strong{modified_memb_list} The vector recording the modified membership values across all samples.
#'     \item \strong{Chla_blend} The final blended Chla result after modification.
#'     \item \strong{new_memb} The new membership values after modification.
#'   }
#' 
#' @importFrom stringr str_split
#' @importFrom magrittr %>%
Chla_algorithms_blend <- function(dt_Chla_opt, memb, Opt_algorithm, Remove_algorithm){
  
  if(!all.equal(
    length(Opt_algorithm),
    ncol(dt_Chla_opt),
    ncol(memb)
  )){
    stop("The input dt_Chla_opt, memb, and Opt_algorithm should have same length or column number!")
  }
  
  if(nrow(dt_Chla_opt) != nrow(memb)) 
    stop("The row number of dt_Chla_opt and memb should be equal!")

  # if(all(Opt_algorithm %in% names(dt_Chla_opt))){
  #   stop("The Opt_algorithm should be parts of dt_Chla_opt colnames")
  # }
  
  if(!(is.list(Remove_algorithm) | is.null(Remove_algorithm))){
    stop("Remove_algorithm should be a list or NULL.")
  }
  
  cluster <- apply(memb, 1, which.max)
  
  dt_Chla_opt <- round(dt_Chla_opt, 4)
  memb        <- round(memb, 4)
  
  dt_Chla_opt_ <- dt_Chla_opt
  memb_     <- memb
  modified_memb_list <- rep(NA, nrow(memb_))
  
  # modify the membership values
  for(i in 1:nrow(memb_)){
    # removed algorithms
    if(!is.null(Remove_algorithm)){
      w_rm <- which(str_split(names(dt_Chla_opt_), "\\.", simplify = TRUE)[, 1] %in%
                      Remove_algorithm[[cluster[i]]])
      modified_memb <- sum(memb_[i, w_rm])
    }else{
      w_rm = NULL
      modified_memb <- 0
    }
    modified_memb_list[i] <- modified_memb
    memb_[i, w_rm] <- 0
    cluster_pos <- as.numeric(cluster[i])
    memb_[i, cluster_pos] <- memb_[i, cluster_pos] + modified_memb
    
    # failure avoid
    w_fail <- which(is.na(dt_Chla_opt[i,] * memb_[i, ]) | (dt_Chla_opt[i,] * memb_[i, ]) <= 0)
    if(length(w_fail) == 0){
      modified_memb <- 0
    }else{
      modified_memb <- sum(memb_[i, w_fail])
    }
    modified_memb_list[i] <- modified_memb_list[i] + modified_memb
    memb_[i, w_fail] <- 0
    memb_[i, cluster_pos] <- memb_[i, cluster_pos] + modified_memb
    
    # compensate the sum of membership to one
    memb_bias = 1 - sum(memb_[i, ])
    memb_[i, cluster_pos] = memb_[i, cluster_pos] + memb_bias
  }
  
  Chla_blend <- apply(dt_Chla_opt_ * memb_, 1, function(x) sum(x, na.rm=TRUE)) 
  
  result <- list(
    modified_memb_list = modified_memb_list,
    Chla_blend         = Chla_blend,
    new_memb           = memb_
  )
  
}  
  
