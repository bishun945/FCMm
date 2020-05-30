#' @name Assessment_via_cluster
#' @title Assessment each algorithm for every cluster
#' @param pred prediciton of Chla
#' @param meas in-situ measurement of Chla
#' @param memb membership value matrix
#' @param metrics metrics need to be calculated
#' @param log10 Should pred and meas be log10-transformed? (Default as FALSE)
#' @param total logical
#' @param hard.mode If \code{FALSE}, the membership values are used to calculate validation metrics
#' @param cal.precision Whether to calculate precision (only support for vectorized metrics), default as \code{FALSE}
#' @param na.process na.process and choose to statistic na value percent
#' @param plot.col option to plot col result for chosed metrics (Default as FALSE)
#' 
#' @note If the cal.precision is running, the \code{hard.mode == TRUE} is used. In that case,
#'   mean and sd calculation is conducted for hard mode based on result from cal.metrics.vector
#' @export
#' @return List
#' @family Algorithm assessment
#' 
#' @import ggplot2  
#' @importFrom stats runif
#' @importFrom reshape2 melt
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
      stop('Not choose to process NA values. Predicted, measured or membership values including NA values!')
    }
  }else{
    if(anyNA(meas)){
      stop('Choose to process NA values. Measured values including NA values!')
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
        tmp <- cbind(x,y) %>% na.omit
        x <- tmp[,1]
        y <- tmp[,2]
        num_new <- length(x)
        if(num_new > num_raw)
          stop("Error! The subseted sample number is smaller the raw.")
      }
      
      # calculate precision
      if(cal.precision == TRUE){
      
        for(metric in metrics){
          metric_value <- cal.metrics.vector(x,y,metric,log10=log10)
          result[[metric]][cluster, model] <- mean(metric_value, na.rm=T)
          precision_name <- paste0(metric, '_p')
          result[[precision_name]][cluster, model] <- trim_sd(metric_value, na.rm=T)
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
        tmp <- cbind(x,y) %>% na.omit
        x <- tmp[,1]
        y <- tmp[,2]
        num_new <- length(x)
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
            result[[metric]]['SUM', model] <- mean(metric_value, na.rm=T)
            precision_name <- paste0(metric, '_p')
            result[[precision_name]]['SUM', model] <- sd(metric_value, na.rm=T)
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
          
          if(na.process){
            w <- which(is.na(y) == FALSE)
            x <- x[w]
            y <- y[w]
            memb_ <- memb[w,]
          }
          
          if(dim(memb_)[1] != length(x))
            stop("The fuzzy metrics are calculated with different rows from memb and pred")
          
          Er <- cal.metrics.vector(x,y,metric,log10)
          result_fz[[metric]][i,j] <- sum(memb_[,i] * Er, na.rm=T)/sum(memb_[,i], na.rm=T)
          
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
    ind = runif(num.model) %>% sort.int(., index.return=T) %>% .$ix
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
#' @title Score a vector of error metrics for algorithms
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
#' @return List
#' @family Algorithm assessment
#' 
#' @importFrom stats qt sd
#' @import ggplot2

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
    u <- mean(x_, na.rm=T)
    sig <- sd(x_, na.rm=T)
  }else{
    u <- mean(x, na.rm=T)
    sig <- sd(x, na.rm=T)
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
#' @param num The number of sampling
#' @param replace The way of sampling. See \code{help(sample)} for more details
#' @export
#' @return List
#' @family Algorithm assessment
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
    
    types_num_final <- types_num
    rand_select <- sample(types, resi, replace=TRUE)
    for(i in rand_select){
      types_num_final[i] <- types_num_final[i] + 1
    }
    
    # If one cluster having sample number larger than the total OR not getting the target number then re-run
    while(!all(types_num_final <= raw_num) | sum(types_num_final) < num){
      types_num_final <- types_num
      rand_select <- sample(types, resi, replace=TRUE)
      for(i in rand_select){
        types_num_final[i] <- types_num[i] + 1
      }
    }
    
  }
  
  # Start to sample
  ind <- seq(length(x))
  
  result <- NULL
  
  for(i in 1:length(types_num_final)){
    
    x_type <- ind[which(x == i)]
    tmp_sample <- sample(x_type, types_num_final[i], replace=replace)
    result <- c(result, tmp_sample)
    
  }
  
  return(result)
  
}


#' @name Getting_Asses_results
#' @title Get the result of function Assessment_via_cluster
#' @description This function mainly use function \code{Assessment_via_cluster} to get
#'   assesments both from fuzzy and hard mode. Specifically, it will return fuzzy_MAE, 
#'   fuzzy_MAPE, hard_Slope, Valid_percent, and hard_Rsquare which should be seemed as the 
#'   input value of function \code{Scoring_system}.
#' @param sample.size Sample size. This supports a bootstrap way to run the function 
#'   \code{Assessment_via_cluster}. The number should not be larger than the row number 
#'   of pred or so.
#' @param replace Logical, replace, default as \code{FALSE}
#' @param pred Prediction matrix
#' @param meas Measured (actual) vector
#' @param memb Membership matrix
#' @param cluster Cluster vector
#' @param seed Seed number for fixing the random process. See \code{help(set.seed)} for more details.
#' @note The row number of \code{pred}, \code{meas}, \code{memb}, and \code{cluster} should be the same. 
#'   This function is designed for bootstrapping process to get Chla algorithms assessment. Therefore, 
#'   parameters of \code{Assessment_via_cluster} is set as fixed such as \code{log10 = TRUE}, 
#'   \code{na.process = TRUE}. Given that, I will not export this function in latter to avoid confuses.
#' @export
#' @return A list containing fuzzy and hard results from \code{Assessment_via_cluster}
#' @family Algorithm assessment
#' 
Getting_Asses_results <- function(sample.size, replace=FALSE,
                                  pred, meas, memb, cluster, 
                                  seed=NULL){
  
  if(sample.size > nrow(pred))
    stop("Enter a smaller sample size to subset the data set.")
  
  if(var(c(nrow(pred), nrow(meas), nrow(memb), nrow(cluster))) != 0)
    stop("Row number of pred, meas, memb and cluster is different.")
  
  # Stratified sampling by cluster
  if(is.null(seed)){
    w <- Sampling_via_cluster(cluster, num=sample.size, replace=replace)
  }else{
    set.seed(as.numeric(seed))
    w <- Sampling_via_cluster(cluster, num=sample.size, replace=replace)
  }
  
  pred_ <- pred[w,]
  meas_ <- meas[w,]
  memb_ <- memb[w,]
  cluster_ <- cluster[w]
  
  Asses_fz <-  Assessment_via_cluster(pred=pred_,
                                      meas=meas_,
                                      memb=memb_,
                                      metrics = c("MAE","SMAPE","BIAS",'SMRPE'),
                                      log10 = T,
                                      hard.mode = F,
                                      na.process = TRUE,
                                      plot.col = F)
  
  Asses_p <-  Assessment_via_cluster(pred=pred_,
                                      meas=meas_,
                                      memb=memb_,
                                      metrics = c("MAE","SMAPE","BIAS",'SMRPE'),
                                      log10 = T,
                                      hard.mode = F,
                                      cal.precision = T,
                                      na.process = TRUE,
                                      plot.col = F)
  
  # NOTE 2020-03-06: As the difination of Accuracy and Precision were changed, 
  #   we deleted the following codes and corresponding contexts in `Scoring_system` function
  # 
  # Asses_hd <-  Assessment_via_cluster(pred=pred_,
  #                                     meas=meas_,
  #                                     memb=memb_,
  #                                     metrics = c("SLOPE", "R2_SMA"),
  #                                     log10 = T,
  #                                     hard.mode = T,
  #                                     na.process = TRUE,
  #                                     plot.col = F)
  
  result <- list(
                  # Asses_hd=Asses_hd, 
                  Asses_fz = Asses_fz,
                  Asses_p  = Asses_p 
                  )
  
  return(result)
}


#' @name Scoring_system
#' @title The main function for algorithms scoring
#' @param Inputs The list returned form function \code{Getting_Asses_results}
#' @param method The method selected to score algorithms: 'sort-based' (default) or 'interval-based'
#' @param param_sort The parameters of function \code{Score_algorithms_sort}
#' @param param_interval The parameters of function \code{Score_algorithms_interval}
#' @param remove.negative Option to replace the negative score as zero (Default as \code{FALSE})
#' @export
#' @return A list including Total_score, Accuracy, Precision, Effectiveness, Accuracy.list, 
#'   Precision.list, and Total_score.melt.
#' @family Algorithm assessment
#' 
#' @importFrom reshape2 melt
Scoring_system <- function(Inputs, 
                           method = 'sort-based',
                           param_sort = list(decreasing = TRUE),
                           param_interval = list(trim=FALSE, reward.punishment=TRUE,
                                             decreasing=TRUE, hundred.percent=FALSE),
                           remove.negative=FALSE){
  
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
  # The accuracy and precision is newly defined in this package (referred by Hooker):
  # Accuracy is the estimation of how close the result of the experiment is to the 
  #   true value.
  # Precision is the estimation of how excatly the result is determined independently 
  #   of any true value.
  # In other words, accuracy is telling a story truthfully and precision is how similarly
  #   the story is represented over and over again.
  # 
  # Here we use AE, a vector for each sample, for instance.
  # Accuracy is the aggregation (no matter mean or median, in fuzzy calculation process), 
  #   we use mean to some extent. 
  # While, precision is acutally the stability of AE (reproducebility) which means the error
  #   produced by the algorithm is under certain control.
  #
  Accuracy_SMRPE <- Accuracy_BIAS <- Accuracy_MAE <- Accuracy_SMAPE <- Asses_fz$MAE * NA
  for(i in 1:nrow(Accuracy_MAE)){
    # Accuracy_MAE[i,]  <- Score_algorithms(Asses_fz$MAE[i,], trim=trim)$score
    # Note 2020-03-06:
    # Bias was banned as it sometimes has the issue that may arise at the aggregation the phase, 
    #   when the positive and negative errors will be cancelling each other. The use of this relatve
    #   (not absolute) often demonstrate a falsely high accuracy.
    # Note 2020-03-07:
    # I think Bias is okay now.
    # Accuracy_BIAS[i,]  <- Score_algorithms(abs(Asses_fz$BIAS[i,]), trim=trim)$score
    # Accuracy_SMAPE[i,] <- Score_algorithms(Asses_fz$SMAPE[i,], trim=trim)$score
    # Accuracy_SMRPE[i,] <- Score_algorithms(abs(Asses_fz$SMRPE[i,]), trim=trim)$score
    
    # Test in 2020-03-08 using sorted score
    Accuracy_MAE[i,]   <- Score_algorithms(Asses_fz$MAE[i,])
    Accuracy_BIAS[i,]  <- Score_algorithms(abs(Asses_fz$BIAS[i,]))
    Accuracy_SMAPE[i,] <- Score_algorithms(Asses_fz$SMAPE[i,])
    Accuracy_SMRPE[i,] <- Score_algorithms(abs(Asses_fz$SMRPE[i,]))
  }
  Accuracy <- Accuracy_MAE + Accuracy_BIAS + Accuracy_SMRPE + Accuracy_SMAPE
  
  
  # Precision part
  Precision_SMRPE <- Precision_BIAS <- Precision_MAE <- Precision_SMAPE <- Asses_p$MAE_p * NA
  for(i in 1:nrow(Precision_MAE)){
    # Precision_MAE[i,]    <- Score_algorithms(Asses_p$MAE_p[i,], trim=trim)$score
    # Precision_BIAS[i,]    <- Score_algorithms(Asses_p$BIAS_p[i,], trim=trim)$score
    # Precision_SMAPE[i,]  <- Score_algorithms(Asses_p$SMAPE_p[i,], trim=trim)$score
    # Precision_SMRPE[i,]  <- Score_algorithms(Asses_p$SMRPE_p[i,], trim=trim)$score
    
    # Test in 2020-03-08 using sorted score
    Precision_MAE[i,]   <- Score_algorithms(Asses_p$MAE_p[i,])
    Precision_BIAS[i,]  <- Score_algorithms(abs(Asses_p$BIAS_p[i,]))
    Precision_SMAPE[i,] <- Score_algorithms(Asses_p$SMAPE_p[i,])
    Precision_SMRPE[i,] <- Score_algorithms(abs(Asses_p$SMRPE_p[i,]))
  }
  Precision <- Precision_MAE + Precision_BIAS + Precision_SMRPE + Precision_SMAPE

  Total_score <- (Accuracy + Precision) * Asses_fz$Valid_percent / 100
  
  if(remove.negative == TRUE){
    w = which(Total_score < 0, arr.ind=T)
    for(i in 1:nrow(w))
      Total_score[w[i,1],w[i,2]] <- 0
  }
  
  Total_score.melt <- melt(cbind(x=rownames(Total_score), Total_score), id='x')
  
  result <- list(
                 Total_score           = Total_score,
                 Accuarcy              = Accuracy,
                 Precision             = Precision,
                 Accuracy.list         = list(Accuracy_MAE = Accuracy_MAE,
                                              Accuracy_SMAPE = Accuracy_SMAPE,
                                              Accuracy_BIAS = Accuracy_BIAS,
                                              Accuracy_SMRPE = Accuracy_SMRPE),
                 Precision.list        = list(Precision_MAE = Precision_MAE,
                                              Precision_SMAPE = Precision_SMAPE,
                                              Precision_BIAS = Precision_BIAS,
                                              Precision_SMRPE = Precision_SMRPE),
                 Total_score.melt      = Total_score.melt,
                 Inputs                = Inputs
                )
  
  return(result)
  
}

