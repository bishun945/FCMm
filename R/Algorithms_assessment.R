#' @title Assessment each algorithm for every cluster
#' @name Assessment_via_cluster
#' @param pred prediciton of Chla
#' @param meas in-situ measurement of Chla
#' @param memb membership value matrix
#' @param metrics metrics need to be calculated
#' @param log10 Should pred and meas be log10-transformed? (Default as FALSE)
#' @param total logical
#' @param hard.mode If \code{FALSE}, the membership values are used to calculate validation metrics
#' @param na.process na.process and choose to statistic na value percent
#' @param plot.col option to plot col result for chosed metrics (Default as FALSE)
#' @export
#' @return List
#' @family Algorithm assessment
Assessment_via_cluster <- function(pred, meas, memb,
                                   metrics = c('MAE','MAPE'),
                                   log10 = FALSE, 
                                   total = TRUE,
                                   hard.mode= TRUE, 
                                   na.process = FALSE,
                                   plot.col = FALSE
){
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
  
  
  # generate the output dataframe  
  model_names <- colnames(pred)
  cluster_names <- colnames(memb)
  cluster_crisp <- apply(memb,1,which.max)
  
  validation <- matrix(data=0,
                       nrow=length(cluster_names), 
                       ncol=length(model_names)) %>% as.data.frame()
  colnames(validation) <- model_names
  rownames(validation) <- cluster_names
  
  # output is a list
  result <- list()
  for(i in 1:length(metrics))
    result[[metrics[i]]] <- validation
  
  if(na.process){
    result[["NA_percent"]] <- validation
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
      
      for(metric in metrics){
        result[[metric]][cluster, model] <- cal.metrics(x,y,metric,log10=log10)
      }
      if(na.process){
        result[["NA_percent"]][cluster, model] <- num_new / num_raw * 100
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
        result[[metric]]['SUM', model] <- cal.metrics(x,y,metric,log10=log10)  
      }
      
      if(na.process){
        result[["NA_percent"]]['SUM', model] <- num_new / num_raw * 100
      }
    }
  }
  
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
  
  if(plot.col == TRUE){
    
    res_plot <- list()
    res_plot_dt <- list()
    
    for(metric in names(result)){
      
      tmp <- result[[metric]]
      
      tmp <- cbind(x=rownames(tmp), tmp)
      
      tmp <- melt(tmp, id="x", variable.name='Models')
      
      p <- ggplot(tmp) + 
        geom_col(aes(x=x,y=value,group=Models,fill=Models),
                 position="dodge") + 
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
      labs(x=NULL, y=NULL) + 
      theme_bw() + 
      facet_wrap(~Metric, scales="free_y")
    
  }else{
    
    result$res_plot <- NULL
    result$res_plot_dt <- NULL
    result$res_plot_facet <- NULL
    
  }
  
  result$input <- list(
    pred=pred,
    meas=meas,
    memb=memb,
    metrics=metrics,
    log10=log10,
    total=total,
    hard.mode=hard.mode,
    na.process=na.process
  )
  
  return(result)
  
}

#' @title Score a vector of error metrics for algorithms
#' @name Score_algorithms
#' @param x Input vector of error metrics (NA values are allowed)
#' @param trim whether to run a trim process to calculate mean and standard deviation of 
#'   input vector x (Default as \code{FALSE})
#' @param reward.punishment wheter to conduct the reward and punishment mechainism in scoring
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

Score_algorithms <- function(x, trim=FALSE, reward.punishment=TRUE, 
                             decreasing=TRUE, hundred.percent=FALSE){
  
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
  
  p <- ggplot(tmp) + geom_point(aes(x,y)) + 
    geom_hline(yintercept=u,color='blue') + 
    geom_hline(yintercept=bds,color='blue',linetype=2) + 
    scale_x_continuous(breaks=seq(length(x))) +
    labs(x='Model ID', y='Error metric') + 
    theme_bw()
  
  score <- x * 0
  
  if(hundred.percent == TRUE){
    if(max(x) > 100 | min(x) < 0){
      stop("The input x has a number larger than 100 in NA_percent metric!")
    }
    score[which(x >= 95)]          = 0
    score[which(x >= 75 & x < 95)] = -0.5
    score[which(x >= 60 & x < 75)] = -1    
    score[which(x >= 50 & x < 60)] = -1.5
    score[which(x <  50)]          = -2
    
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
                 score = score)
  
  return(result)
}


#' @title Stratified sampling by clusters or types
#' @name Sampling_via_cluster
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
  
  types <- unique(x)
  types_num <- as.matrix(table(x)) / length(x) * num
  types_num <- floor(types_num)
  resi <- num - sum(types_num)
  
  if(resi < 0){
    stop('Stratifed sampling total number is greater than the sample number.')
  }else if(resi == 0){
    types_num <- types_num
  }else if(resi > 0){
    
    rand_select <- sample(types, resi, replace=TRUE)
    for(i in rand_select){
      types_num[i] <- types_num[i] + 1
    }
  }
  
  ind <- seq(length(x))
  
  result <- NULL
  
  for(i in 1:length(types_num)){
    
    x_type <- ind[which(x == i)]
    tmp_sample <- sample(x_type, types_num[i], replace=replace)
    result <- c(result, tmp_sample)
    
  }
  
  return(result)
  
}


#' @title Get the result of function Assessment_via_cluster
#' @description This function mainly use function \code{Assessment_via_cluster} to get
#'   assesments both from fuzzy and hard mode. Specifically, it will return fuzzy_MAE, 
#'   fuzzy_MAPE, hard_Slope, NA_percent, and hard_Rsquare which should be seemed as the 
#'   input value of function \code{Scoring_system}.
#' @name Getting_Asses_result
#' @param sample.size Sample size. This supports a bootstrap way to run the function 
#'   \code{Assessment_via_cluster}. The number should not be larger than the row number 
#'   of pred or so.
#' @param pred Prediction matrix
#' @param meas Measured (actual) vector
#' @param memb Membership matrix
#' @param cluster Cluster vector
#' @param seed Seed number for fixing the random process. See \code{help(set.seed)} for more details.
#' @note The row number of \code{pred}, \code{meas}, \code{memb}, and \code{cluster} should be the same. 
#'   This function is designed for bootstraping process to get Chla algorithms assessment. Therefore, 
#'   parameters of \code{Assessment_via_cluster} is set as fixed such as \code{log10 = TRUE}, 
#'   \code{na.process = TRUE}. Given that, I will not export this function in latter to avoid confuses.
#' @export
#' @return A list containing fuzzy and hard results from \code{Assessment_via_cluster}
#' @family Algorithm assessment
#' 
Getting_Asses_results <- function(sample.size, pred, meas, memb, cluster, seed=NULL){
  
  if(sample.size > nrow(pred))
    stop("Enter a smaller sample size to subset the data set.")
  
  if(var(c(nrow(pred), nrow(meas), nrow(memb), nrow(cluster))) != 0)
    stop("Row number of pred, meas, memb and cluster is different.")
  
  # Stratified sampling by cluster
  if(is.null(seed)){
    w <- Sampling_via_cluster(cluster, num=sample.size)
  }else{
    set.seed(as.numeric(seed))
    w <- Sampling_via_cluster(cluster, num=sample.size)
  }
  
  pred_ <- pred[w,]
  meas_ <- meas[w,]
  memb_ <- memb[w,]
  cluster_ <- cluster[w]
  
  Asses_fz <-  Assessment_via_cluster(pred=pred_,
                                      meas=meas_,
                                      memb=memb_,
                                      metrics = c("MAE","MAPE"),
                                      log10 = T,
                                      hard.mode = F,
                                      na.process = TRUE,
                                      plot.col = F)
  
  Asses_hd <-  Assessment_via_cluster(pred=pred_,
                                      meas=meas_,
                                      memb=memb_,
                                      metrics = c("SLOPE", "R2_SMA"),
                                      log10 = T,
                                      hard.mode = T,
                                      na.process = TRUE,
                                      plot.col = F)
  
  result <- list(Asses_fz=Asses_fz,
                 Asses_hd=Asses_hd)
  
  return(result)
}


#' @title The main function for algorithms scoring
#' @name Scoring_system
#' @param Inputs The list returned form function \code{Getting_Asses_results}
#' @param trim The input parameter of function \code{Score_algorithms} (Default as \code{FALSE})
#' @param remove.negative Option to repalce the negative score as zero (Default as \code{FALSE})
#' @export
#' @return A list including Total_score, Accuracy, Precision, Effectiveness, Accuracy.list, 
#'   Precision.list, and Total_score.melt.
#' @family Algorithm assessment
#' 

Scoring_system <- function(Inputs, trim=FALSE, remove.negative=FALSE){
  
  Asses_fz <- Inputs$Asses_fz
  Asses_hd <- Inputs$Asses_hd
  
  Accuracy_R2 <- Accuracy_Ratio <- Asses_hd$R2_SMA * NA
  for(i in 1:nrow(Accuracy_R2)){
    Accuracy_R2[i,]    <- Score_algorithms(Asses_hd$R2_SMA[i,], decreasing=F, trim=trim)$score
    Accuracy_Ratio[i,] <- Score_algorithms(abs(Asses_hd$SLOPE-1)[i,], trim=trim)$score
  }
  Accuracy <- Accuracy_R2 + Accuracy_Ratio
  
  Precision_MAE <- Precision_MAPE <- Asses_fz$MAE * NA
  for(i in 1:nrow(Precision_MAE)){
    Precision_MAE[i,]  <- Score_algorithms(Asses_fz$MAE[i,], trim=trim)$score
    Precision_MAPE[i,]  <- Score_algorithms(Asses_fz$MAPE[i,], trim=trim)$score
  }
  Precision <- Precision_MAE + Precision_MAPE
  
  Effectiveness <- Asses_fz$NA_percent * NA
  for(i in 1:nrow(Effectiveness)){
    Effectiveness[i,] <- Score_algorithms(Asses_fz$NA_percent[i,],
                                          trim=trim, 
                                          hundred.percent = TRUE)$score
  }
  
  Total_score <- Accuracy + Precision + Effectiveness
  
  if(remove.negative == TRUE){
    w = which(Total_score < 0, arr.ind=T)
    for(i in 1:nrow(w))
      Total_score[w[i,1],w[i,2]] <- 0
  }
  
  Total_score.melt <- melt(cbind(x=rownames(Total_score), Total_score), id='x')
  
  result <- list(Total_score      = Total_score,
                 Accuarcy         = Accuracy,
                 Precision        = Precision,
                 Effectiveness    = Effectiveness,
                 Accuracy.list    = list(Accuracy_R2, Accuracy_Ratio),
                 Precision.list   = list(Precision_MAE, Precision_MAPE),
                 Total_score.melt = Total_score.melt
  )
  
  return(result)
  
}

