#' @title Estimate Chla concentration by algorithms blending via membership values
#' @name FCM_m_Chla_estimation
#' @description
#' To calculate the Chla concentration via blending remote sensing algorithms.
#'
#' @usage FCM_m_Chla_estimation(Rrs, U)
#'
#' @param Rrs Data.frame of Remote sensing reflectance which should
#'   contain colname with (also at) \code{Rrs665}, \code{Rrs709},
#'   and \code{Rrs754} with unit sr^-1.
#' @param U Data.frame of membership values which should have seven
#'   columns presenting seven membership values from each cluster produced by FCM-m
#'
#' @return Each column presents Chla concentration with unit ug/L or mg/m^3 by model
#'   BR, TBA, C6 and blending result. Also with the membership values of each cluster.
#'
#' @note The input of \code{Rrs} must have bands with wavelength at 665, 709 and 754 nm.
#'   See examples of using this function.
#'
#' @export
#' 
#' @examples
#' \dontrun{
#' library(FCMm)
#' data("WaterSpec35")
#' data("Bi_clusters")
#' Rrs <- WaterSpec35[,3:17]
#' result <- apply_FCM_m(Rrs=Rrs, option.plot=TRUE)
#' dt_Chla <- FCM_m_Chla_estimation(Rrs=data.frame(Rrs665=Rrs$`665`, 
#'   Rrs709=Rrs$`708.75`,Rrs754=Rrs$`753.75`),U=result$u)
#' }
#' 
#' @references
#' \itemize{
#'   \item Bi S, Li Y, Xu J, et al. Optical classification of inland waters based on
#'     an improved Fuzzy C-Means method[J]. Optics Express, 2019, 27(24): 34838-34856.
#'   \item Gilerson A A, Gitelson A A, Zhou J, et al. Algorithms for remote estimation 
#'     of chlorophyll-a in coastal and inland waters using red and near infrared bands[J].
#'     Optics Express, 2010, 18(23): 24109-24125.
#' }
#' 
#' 

FCM_m_Chla_estimation <- function(Rrs, U){
  if(missing(Rrs)|missing(U))
    stop('Missing input variables: Rrs or U')
  if(!is.matrix(Rrs)&!is.data.frame(Rrs))
    stop('Input Rrs must be a matrix or data.frame!')
  if(!is.matrix(U)&!is.data.frame(U))
    stop('Input Membership U must be a matrix or data.frame!')
  if(nrow(Rrs)!=nrow(U))
    stop('The line of Rrs and U must match!')
  if(ncol(Rrs)!=3)
    stop('Only three bands needed in Rrs!')
  n <- c("Rrs665","Rrs709","Rrs754")
  for(i in names(Rrs)){
    if(sum(n %in% i) != 1)
      stop(paste0('Cannot find band ',i))
  }
  k <- ncol(U)
  if(k != 7)
    stop('The default cluster number should be set as seven.')
  if(anyNA(Rrs)|anyNA(U))
    stop('Please clean the NA vlaues in Rrs or U.')

  message("Note: this function is designed for OLCI or MERIS band settings!")
  message("The following bands (also as their names) must be contained in Rrs:")
  message("Rrs665 Rrs709 Rrs754")

  Rrs665 <- Rrs[, n[1]]
  Rrs709 <- Rrs[, n[2]]
  Rrs754 <- Rrs[, n[3]]

  names(U) <- seq(1,k) %>% as.character
  bind.Chla <- data.frame(BR=BR_Gil10(Rrs709,Rrs665),
                          TBA=TBA_Gil10(Rrs665,Rrs709,Rrs754),
                          C6=C6(Rrs665,Rrs754),
                          M=round(U,4))
  # TBA: C1 C2 C5
  # BR: C3 C4 C7
  # C6: C6
  U3 <- data.frame(u.TBA=bind.Chla[,c('M.1','M.2','M.5')] %>% apply(.,1,sum),
                   u.BR =bind.Chla[,c('M.3','M.4','M.7')] %>% apply(.,1,sum),
                   u.BR =bind.Chla[,c('M.6')])
  CONC3 <- bind.Chla[,c('TBA','BR','C6')]
  bind.Chla$conc.Blend <- apply(U3*CONC3, 1, sum)

  return(bind.Chla)
}

#' @title BR_Gil10
#' @param Rrs709 Rrs709
#' @param Rrs665 Rrs665
#' @export
BR_Gil10 <- function(Rrs709, Rrs665){
  return(abs(35.75*Rrs709/Rrs665-19.3)^1.124)
}

#' @title TBA_Gil10
#' @param Rrs665 Rrs665
#' @param Rrs709 Rrs709
#' @param Rrs754 Rrs754
#' @export
TBA_Gil10 <- function(Rrs665, Rrs709, Rrs754){
  return(abs(113.36*(1/Rrs665-1/Rrs709)*Rrs754+16.45)^1.124)
}

#' @title C6
#' @param Rrs665 Rrs665
#' @param Rrs754 Rrs754
#' @export
C6 <- function(Rrs665, Rrs754){
  return(10^( 1/Rrs665*Rrs754 * 0.14 + 2.11))
}

#' @title OC4E
#' @param Rrs443 Rrs443
#' @param Rrs490 Rrs490
#' @param Rrs510 Rrs510
#' @param Rrs555 Rrs555
#' @export
OC4E <- function(Rrs443, Rrs490, Rrs510, Rrs555){
  X <- apply(cbind(Rrs443,Rrs490,Rrs510),1,max)/Rrs555 %>% log10
  return(10^(0.3255-2.7677*X+2.4409*X^2-1.1288*X^3-0.4990*X^4))
}

#' @title BR_Git11
#' @param Rrs709 Rrs709
#' @param Rrs665 Rrs665
#' @export
BR_Git11 <- function(Rrs709, Rrs665){
  return(72.66*Rrs709/Rrs665-46.535)
}

#' @title TBA_Git11
#' @param Rrs665 Rrs665
#' @param Rrs709 Rrs709
#' @param Rrs754 Rrs754
#' @export
TBA_Git11 <- function(Rrs665, Rrs709, Rrs754){
  return(243.86*(1/Rrs665-1/Rrs709)*Rrs754+23.17)
}

#' @title NDCI_Mi12
#' @param Rrs665 Rrs665
#' @param Rrs708 Rrs708
#' @export
NDCI_Mi12 <- function(Rrs665,Rrs708){
  ind <- (Rrs708-Rrs665)/(Rrs708+Rrs665)
  return(14.039+86.115*ind+194.325*ind^2)
}

#' @title FBA_Le13
#' @param Rrs665 Rrs665
#' @param Rrs682 Rrs682
#' @param Rrs709 Rrs709
#' @export
FBA_Le13 <- function(Rrs665, Rrs682, Rrs709){
  return(18.492*(1/Rrs665-1/Rrs682)/(1/Rrs709-1/Rrs682)+6.1302)
}

#' @title FBA_Yang10
#' @param Rrs665 Rrs665
#' @param Rrs708 Rrs708
#' @param Rrs753 Rrs753
#' @export
FBA_Yang10 <- function(Rrs665, Rrs708, Rrs753){
  return(161.24*(1/Rrs665-1/Rrs708)/(1/Rrs753-1/Rrs708)+28.04)
}

#' @title SCI_Shen10
#' @param Rrs560 Rrs560
#' @param Rrs620 Rrs620
#' @param Rrs665 Rrs665
#' @param Rrs681 Rrs681
#' @export
SCI_Shen10 <- function(Rrs560, Rrs620, Rrs665, Rrs681){
  H_Chl <- (Rrs681+(681-665)/(681-620)*(Rrs620-Rrs681)) - Rrs665
  H_deta <- Rrs620 - (Rrs681+(681-620)/(681-560)*(Rrs560-Rrs681))
  SCI <- H_Chl - H_deta
  return(0.057*SCI^(-0.6327))
}

#' @title Gons08
#' @param Rrs665 Rrs665
#' @param Rrs709 Rrs709
#' @param Rrs779 Rrs779
#' @export
Gons08 <- function(Rrs665, Rrs709, Rrs779){
  bbp <- 1.61 * Rrs779 / (0.082 - 0.6*Rrs779)
  return((Rrs709/Rrs665*(0.7+bbp)-0.4-bbp^1.06)/0.016)
}

#' @title Chla_algorithms_name
#' @export
Chla_algorithms_name <- function(){
  return(c('BR_Gil10','BR_Git11',
           'TBA_Gil10','TBA_Git11',
           'C6','OC4E','NDCI_Mi12', 'Gons08',
           'FBA_Le13','FBA_Yang10',
           'SCI_Shen10'))
}

#' @title Assessment each algorithm for every cluster
#' @name Assessment_via_cluster
#' @param pred prediciton of Chla
#' @param meas in-situ measurement of Chla
#' @param memb membership value matrix
#' @param metrics metrics need to be calculated
#' @param total logical
#' @param hard.mode hard.mode
#' @export
#' @return List
Assessment_via_cluster <- function(pred, meas, memb,
                                   metrics = c('MAE2','MAPE'),
                                   total = TRUE,
                                   hard.mode= TRUE){
  if(nrow(pred) != length(meas) | nrow(pred) != nrow(memb))
    stop('Rows of input are different!')
  if(anyNA(pred) | anyNA(meas) | anyNA(memb))
    stop('Including NA values!')
  
  for(i in 1:length(metrics))
    metrics[i] <- match.arg(metrics[i], cal.metrics.names())
  
  # generate the output dataframe  
  model_names <- colnames(pred)
  cluster_names <- colnames(memb)
  cluster_crisp <- apply(memb,1,which.max) %>% sprintf('M%s',.)
  
  validation <- matrix(data=0,
                   nrow=length(cluster_names), 
                   ncol=length(model_names)) %>% as.data.frame()
  colnames(validation) <- model_names
  rownames(validation) <- cluster_names
  
  # output is a list
  result <- list()
  for(i in 1:length(metrics))
    result[[metrics[i]]] <- validation
  
  # strat loop via model and cluster
  for(model in model_names){
    for(cluster in cluster_names){
      
      x <- meas[cluster_crisp %in% cluster] # true
      y <- pred[cluster_crisp %in% cluster, model] # pred
      
      for(metric in names(result)){
        result[[metric]][which(cluster_names == cluster),
                         which(model_names == model)] <- cal.metrics(x,y,metric)
      }
    }
    if(total == TRUE){
      
      x <- meas
      y <- pred[, model]
      
      for(metric in names(result)){
        
        if(rownames(result[[metric]])[nrow(result[[metric]])] != 'SUM'){
          result[[metric]] %<>% rbind(.,NA)
          rownames(result[[metric]])[nrow(result[[metric]])] <- 'SUM'
        }
        result[[metric]]['SUM',
                         which(model_names == model)] <- cal.metrics(x,y,metric)  
      }
    }
  }
  return(result)
}






