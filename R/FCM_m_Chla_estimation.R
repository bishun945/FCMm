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
#' @family Algorithms: Chla concentration

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
#' @param Rrs665 Rrs665
#' @param Rrs709 Rrs709
#' @export
#' @family Algorithms: Chla concentration
BR_Gil10 <- function(Rrs665, Rrs709){
  return(abs(35.75*Rrs709/Rrs665-19.3)^1.124)
}

#' @title TBA_Gil10
#' @param Rrs665 Rrs665
#' @param Rrs709 Rrs709
#' @param Rrs754 Rrs754
#' @export
#' @family Algorithms: Chla concentration
TBA_Gil10 <- function(Rrs665, Rrs709, Rrs754){
  return(abs(113.36*(1/Rrs665-1/Rrs709)*Rrs754+16.45)^1.124)
}

#' @title C6
#' @param Rrs665 Rrs665
#' @param Rrs754 Rrs754
#' @export
#' @family Algorithms: Chla concentration
C6 <- function(Rrs665, Rrs754){
  return(10^( 1/Rrs665*Rrs754 * 0.14 + 2.11))
}

#' @title OC4E
#' @param Rrs443 Rrs443
#' @param Rrs490 Rrs490
#' @param Rrs510 Rrs510
#' @param Rrs560 Rrs560
#' @export
#' @family Algorithms: Chla concentration
OC4E <- function(Rrs443, Rrs490, Rrs510, Rrs560){
  X <- apply(cbind(Rrs443,Rrs490,Rrs510),1,max)/Rrs560 %>% log10
  return(10^(0.3255-2.7677*X+2.4409*X^2-1.1288*X^3-0.4990*X^4))
}

#' @title BR_Git11
#' @param Rrs665 Rrs665
#' @param Rrs709 Rrs709
#' @export
#' @family Algorithms: Chla concentration
BR_Git11 <- function(Rrs665, Rrs709){
  return(72.66*Rrs709/Rrs665-46.535)
}

#' @title TBA_Git11
#' @param Rrs665 Rrs665
#' @param Rrs709 Rrs709
#' @param Rrs754 Rrs754
#' @export
#' @family Algorithms: Chla concentration
TBA_Git11 <- function(Rrs665, Rrs709, Rrs754){
  return(243.86*(1/Rrs665-1/Rrs709)*Rrs754+23.17)
}

#' @title NDCI_Mi12
#' @param Rrs665 Rrs665
#' @param Rrs709 Rrs709
#' @export
#' @family Algorithms: Chla concentration
NDCI_Mi12 <- function(Rrs665, Rrs709){
  ind <- (Rrs709-Rrs665)/(Rrs709+Rrs665)
  return(14.039+86.115*ind+194.325*ind^2)
}

#' @title FBA_Le13
#' @param Rrs665 Rrs665
#' @param Rrs681 Rrs681
#' @param Rrs709 Rrs709
#' @export
#' @family Algorithms: Chla concentration
FBA_Le13 <- function(Rrs665, Rrs681, Rrs709){
  return(18.492*(1/Rrs665-1/Rrs681)/(1/Rrs709-1/Rrs681)+6.1302)
}

#' @title FBA_Yang10
#' @param Rrs665 Rrs665
#' @param Rrs709 Rrs709
#' @param Rrs754 Rrs754
#' @export
#' @family Algorithms: Chla concentration
FBA_Yang10 <- function(Rrs665, Rrs709, Rrs754){
  return(161.24*(1/Rrs665-1/Rrs709)/(1/Rrs754-1/Rrs709)+28.04)
}

#' @title SCI_Shen10
#' @param Rrs560 Rrs560
#' @param Rrs620 Rrs620
#' @param Rrs665 Rrs665
#' @param Rrs681 Rrs681
#' @export
#' @family Algorithms: Chla concentration
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
#' @family Algorithms: Chla concentration
Gons08 <- function(Rrs665, Rrs709, Rrs779){
  bbp <- 1.61 * Rrs779 / (0.082 - 0.6*Rrs779)
  return((Rrs709/Rrs665*(0.7+bbp)-0.4-bbp^1.06)/0.016)
}

#' @title TC2
#' @param Rrs443 Rrs443
#' @param Rrs560 Rrs560
#' @param Rrs665 Rrs665
#' @param Rrs709 Rrs709
#' @param Rrs754 Rrs754
#' @export
#' @return A list with names as Chla_final, Chla_clean, Chla_turbid, and flag
#' @family Algorithms: Chla concentration
TC2 <- function(Rrs443, Rrs560, Rrs665, Rrs709, Rrs754){

  Chla_final <- Chla_clean <- TC2_clean(Rrs443, Rrs560, Rrs665, Rrs709)
  flag <- rep('Clean',length(Chla_final))
  Chla_turbid <- TC2_turbid(Rrs443, Rrs560, Rrs665, Rrs754)
  MCI <- Rrs709 - Rrs665 - (Rrs754-Rrs665) * (709-665) / (754-665)
  w <- which(MCI > 0.0016)
  Chla_final[w] <- Chla_turbid[w]
  flag[w] <- 'Turbid'
  return(list(Chla_final=Chla_final,
              Chla_clean=Chla_clean,
              Chla_turbid=Chla_turbid,
              flag=flag))
}  

#' @title TC2_clean
#' @param Rrs443 Rrs443
#' @param Rrs560 Rrs560
#' @param Rrs665 Rrs665
#' @param Rrs709 Rrs709
#' @export
#' @family Algorithms: Chla concentration
TC2_clean <- function(Rrs443, Rrs560, Rrs665, Rrs709){
  g0 = 0.084
  g1 = 0.170
  lambda_0 <- 709
  yita <- 0.17
  aChla_star <- 0.017
  Rrs <- cbind(Rrs443, Rrs560, Rrs665, Rrs709)
  rrs <- Rrs / (0.52 + 1.7 * Rrs)
  u <- (-g0 + sqrt(4 * g1 * rrs)) / (2 * g1)
  bbp_0 <- u[,4] * dt_water$aw[dt_water$nm == lambda_0] / (1-u[,4]) -
    dt_water$bbw[dt_water$nm == lambda_0]
  Y <- 2.0 * (1 - 1.2 * exp(-0.9 * rrs[,1] / rrs[,2]))
  bb560 <- bbp_0 * (lambda_0/560)^Y + dt_water$bbw[dt_water$nm == 560]
  bb665 <- bbp_0 * (lambda_0/665)^Y + dt_water$bbw[dt_water$nm == 665]
  anw560 <- (1-u[,2])*bb560 / u[,2] - dt_water$aw[dt_water$nm == 560]
  anw665 <- (1-u[,3])*bb665 / u[,3] - dt_water$aw[dt_water$nm == 665]
  aph665 <- yita * anw560 + (1-yita) * anw665
  Chla <- aph665 / aChla_star
  return(Chla)
}


#' @title TC2_turbid
#' @param Rrs443 Rrs443
#' @param Rrs560 Rrs560
#' @param Rrs665 Rrs665
#' @param Rrs754 Rrs754
#' @export
#' @family Algorithms: Chla concentration 
TC2_turbid <- function(Rrs443, Rrs560, Rrs665, Rrs754){
  g0 = 0.084
  g1 = 0.170
  lambda_0 <- 754
  yita <- 0.17
  aChla_star <- 0.017
  Rrs <- cbind(Rrs443, Rrs560, Rrs665, Rrs754)
  rrs <- Rrs / (0.52 + 1.7 * Rrs)
  u <- (-g0 + sqrt(4 * g1 * rrs)) / (2 * g1)
  bbp_0 <- u[,4] * dt_water$aw[dt_water$nm == lambda_0] / (1-u[,4]) -
    dt_water$bbw[dt_water$nm == lambda_0]
  Y <- 2.0 * (1 - 1.2 * exp(-0.9 * rrs[,1] / rrs[,2]))
  bb560 <- bbp_0 * (lambda_0/560)^Y + dt_water$bbw[dt_water$nm == 560]
  bb665 <- bbp_0 * (lambda_0/665)^Y + dt_water$bbw[dt_water$nm == 665]
  anw560 <- (1-u[,2])*bb560 / u[,2] - dt_water$aw[dt_water$nm == 560]
  anw665 <- (1-u[,3])*bb665 / u[,3] - dt_water$aw[dt_water$nm == 665]
  aph665 <- yita * anw560 + (1-yita) * anw665
  Chla <- aph665 / aChla_star
  return(Chla)
}


#' @title Chla_algorithms_name
#' @export
#' @family Algorithms: Chla concentration
Chla_algorithms_name <- function(){
  return(c('BR_Gil10','BR_Git11',
           'TBA_Gil10','TBA_Git11',
           'C6','OC4E','NDCI_Mi12', 'Gons08',
           'FBA_Le13','FBA_Yang10',
           'SCI_Shen10', 'TC2','TC2_turbid','TC2_clean'))
}

#' @title run_all_Chla_algorithms
#' @param Rrs Dataframe that should with required colnames 443, 490, 510, 560, 620, 665, 681, 709, 754, 779
#' @param wv_range Number that used to define the range of wavelength to capture
#'   the center wavelength of required band
#' @details Please type \code{Chal_algorithm_name()} to see all Chla algorithms
#' @export
#' @return A list
#' @family Algorithms: Chla concentration
#' 
run_all_Chla_algorithms <- function(Rrs, wv_range=3){
  wv <- str_subset(names(Rrs), '\\d') %>% as.numeric
  wv_need <- c(443, 490, 510, 560, 620, 665, 681, 709, 754, 779)
  for(i in wv_need){
    text <- sprintf('Rrs%s <- Rrs[,which(abs(wv-%s) <= wv_range)]',i,i)
    # print(text)
    eval(parse(text=text))
  }
  result <- list(
    BR_Gil10 = BR_Gil10(Rrs709, Rrs665),
    TBA_Gil10 = TBA_Gil10(Rrs665, Rrs709, Rrs754),
    C6 = C6(Rrs665, Rrs754),
    OC4E = OC4E(Rrs443, Rrs490, Rrs510, Rrs560),
    BR_Git11 = BR_Git11(Rrs665, Rrs709),
    TBA_Git11 = TBA_Git11(Rrs665, Rrs709, Rrs754),
    NDCI_Mi12 = NDCI_Mi12(Rrs665, Rrs709),
    FBA_Le13 = FBA_Le13(Rrs665, Rrs681, Rrs709),
    FBA_Yang10 = FBA_Yang10(Rrs665, Rrs709, Rrs754),
    SCI_Shen10 = SCI_Shen10(Rrs560, Rrs620, Rrs665, Rrs681),
    Gons08 = Gons08(Rrs665, Rrs709, Rrs779),
    TC2_final = TC2(Rrs443, Rrs560, Rrs665, Rrs709, Rrs754)$Chla_final,
    TC2_clean = TC2(Rrs443, Rrs560, Rrs665, Rrs709, Rrs754)$Chla_clean,
    TC2_turbid = TC2(Rrs443, Rrs560, Rrs665, Rrs709, Rrs754)$Chla_turbid)
  w_zero <- lapply(result, length) %>% as.matrix %>% {which(. == 0)}
  for(i in w_zero){
    result[[i]] <- NA
  }
  return(result)
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
#' @family Algorithm assessment
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
      
      w <- str_extract(cluster_crisp,'\\d') %in% str_extract(cluster,"\\d")
      x <- meas[w] # true
      y <- pred[w, model] # pred
      
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






