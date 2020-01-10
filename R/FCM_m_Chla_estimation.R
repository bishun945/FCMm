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
  bind.Chla <- data.frame(BR=.BR(Rrs709,Rrs665),
                          TBA=.TBA(Rrs665,Rrs709,Rrs754),
                          C6=.C6(Rrs665,Rrs754),
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

.BR <- function(Rrs709, Rrs665){
  return(abs(35.75*Rrs709/Rrs665-19.3)^1.124)
}

.TBA <- function(Rrs665, Rrs709, Rrs754){
  return(abs(113.36*(1/Rrs665-1/Rrs709)*Rrs754+16.45)^1.124)
}

.C6 <- function(Rrs665, Rrs754){
  return(10^( 1/Rrs665*Rrs754 * 0.14 + 2.11))
}
