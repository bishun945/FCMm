#' @name FCM_m_Chla_estimation
#' @title Estimate Chla concentration by algorithms blending via membership values
#' @description
#' To calculate the Chla concentration via blending remote sensing algorithms.
#' @usage FCM_m_Chla_estimation(Rrs, U)
#' @param Rrs Data.frame of Remote sensing reflectance which should
#'   contain colname with (also at) \code{Rrs665}, \code{Rrs709},
#'   and \code{Rrs754} with unit sr^-1.
#' @param U Data.frame of membership values which should have seven
#'   columns presenting seven membership values from each cluster produced by FCM-m
#' @return Each column presents Chla concentration with unit ug/L or mg/m^3 by model
#'   BR, TBA, C6 and blending result. Also with the membership values of each cluster.
#' @note The input of \code{Rrs} must have bands with wavelength at 665, 709 and 754 nm.
#'   See examples of using this function.
#'   (2020-02-28) The \code{C6} model was replaced by \code{Bloom} model.
#' @export
#' @examples
#' 
#' library(FCMm)
#' data("WaterSpec35")
#' data("Bi_clusters")
#' Rrs <- WaterSpec35[,3:17]
#' result <- apply_FCM_m(Rrs=Rrs, option.plot=TRUE)
#' dt_Chla <- FCM_m_Chla_estimation(Rrs=data.frame(Rrs665=Rrs$`665`, 
#'   Rrs709=Rrs$`708.75`,Rrs754=Rrs$`753.75`),U=result$u)
#' 
#' @references
#' \itemize{
#'   \item Bi S, Li Y, Xu J, et al. Optical classification of inland waters based on
#'     an improved Fuzzy C-Means method[J]. Optics Express, 2019, 27(24): 34838-34856.
#'   \item Gilerson A A, Gitelson A A, Zhou J, et al. Algorithms for remote estimation 
#'     of chlorophyll-a in coastal and inland waters using red and near infrared bands[J].
#'     Optics Express, 2010, 18(23): 24109-24125.
#' }
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
  bind.Chla <- data.frame(
    BR=BR_Gil10(Rrs665,Rrs709)$Chla,
    TBA=TBA_Gil10(Rrs665,Rrs709,Rrs754)$Chla,
    # BR=BR_Git11(Rrs709,Rrs665)$Chla,
    # TBA=TBA_Git11(Rrs665,Rrs709,Rrs754)$Chla,
    Bloom=Bloom(Rrs665,Rrs754)$Chla,
    M=round(U,4)
    )
  # TBA: C1 C2 C5
  # BR: C3 C4 C7
  # Bloom: C6
  U3 <- data.frame(u.TBA=bind.Chla[,c('M.1','M.2','M.5')] %>% apply(.,1,sum),
                   u.BR =bind.Chla[,c('M.3','M.4','M.7')] %>% apply(.,1,sum),
                   u.Bloom =bind.Chla[,c('M.6')])
  CONC3 <- bind.Chla[,c('TBA','BR','Bloom')]
  bind.Chla$conc.Blend <- apply(U3*CONC3, 1, sum)

  return(bind.Chla)
}



#' @name Chla_Gil10_Git11
#' @title Band-Ratio and Three-bands Chla algorithm by Gilerson et al. (2010) and Gitelson et al. (2011)
#' @param Rrs665 Rrs665
#' @param Rrs709 Rrs709
#' @param Rrs754 Rrs754
#' @export
#' @return A list including the algorithm index and Chla concentration.
#' @family Algorithms: Chla concentration
#' @examples 
#' data(WaterSpec35)
#' res_TBA_Gil10 = TBA_Gil10(WaterSpec35$`665`, WaterSpec35$`708.75`, WaterSpec35$`753.75`)
#' res_BR_Gil10  = BR_Gil10(WaterSpec35$`665`, WaterSpec35$`708.75`)
#' res_TBA_Git11 = TBA_Git11(WaterSpec35$`665`, WaterSpec35$`708.75`, WaterSpec35$`753.75`)
#' res_BR_Git11 = BR_Git11(WaterSpec35$`665`, WaterSpec35$`708.75`)
#' @references 
#' \itemize{
#'   \item Gilerson A A, Gitelson A A, Zhou J, et al. Algorithms for remote estimation of 
#'     chlorophyll-a in coastal and inland waters using red and near infrared bands[J].
#'     Optics Express, 2010, 18(23): 24109-24125.
#'   \item Gitelson A A, Gurlin D, Moses W J, et al. Remote Estimation of Chlorophyll-a 
#'     Concentration[J]. Advances in environmental remote sensing: sensors, algorithms, 
#'     and applications, 2011: 439.
#' }
#'   
TBA_Gil10 <- function(Rrs665, Rrs709, Rrs754){
  ind = (1/Rrs665-1/Rrs709)*Rrs754
  Chla = (113.36*ind+16.45)^1.124
  result <- list(ind=ind,
                 Chla=Chla)
  return(result)
}

#' @rdname Chla_Gil10_Git11
#' @export
BR_Gil10 <- function(Rrs665, Rrs709){
  ind = (Rrs709/Rrs665)
  Chla = (35.75*ind-19.3)^1.124
  result <- list(ind=ind,
                 Chla=Chla)
  return(result)
}

#' @rdname Chla_Gil10_Git11
#' @export
TBA_Git11 <- function(Rrs665, Rrs709, Rrs754){
  ind = (1/Rrs665-1/Rrs709)*Rrs754
  Chla = 243.86*ind+23.17
  result <- list(ind=ind,
                 Chla=Chla)
  return(result)
}

#' @rdname Chla_Gil10_Git11
#' @export
BR_Git11 <- function(Rrs665, Rrs709){
  ind = Rrs709/Rrs665
  Chla = 72.66*ind-46.535
  result <- list(ind=ind,
                 Chla=Chla)
  return(result)
}



#' @name C6
#' @title C6 (Being deprecated)
#' @param Rrs665 Rrs665
#' @param Rrs754 Rrs754
#' @export
#' @return A list including the algorithm index and Chla concentration.
#' @family Algorithms: Chla concentration
#' @examples 
#' data(WaterSpec35)
#' res = Bloom(WaterSpec35$`665`, WaterSpec35$`753.75`)
C6 <- function(Rrs665, Rrs754){
  ind <- 1/Rrs665*Rrs754
  Chla <- 10^( ind * 0.14 + 2.11)
  result <- list(ind=ind,
                 Chla=Chla)
  return(result)
}



#' @name Bloom
#' @title Chla estimation algorithm for the bloom water type
#' @param Rrs665 Rrs665
#' @param Rrs754 Rrs754
#' @export
#' @return A list including the algorithm index and Chla concentration.
#' @family Algorithms: Chla concentration
#' @examples 
#' data(WaterSpec35)
#' res = Bloom(WaterSpec35$`665`, WaterSpec35$`753.75`)
Bloom <- function(Rrs665, Rrs754){
  ind <- 1/Rrs665*Rrs754
  Chla <- 257.740*exp(0.279*ind)
  result <- list(ind=ind,
                 Chla=Chla)
  return(result)
}



#' @name Chla_OC_OLCI
#' @title NASA standard ocean color algorithm (version 4, 5, and 6) for the Ocean and Land 
#'   Color Instrument (OLCI) bands
#' @param Rrs412 Rrs412
#' @param Rrs443 Rrs443
#' @param Rrs490 Rrs490
#' @param Rrs510 Rrs510
#' @param Rrs560 Rrs560
#' @param Rrs665 Rrs665
#' @export
#' @return A list including the algorithm index and Chla concentration.
#' @family Algorithms: Chla concentration
#' @examples 
#' data(WaterSpec35)
#' res_OC4 = OC4_OLCI(WaterSpec35$`442.5`, WaterSpec35$`490`, 
#' WaterSpec35$`510`, WaterSpec35$`560`)
#' res_OC5 = OC5_OLCI(WaterSpec35$`412.5`, WaterSpec35$`442.5`, 
#' WaterSpec35$`490`, WaterSpec35$`510`, WaterSpec35$`560`)
#' res_OC6 = OC6_OLCI(WaterSpec35$`412.5`, WaterSpec35$`442.5`, 
#' WaterSpec35$`490`, WaterSpec35$`510`, WaterSpec35$`560`, WaterSpec35$`665`)
#' 
#' @references O'Reilly J E, Werdell P J. Chlorophyll algorithms for ocean color
#'   sensors-OC4, OC5 & OC6[J]. Remote sensing of environment, 2019, 229: 32-47.
#'   
OC4_OLCI <- function(Rrs443, Rrs490, Rrs510, Rrs560){
  X <- apply(cbind(Rrs443,Rrs490,Rrs510),1,max)/Rrs560
  X <- log10(X)
  Chla <- 10^(0.42540-3.21679*X+2.86907*X^2-0.62628*X^3-1.09333*X^4)
  result <- list(X=X,
                 Chla=Chla)
  return(result)
}

#' @rdname Chla_OC_OLCI
#' @export
OC5_OLCI <- function(Rrs412, Rrs443, Rrs490, Rrs510, Rrs560){
  X <- apply(cbind(Rrs412,Rrs443,Rrs490,Rrs510),1,max)/Rrs560
  X <- log10(X)
  Chla <- 10^(0.43213-3.13001*X+3.05479*X^2-1.45176*X^3-0.24947*X^4) 
  result <- list(X=X,
                 Chla=Chla)
  return(result)
}

#' @rdname Chla_OC_OLCI
#' @export
OC6_OLCI <- function(Rrs412, Rrs443, Rrs490, Rrs510, Rrs560, Rrs665){
  X <- apply(cbind(Rrs412,Rrs443,Rrs490,Rrs510),1,max)/apply(cbind(Rrs560,Rrs665),1,mean)
  X <- log10(X)
  Chla <- 10^(0.95039-3.05404*X+2.17992*X^2-1.12097*X^3-0.15262*X^4)
  result <- list(X=X,
                 Chla=Chla)
  return(result)
}

#' @name OCI_Hu12
#' @title Ocean Color Index (OCI) algorithm by Hu et al. (2012)
#' @param Rrs443 Rrs443
#' @param Rrs490 Rrs490
#' @param Rrs510 Rrs510
#' @param Rrs560 Rrs560
#' @param Rrs665 Rrs665
#' @export
#' @return A list including the algorithm indices and Chla concentration.
#' @family Algorithms: Chla concentration
#' @examples 
#' data(WaterSpec35)
#' res = OCI_Hu12(WaterSpec35$`442.5`, WaterSpec35$`490`, WaterSpec35$`510`, 
#' WaterSpec35$`560`, WaterSpec35$`665`)
#' @references Hu C, Lee Z, Franz B. Chlorophyll a algorithms for oligotrophic oceans: 
#'   A novel approach based on three band reflectance difference[J]. Journal of Geophysical 
#'   Research: Oceans, 2012, 117(C1).
OCI_Hu12 <- function(Rrs443, Rrs490, Rrs510, Rrs560, Rrs665){
  CI <- Rrs560 - (Rrs443+(560-443)/(665-443)*(Rrs665-Rrs443))
  Chl_CI <- 10^(-0.4909+191.6590*CI)
  X <- log10(apply(cbind(Rrs443, Rrs490, Rrs510), 1, max)/Rrs560)
  Chl_OC4 <- 10^(0.3272-2.9940*X+2.7218*X^2-1.2259*X^3-0.5683*X^4)
  a <- (Chl_CI - 0.25)/(0.3 - 0.25)
  b <- (0.3 - Chl_CI)/(0.3 - 0.25)
  Chl_OCI <- Chl_CI
  Chl_OCI[Chl_CI <= 0.25] <- Chl_CI[Chl_CI <= 0.25]
  Chl_OCI[Chl_CI >  0.3 ] <- Chl_OC4[Chl_CI >  0.3 ]
  w <- which(Chl_CI > 0.25 & Chl_CI <= 0.3)
  tmp <- a*Chl_OC4 + b*Chl_CI
  Chl_OCI[w] <- tmp[w]
  
  result <- list(
    CI      = CI,
    X       = X,
    Chl_CI  = Chl_CI,
    Chl_OC4 = Chl_OC4,
    Chl_OCI = Chl_OCI
  )
  
  return(result)
}


#' @name NDCI_Mi12
#' @title Normalized Difference Chlorophyll Index (NDCI) algorithm by Mishra et al. (2012)
#' @param Rrs665 Rrs665
#' @param Rrs709 Rrs709
#' @export
#' @return A list including the algorithm index and Chla concentration.
#' @family Algorithms: Chla concentration
#' @note Chla_test as a test version of NDCI_Mi12 results parameterized by Bi
#' @examples
#' data(WaterSpec35)
#' res = NDCI_Mi12(WaterSpec35$`665`, WaterSpec35$`708.75`)
#' @references Mishra S, Mishra D R. Normalized difference chlorophyll index:
#'   A novel model for remote estimation of chlorophyll-a concentration in turbid 
#'   productive waters[J]. Remote Sensing of Environment, 2012, 117: 394-406.
NDCI_Mi12 <- function(Rrs665, Rrs709){
  NDCI <- (Rrs709-Rrs665)/(Rrs709+Rrs665)
  Chla <- 14.039+86.115*NDCI+194.325*NDCI^2
  Chla_test <- 10^(1.3636 + 3.1257*NDCI - 3.0591*NDCI^2 + 3.6429*NDCI^3)
  result <- list(
    NDCI = NDCI,
    Chla = Chla,
    Chla_test = Chla_test
  )
  return(result)
}



#' @name Chla_FBA
#' @title Four-bands Chla algorithm by Yang et al. (2010) and Le et al. (2013)
#' @param Rrs665 Rrs665
#' @param Rrs681 Rrs681
#' @param Rrs709 Rrs709
#' @param Rrs754 Rrs754
#' @export
#' @return A list including the algorithm index and Chla concentration.
#' @family Algorithms: Chla concentration
#' @examples 
#' data(WaterSpec35)
#' res_L = FBA_Le13(WaterSpec35$`665`, WaterSpec35$`681.25`, WaterSpec35$`708.75`)
#' res_Y = FBA_Yang10(WaterSpec35$`665`, WaterSpec35$`708.75`, WaterSpec35$`753.75`)
#' @references 
#' \itemize{
#'   \item Le C, Hu C, Cannizzaro J, et al. Evaluation of chlorophyll-a remote
#'   sensing algorithms for an optically complex estuary[J]. Remote Sensing of Environment, 
#'   2013, 129: 75-89.
#'   \item Yang W, Matsushita B, Chen J, et al. An enhanced three-band index for estimating 
#'   chlorophyll-a in turbid case-II waters: case studies of Lake Kasumigaura, Japan, and Lake 
#'   Dianchi, China[J]. IEEE Geoscience and Remote Sensing Letters, 2010, 7(4): 655-659.
#' }
#'   
FBA_Le13 <- function(Rrs665, Rrs681, Rrs709){
  ind = (1/Rrs665-1/Rrs681)*(1/Rrs709-1/Rrs681)
  Chla = 18.492*ind+6.1302
  result <- list(ind=ind,
                 Chla=Chla)
  return(result)
}

#' @rdname Chla_FBA
#' @export
FBA_Yang10 <- function(Rrs665, Rrs709, Rrs754){
  ind = (1/Rrs665-1/Rrs709)/(1/Rrs754-1/Rrs709)
  Chla = 161.24*ind+28.04
  result <- list(ind=ind,
                 Chla=Chla)
  return(result)
}



#' @name SCI_Shen10
#' @title Synthetic Chlorophyll Index (SCI) algorithm by Shen et al. (2010)
#' @param Rrs560 Rrs560
#' @param Rrs620 Rrs620
#' @param Rrs665 Rrs665
#' @param Rrs681 Rrs681
#' @export
#' @return A list including the algorithm index and Chla concentration.
#' @family Algorithms: Chla concentration
#' @examples 
#' data(WaterSpec35)
#' res = SCI_Shen10(WaterSpec35$`560`, WaterSpec35$`620`, 
#'   WaterSpec35$`665`, WaterSpec35$`681.25`)
#' @references Shen F, Zhou Y X, Li D J, et al. Medium resolution imaging spectrometer 
#'   (MERIS) estimation of chlorophyll-a concentration in the turbid sediment-laden waters 
#'   of the Changjiang (Yangtze) Estuary[J]. International Journal of Remote Sensing, 2010, 
#'   31(17-18): 4635-4650.
SCI_Shen10 <- function(Rrs560, Rrs620, Rrs665, Rrs681){
  H_Chl <- (Rrs681+(681-665)/(681-620)*(Rrs620-Rrs681)) - Rrs665
  H_deta <- Rrs620 - (Rrs681+(681-620)/(681-560)*(Rrs560-Rrs681))
  SCI <- H_Chl - H_deta
  result <- list(
    H_Chl      = H_Chl,
    H_deta     = H_deta,
    SCI        = SCI,
    Chla_EPCHC = 0.057*SCI^(-0.6327), # EPCHC
    Chla_SPR   = 179378*SCI^2+92.934*SCI+0.2736*SCI, # Shen in Spring
    Chla_SUM   = 550383*SCI^2+2769*SCI+4.3866
  )
  return(result)
}


#' @noRd
MPH_Mat14 <- function(Rrs620, Rrs665, Rrs681, Rrs709, Rrs754, Rrs885) {
  
  dt_MPH <- data.frame(Rrs665, Rrs681, Rrs709, Rrs754, Rrs885) %>%
    setNames(., c(665, 681, 709, 754, 885))
  IND_MPH <- apply(dt_MPH, 1, function(x) {
    tmp <- x[c("681", "709", "754")]
    w_max <- which.max(x[c("681", "709", "754")])
    lamb_max <- as.numeric(names(tmp)[w_max])
    Rrs_max  <- as.numeric(tmp[w_max])
    ind <- Rrs_max - x["665"] - ((x["885"]-x["665"])*(lamb_max-665)/(885-665))
    return(ind)
  })
  SICF_peak <- Rrs681 - Rrs665 - ((Rrs709-Rrs665)*(681-665)/(709-665))
  SIPAF_peak <- Rrs665 - Rrs620 - ((Rrs681-Rrs620)*(665-620)/(681-620))
  cyano_flag <- rep(FALSE, length(SICF_peak))  
  cyano_flag[SICF_peak < 0 & SIPAF_peak > 0] <- TRUE
  return(list(
    IND_MPH = IND_MPH,
    SICF_peak = SICF_peak,
    SIPAF_peak = SIPAF_peak,
    cyano_flag = cyano_flag
  ))
}



#' @name Gons08
#' @title Gons algorithm by Gons et al. (2008)
#' @param Rrs665 Rrs665
#' @param Rrs709 Rrs709
#' @param Rrs779 Rrs779
#' @export
#' @return A list including the algorithm estimated bbp and Chla concentration.
#' @family Algorithms: Chla concentration
#' @examples 
#' data(WaterSpec35)
#' res = Gons08(WaterSpec35$`665`, WaterSpec35$`708.75`, WaterSpec35$`778.75`)
#' @references Gons H J, Auer M T, Effler S W. MERIS satellite chlorophyll mapping of 
#'   oligotrophic and eutrophic waters in the Laurentian Great Lakes[J]. Remote Sensing 
#'   of Environment, 2008, 112(11): 4098-4106.
Gons08 <- function(Rrs665, Rrs709, Rrs779){
  bbp <- 1.61 * Rrs779 / (0.082 - 0.6*Rrs779)
  Chla <- (Rrs709/Rrs665*(0.7+bbp)-0.4-bbp^1.06)/0.016
  result <- list(
    bbp  = bbp,
    Chla = Chla
  )
  return(result)
}



#' @name Smith18
#' @title Smith algorithm by Smith et al. (2018)
#' @param Rrs443 Rrs443
#' @param Rrs490 Rrs490
#' @param Rrs510 Rrs510
#' @param Rrs560 Rrs560
#' @param Rrs665 Rrs665
#' @param Rrs709 Rrs709
#' @export
#' @return A list including the algorithm estimated CHL_G2B, CHL_OCI, 
#'   CHL_BLEND concentration. PHI index. a1 and a2 coefficients.
#' @family Algorithms: Chla concentration
#' @examples 
#' data(WaterSpec35)
#' res = Smith18(WaterSpec35$`442.5`, WaterSpec35$`490`, WaterSpec35$`510`, 
#' WaterSpec35$`560`, WaterSpec35$`665`, WaterSpec35$`708.75`)
#' @references Smith M E, Lain L R, Bernard S. An optimized chlorophyll a switching 
#'   algorithm for MERIS and OLCI in phytoplankton-dominated waters[J]. 
#'   Remote Sensing of Environment, 2018, 215: 217-227.
Smith18 <- function(Rrs443, Rrs490, Rrs510, Rrs560, Rrs665, Rrs709) {
  CHL_G2B = BR_Gil10(Rrs665, Rrs709)$Chla
  CHL_OCI = OCI_Hu12(Rrs443, Rrs490, Rrs510, Rrs560, Rrs665)$Chl_OCI
  r   = Rrs709 / Rrs665
  PHI = r * NA
  PHI[r < 0.75] = 0.75
  PHI[r > 1.15] = 1.15
  PHI[r <= 1.15 & r >= 0.75] = r[r <= 1.15 & r >= 0.75]
  a1 = (PHI - 0.75) / (1.15 - 0.75)
  a2 = abs(a1 - 1)
  CHL_BLEND <- a1 * CHL_G2B + a2 * CHL_OCI
  result <- list(
    CHL_G2B   = CHL_G2B,
    CHL_OCI   = CHL_OCI,
    CHL_BLEND = CHL_BLEND,
    PHI = PHI,
    a1  = a1,
    a2  = a2
  )
  return(result)
}




#' @name TC2
#' @title Chla concentration algorithm for Turbid Case-2 (TC2) waters by Liu et al. (2020)
#' @param Rrs443 Rrs443
#' @param Rrs560 Rrs560
#' @param Rrs665 Rrs665
#' @param Rrs709 Rrs709
#' @param Rrs754 Rrs754
#' @export
#' @return A list including MCI index, Chla_final, Chla_clean, Chla_turbid, and flag.
#' @family Algorithms: Chla concentration
#' @examples 
#' data(WaterSpec35)
#' res = TC2(WaterSpec35$`442.5`, WaterSpec35$`560`, WaterSpec35$`665`,
#' WaterSpec35$`708.75`, WaterSpec35$`753.75`)
#' 
#' @references Liu G, Li L, Song K, et al. An OLCI-based algorithm for semi-empirically 
#'   partitioning absorption coefficient and estimating chlorophyll a concentration in 
#'   various turbid case-2 waters[J]. Remote Sensing of Environment, 2020, 239: 111648.
TC2 <- function(Rrs443, Rrs560, Rrs665, Rrs709, Rrs754){
  
  Chla_final <- Chla_clean <- TC2_clean(Rrs443, Rrs560, Rrs665, Rrs709)
  flag <- rep('Clean',length(Chla_final))
  Chla_turbid <- TC2_turbid(Rrs443, Rrs560, Rrs665, Rrs709, Rrs754)
  MCI <- Rrs709 - Rrs665 - (Rrs754-Rrs665) * (709-665) / (754-665)
  w <- which(MCI > 0.0016)
  Chla_final[w] <- Chla_turbid[w]
  flag[w] <- 'Turbid'
  return(list(MCI=MCI,
              Chla_final=Chla_final,
              Chla_clean=Chla_clean,
              Chla_turbid=Chla_turbid,
              flag=flag))
}  

#' @name TC2_clean
#' @title TC2_clean
#' @param Rrs443 Rrs443
#' @param Rrs560 Rrs560
#' @param Rrs665 Rrs665
#' @param Rrs709 Rrs709
#' @return Chla concentration.
#' @family Algorithms: Chla concentration
#' @noRd
TC2_clean <- function(Rrs443, Rrs560, Rrs665, Rrs709){
  g0 = 0.089
  g1 = 0.125
  lambda_0 <- 709
  yita <- 0.17
  aChla_star <- 0.017
  aw  = c(0.006963, 0.06209, 0.4285, 0.8013, 2.8666) # 443, 560, 665, 709, 754
  bbw = c(0.002178, 0.0008057, 0.0003937, 0.0003023, 0.0002343)
  # aw  = dt_water$aw[dt_water$nm %in% c(443,560,665,709)]
  # bbw = dt_water$bbw[dt_water$nm %in% c(443,560,665,709)]
  Rrs <- cbind(Rrs443, Rrs560, Rrs665, Rrs709)
  rrs <- Rrs / (0.52 + 1.7 * Rrs)
  u <- (-g0 + sqrt(g0^2 + 4 * g1 * rrs)) / (2 * g1)
  bbp_0 <- u[,4] * aw[4] / (1-u[,4]) - bbw[4]
  Y <- 2.0 * (1 - 1.2 * exp(-0.9 * rrs[,1] / rrs[,2]))
  bb560 <- bbp_0 * (lambda_0/560)^Y + bbw[2]
  bb665 <- bbp_0 * (lambda_0/665)^Y + bbw[3]
  bb709 <- bbp_0 * (lambda_0/709)^Y + bbw[3]
  anw560 <- (1-u[,2])*bb560 / u[,2] - aw[2]
  anw665 <- (1-u[,3])*bb665 / u[,3] - aw[3]
  anw709 <- (1-u[,4])*bb709 / u[,4] - aw[4]
  aph665 <- anw665 - yita * anw560 - (1-yita) * anw709
  Chla <- aph665 / aChla_star
  return(Chla)
}


#' @name TC2_turbid
#' @title TC2_turbid
#' @param Rrs443 Rrs443
#' @param Rrs560 Rrs560
#' @param Rrs665 Rrs665
#' @param Rrs709 Rrs709
#' @param Rrs754 Rrs754
#' @return Chla concentration.
#' @family Algorithms: Chla concentration 
#' @noRd
TC2_turbid <- function(Rrs443, Rrs560, Rrs665, Rrs709, Rrs754){
  g0 = 0.089
  g1 = 0.125
  lambda_0 <- 754
  yita <- 0.17
  aChla_star <- 0.017
  aw  = c(0.006963, 0.06209, 0.4285, 0.8013, 2.8666) # 443, 560, 665, 709, 754
  bbw = c(0.002178, 0.0008057, 0.0003937, 0.0003023, 0.0002343)
  # aw  = dt_water$aw[dt_water$nm %in% c(443,560,665,709)]
  # bbw = dt_water$bbw[dt_water$nm %in% c(443,560,665,709)]
  Rrs <- cbind(Rrs443, Rrs560, Rrs665, Rrs709, Rrs754)
  rrs <- Rrs / (0.52 + 1.7 * Rrs)
  u <- (-g0 + sqrt(g0^2 + 4 * g1 * rrs)) / (2 * g1)
  bbp_0 <- u[,5] * aw[5] / (1-u[,5]) - bbw[5]
  Y <- 2.0 * (1 - 1.2 * exp(-0.9 * rrs[,1] / rrs[,2]))
  bb560 <- bbp_0 * (lambda_0/560)^Y + bbw[2]
  bb665 <- bbp_0 * (lambda_0/665)^Y + bbw[3]
  bb709 <- bbp_0 * (lambda_0/709)^Y + bbw[4]
  anw560 <- (1-u[,2])*bb560 / u[,2] - aw[2]
  anw665 <- (1-u[,3])*bb665 / u[,3] - aw[3]
  anw709 <- (1-u[,4])*bb709 / u[,4] - aw[4]
  aph665 <- anw665 - yita * anw560 - (1-yita) * anw709
  Chla <- aph665 / aChla_star
  return(Chla)
}

#' @name QAA_v5
#' @title Quasi-Analytical Algorithm version 5 (QAA_v5) by Lee et al. (2002)
#' @param wv Wavelength vector
#' @param Rrs Rrs matrix
#' @param wv412 Wavelength index of 412 nm
#' @param wv443 Wavelength index of 443 nm
#' @param wv490 Wavelength index of 490 nm
#' @param wv560 Wavelength index of 560 nm
#' @param wv667 Wavelength index of 667 nm
#' @param verbose Print process information (default as FASLE)
#' @export
#' @return Results of \code{QAA_v5} are the list including:
#'   \item{input}{input list including \code{wv}, \code{Rrs}, \code{wv412},
#'     \code{wv443}, \code{wv490}, \code{wv560}, and \code{wv667}}
#'   \item{const}{constant values of \code{g0}, \code{g1}, and \code{lam0}}
#'   \item{Var}{Indirect variables including \code{rrs}, \code{u}, \code{Chi},
#'     \code{Yita}, \code{Zeta}, \code{S}, and \code{Xi}}
#'   \item{IOP}{Inherent optical properties including \code{a0}, \code{bbp0},
#'     \code{ag443}, \code{a}, \code{adg}, \code{aph}, and \code{bbp}}
#'   \item{Chl}{Results of Chl concentration estimation including \code{Chl_Bricaud}, 
#'     \code{Bricaud_A}, \code{Bricaud_B}, \code{Chl_Liu}, \code{aChla665_star}}
#'   \item{residu}{Residual constrains including \code{Upper_lim_for_Rrs667}, \code{Lower_lim_for_Rrs667},
#'     \code{Rrs667_revised}, and \code{Rrs667_final}}
#' @family Algorithms: Chla concentration 
#' @examples 
#' data(WaterSpec35)
#' Rrs = WaterSpec35[,-c(1,2)]
#' wv = as.numeric(names(Rrs))
#' print(wv)
#' res = QAA_v5(wv, Rrs, wv412 = 2, wv443 = 3, wv490 = 4, wv560 = 6, wv667 = 9)
#' 
#' @references Lee Z P, Carder K L, Arnone R A. Deriving inherent optical properties from water color: 
#'   a multiband quasi-analytical algorithm for optically deep waters[J]. Applied optics, 2002, 41(27): 
#'   5755-5772.
QAA_v5 <- function(wv, Rrs,
                   wv412=NULL, wv443=NULL, wv490=NULL, wv560=NULL, wv667=NULL,
                   verbose=FALSE){
  
  if(anyNA(Rrs) | anyNA(wv))
    stop("NA values in wv or Rrs as input, please check!")
  
  if(length(wv) != ncol(Rrs))
    stop("Input parameter wv has different length with column of Rrs!")
  
  if(verbose)
    message("QAA_v5 is running!")
  
  if(is.null(wv443) | is.null(wv490) | is.null(wv560))
    stop("Please assign the required bands in Rrs")
  
  # Constants
  g0   = 0.089
  g1   = 0.125
  lam0 = 560
  
  # Variable
  rrs   = Rrs / (0.52 + 1.72 * Rrs)
  u     = (-g0 + sqrt(g0^2 + 4 * g1 * rrs)) / (2 * g1)
  Chi   = log10( (rrs[, wv443] + rrs[, wv490]) /
                   (rrs[, wv560] + 5 * rrs[, wv667] / rrs[, wv490] * rrs[, wv667]) )
  a0    = dt_water$aw[dt_water$nm == lam0] + 10^(-1.146 - 1.366*Chi - 0.469*Chi^2)
  bbp0  = u[,wv560] * a0 / (1-u[,wv560]) - dt_water$bbw[dt_water$nm == lam0]
  Yita  = 2.0 * (1 - 1.2*exp(-0.9 * rrs[, wv443] / rrs[, wv560])) # Exponent of bbp
  Zeta  = 0.74 + 0.2 / (0.8 + rrs[, wv443] / rrs[, wv560]) # aph411/aph443
  S     = 0.015 + 0.002 / (0.6 + rrs[, wv443] / rrs[, wv560])
  Xi    = exp(S * (443 - 411))
  
  # Outputs
  bbp   = rrs
  for(i in 1:length(wv))
    bbp[, i] = bbp0 * (lam0/wv[i]) ^ Yita
  
  a     = rrs
  for(i in 1:length(wv))
    a[, i]   = (1-u[,i])/u[,i] * (dt_water$bbw[dt_water$nm==wv[i]]+bbp[,i])
  
  ag443 = (a[,wv412]-Zeta*a[,wv443])/(Xi-Zeta) - 
    (dt_water$aw[dt_water$nm==412]-Zeta*dt_water$aw[dt_water$nm==443])/(Xi-Zeta)
  
  adg   = rrs
  for(i in 1:length(wv))
    adg[, i] = ag443 * exp(-S * (wv[i]-443))
  
  aph   = rrs
  for(i in 1:length(wv))
    aph[, i] = a[, i] - adg[, i] - dt_water$aw[dt_water$nm==wv[i]]
  
  # Limits
  Upper_lim_for_Rrs667 = 20.0 * (Rrs[, wv560]) ^ 1.5
  Lower_lim_for_Rrs667 = 0.9  * (Rrs[, wv560]) ^ 1.7
  Rrs667_revised = 1.27 * (Rrs[,wv560]) ^ 1.47 +
    0.00018 * (Rrs[,wv490] / Rrs[,wv560])^3.19
  
  if(is.null(wv667)){
    Rrs667_final = Rrs667_revised
  }else{
    Rrs667_final = Rrs[,wv667] 
    w = which(Rrs[,wv667] < Lower_lim_for_Rrs667 | Rrs[,wv667] > Upper_lim_for_Rrs667)
    Rrs667_final[w] = Rrs667_revised[w]
  }
  
  # Chl estimation by Bricaud et al. (1995)
  Bricaud_A = 0.0497 # coefficient from Brewin et al. (2015)
  Bricaud_B = 0.7575 # coefficient from Brewin et al. (2015)
  Chl_Bricaud = (aph[,wv443]/Bricaud_A)^(1/Bricaud_B)
  
  # Chl estimation by Liu et al. (2020)
  aChla665_star = 0.017
  Chl_Liu = aph[,wv667] / aChla665_star
  
  # return results
  result = list()
  result$input  = list(wv    = wv,
                       Rrs   = Rrs,
                       wv412 = wv412,
                       wv443 = wv443,
                       wv490 = wv490,
                       wv560 = wv560,
                       wv667 = wv667)
  result$const  = list(g0    = g0,
                       g1    = g1,
                       lam0  = lam0)
  result$Var    = list(rrs   = rrs,
                       u     = u,
                       Chi   = Chi,
                       Yita  = Yita,
                       Zeta  = Zeta,
                       S     = S,
                       Xi    = Xi)
  result$IOP    = list(a0    = a0,
                       bbp0  = bbp0,
                       ag443 = ag443,
                       a     = a,
                       adg   = adg,
                       aph   = aph,
                       bbp   = bbp)
  
  result$Chl    = list(Chl_Bricaud          = Chl_Bricaud,
                       Bricaud_A            = Bricaud_A,
                       Bricaud_B            = Bricaud_B,
                       Chl_Liu              = Chl_Liu,
                       aChla665_star        = 0.017)
  result$residu = list(Upper_lim_for_Rrs667 = Upper_lim_for_Rrs667,
                       Lower_lim_for_Rrs667 = Lower_lim_for_Rrs667,
                       Rrs667_revised       = Rrs667_revised,
                       Rrs667_final         = Rrs667_final)
  
  return(result)
  
}



#' @name run_all_Chla_algorithms
#' @title Run all Chla algorithms
#' @param Rrs Dataframe that should with required colnames 443, 490, 510, 560, 620, 665, 681, 709, 754, 779
#' @param wv_range Number that used to define the range of wavelength to capture
#'   the center wavelength of required band
#' @details Run \link{Chla_algorithms_name} to see all Chla algorithms
#' @export
#' @return A list including Chla results by supported algorithms if Rrs includes the required wavelength.
#'   Retrun \code{NA} values for algorithms of which required wavelength cannot be found.
#' @family Algorithms: Chla concentration
#' 
#' @examples 
#' library(FCMm)
#' data(WaterSpec35)
#' out = run_all_Chla_algorithms(WaterSpec35[,-c(1,2)])
#' 
#' @importFrom stringr str_subset
run_all_Chla_algorithms <- function(Rrs, wv_range=3){
  wv <- str_subset(names(Rrs), '\\d') %>% as.numeric # to be updated from `plot_spec_from_df`
  wv_need <- c(443, 490, 510, 560, 620, 665, 681, 709, 754, 779)
  RrsNA <- rep(NA, nrow(Rrs))
  for(i in wv_need){
    w = which(abs(wv-i) <= wv_range)
    if(length(w) == 0){
      text <- sprintf('Rrs%s <- RrsNA', i) # if cannot find this wavelength, give it NA values
      eval(parse(text=text))
    }else{
      text <- sprintf('Rrs%s <- Rrs[,%s]', i, w)
      eval(parse(text=text))
    }
  }
  result <- list(
    BR_Gil10 = BR_Gil10(Rrs665=Rrs665, Rrs709=Rrs709)$Chla,
    TBA_Gil10 = TBA_Gil10(Rrs665, Rrs709, Rrs754)$Chla,
    C6 = C6(Rrs665, Rrs754)$Chla,
    Bloom = Bloom(Rrs665, Rrs754)$Chla,
    OC4_OLCI = OC4_OLCI(Rrs443, Rrs490, Rrs510, Rrs560)$Chla,
    OCI_Hu12 = OCI_Hu12(Rrs443, Rrs490, Rrs510, Rrs560, Rrs665)$Chl_OCI,
    BR_Git11 = BR_Git11(Rrs665, Rrs709)$Chla,
    TBA_Git11 = TBA_Git11(Rrs665, Rrs709, Rrs754)$Chla,
    NDCI_Mi12 = NDCI_Mi12(Rrs665, Rrs709)$Chla,
    FBA_Le13 = FBA_Le13(Rrs665, Rrs681, Rrs709)$Chla,
    FBA_Yang10 = FBA_Yang10(Rrs665, Rrs709, Rrs754)$Chla,
    SCI_Shen10 = SCI_Shen10(Rrs560, Rrs620, Rrs665, Rrs681)$Chla_EPCHC,
    Gons08 = Gons08(Rrs665, Rrs709, Rrs779)$Chla,
    TC2_final = TC2(Rrs443, Rrs560, Rrs665, Rrs709, Rrs754)$Chla_final,
    TC2_clean = TC2(Rrs443, Rrs560, Rrs665, Rrs709, Rrs754)$Chla_clean,
    TC2_turbid = TC2(Rrs443, Rrs560, Rrs665, Rrs709, Rrs754)$Chla_turbid)
  w_zero <- lapply(result, length) %>% as.matrix %>% {which(. == 0)}
  for(i in w_zero){
    result[[i]] <- NA
  }
  return(result)
}




#' @export
#' @rdname run_all_Chla_algorithms
#' 
Chla_algorithms_name <- function(){
  
  message('Please run QAA_v5 alone if required.')
  return(c('BR_Gil10', 'BR_Git11',
           'TBA_Gil10', 'TBA_Git11',
           'C6', 'Bloom',
           'OC4_OLCI', 'OC5_OCLI', 'OC6_OLCI',
           'OCI_Hu12',
           'NDCI_Mi12', 'Gons08',
           'FBA_Le13', 'FBA_Yang10',
           'SCI_Shen10',
           'TC2', 'TC2_turbid', 'TC2_clean',
           'QAA_v5'))
}


