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
#'   (2020-02-28) The \code{C6} model was replaced by \code{Bloom} model.
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

#' @title BR_Gil10
#' @param Rrs665 Rrs665
#' @param Rrs709 Rrs709
#' @export
#' @family Algorithms: Chla concentration
#' @references Gilerson A A, Gitelson A A, Zhou J, et al. Algorithms for remote estimation of 
#'   chlorophyll-a in coastal and inland waters using red and near infrared bands[J]. 
#'   Optics Express, 2010, 18(23): 24109-24125.
BR_Gil10 <- function(Rrs665, Rrs709){
  ind = (Rrs709/Rrs665)
  Chla = (35.75*ind-19.3)^1.124
  result <- list(ind=ind,
                 Chla=Chla)
  return(result)
}

#' @title TBA_Gil10
#' @param Rrs665 Rrs665
#' @param Rrs709 Rrs709
#' @param Rrs754 Rrs754
#' @export
#' @family Algorithms: Chla concentration
#' @references Gilerson A A, Gitelson A A, Zhou J, et al. Algorithms for remote estimation of 
#'   chlorophyll-a in coastal and inland waters using red and near infrared bands[J]. 
#'   Optics Express, 2010, 18(23): 24109-24125.
TBA_Gil10 <- function(Rrs665, Rrs709, Rrs754){
  ind = (1/Rrs665-1/Rrs709)*Rrs754
  Chla = (113.36*ind+16.45)^1.124
  result <- list(ind=ind,
                 Chla=Chla)
  return(result)
}

#' @title C6 (Being deprecated)
#' @param Rrs665 Rrs665
#' @param Rrs754 Rrs754
#' @export
#' @family Algorithms: Chla concentration
C6 <- function(Rrs665, Rrs754){
  ind <- 1/Rrs665*Rrs754
  Chla <- 10^( ind * 0.14 + 2.11)
  result <- list(ind=ind,
                 Chla=Chla)
  return(result)
}

#' @title Bloom
#' @param Rrs665 Rrs665
#' @param Rrs754 Rrs754
#' @export
#' @family Algorithms: Chla concentration
Bloom <- function(Rrs665, Rrs754){
  ind <- 1/Rrs665*Rrs754
  Chla <- 257.740*exp(0.279*ind)
  result <- list(ind=ind,
                 Chla=Chla)
  return(result)
}

#' @title OC4_OLCI
#' @param Rrs443 Rrs443
#' @param Rrs490 Rrs490
#' @param Rrs510 Rrs510
#' @param Rrs560 Rrs560
#' @export
#' @family Algorithms: Chla concentration
#' @references O'Reilly J E, Werdell P J. Chlorophyll algorithms for ocean color
#'   sensors-OC4, OC5 & OC6[J]. Remote sensing of environment, 2019, 229: 32-47.
OC4_OLCI <- function(Rrs443, Rrs490, Rrs510, Rrs560){
  X <- apply(cbind(Rrs443,Rrs490,Rrs510),1,max)/Rrs560
  X <- log10(X)
  Chla <- 10^(0.42540-3.21679*X+2.86907*X^2-0.62628*X^3-1.09333*X^4)
  result <- list(X=X,
                 Chla=Chla)
  return(result)
}

#' @title OC5_OLCI
#' @param Rrs412 Rrs412
#' @param Rrs443 Rrs443
#' @param Rrs490 Rrs490
#' @param Rrs510 Rrs510
#' @param Rrs560 Rrs560
#' @export
#' @family Algorithms: Chla concentration
#' @references O'Reilly J E, Werdell P J. Chlorophyll algorithms for ocean color
#'   sensors-OC4, OC5 & OC6[J]. Remote sensing of environment, 2019, 229: 32-47.
OC5_OLCI <- function(Rrs412, Rrs443, Rrs490, Rrs510, Rrs560){
  X <- apply(cbind(Rrs412,Rrs443,Rrs490,Rrs510),1,max)/Rrs560
  X <- log10(X)
  Chla <- 10^(0.43213-3.13001*X+3.05479*X^2-1.45176*X^3-0.24947*X^4) 
  result <- list(X=X,
                 Chla=Chla)
  return(result)
}


#' @title OC6_OLCI
#' @param Rrs412 Rrs412
#' @param Rrs443 Rrs443
#' @param Rrs490 Rrs490
#' @param Rrs510 Rrs510
#' @param Rrs560 Rrs560
#' @param Rrs665 Rrs665
#' @export
#' @family Algorithms: Chla concentration
#' @references O'Reilly J E, Werdell P J. Chlorophyll algorithms for ocean color
#'   sensors-OC4, OC5 & OC6[J]. Remote sensing of environment, 2019, 229: 32-47.
OC6_OLCI <- function(Rrs412, Rrs443, Rrs490, Rrs510, Rrs560, Rrs665){
  X <- apply(cbind(Rrs412,Rrs443,Rrs490,Rrs510),1,max)/apply(cbind(Rrs560,Rrs665),1,mean)
  X <- log10(X)
  Chla <- 10^(0.95039-3.05404*X+2.17992*X^2-1.12097*X^3-0.15262*X^4)
  result <- list(X=X,
                 Chla=Chla)
  return(result)
}


#' @title OCI_Hu12
#' @param Rrs443 Rrs443
#' @param Rrs490 Rrs490
#' @param Rrs510 Rrs510
#' @param Rrs560 Rrs560
#' @param Rrs665 Rrs665
#' @export
#' @family Algorithms: Chla concentration
#' @references Hu C, Lee Z, Franz B. Chlorophyll aalgorithms for oligotrophic oceans: 
#'   A novel approach based on three‐band reflectance difference[J]. Journal of Geophysical 
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

#' @title BR_Git11
#' @param Rrs665 Rrs665
#' @param Rrs709 Rrs709
#' @export
#' @family Algorithms: Chla concentration
#' @references Gitelson A A, Gurlin D, Moses W J, et al. Remote Estimation of 
#'   Chlorophyll-a Concentration[J]. Advances in environmental remote sensing: sensors, 
#'   algorithms, and applications, 2011: 439.
BR_Git11 <- function(Rrs665, Rrs709){
  ind = Rrs709/Rrs665
  Chla = 72.66*ind-46.535
  result <- list(ind=ind,
                 Chla=Chla)
  return(result)
}

#' @title TBA_Git11
#' @param Rrs665 Rrs665
#' @param Rrs709 Rrs709
#' @param Rrs754 Rrs754
#' @export
#' @family Algorithms: Chla concentration
#' @references Gitelson A A, Gurlin D, Moses W J, et al. Remote Estimation of 
#'   Chlorophyll-a Concentration[J]. Advances in environmental remote sensing: sensors, 
#'   algorithms, and applications, 2011: 439.
TBA_Git11 <- function(Rrs665, Rrs709, Rrs754){
  ind = (1/Rrs665-1/Rrs709)*Rrs754
  Chla = 243.86*ind+23.17
  result <- list(ind=ind,
                 Chla=Chla)
  return(result)
}

#' @title NDCI_Mi12
#' @param Rrs665 Rrs665
#' @param Rrs709 Rrs709
#' @export
#' @family Algorithms: Chla concentration
#' @note Chla_test as a test verison of NDCI_Mi12 results parameterized by Bi
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

#' @title FBA_Le13
#' @param Rrs665 Rrs665
#' @param Rrs681 Rrs681
#' @param Rrs709 Rrs709
#' @export
#' @family Algorithms: Chla concentration
#' @references Le C, Hu C, Cannizzaro J, et al. Evaluation of chlorophyll-a remote
#'   sensing algorithms for an optically complex estuary[J]. Remote Sensing of Environment, 
#'   2013, 129: 75-89.
FBA_Le13 <- function(Rrs665, Rrs681, Rrs709){
  ind = (1/Rrs665-1/Rrs681)*(1/Rrs709-1/Rrs681)
  Chla = 18.492*ind+6.1302
  result <- list(ind=ind,
                 Chla=Chla)
  return(result)
}

#' @title FBA_Yang10
#' @param Rrs665 Rrs665
#' @param Rrs709 Rrs709
#' @param Rrs754 Rrs754
#' @export
#' @family Algorithms: Chla concentration
#' @references Yang W, Matsushita B, Chen J, et al. An enhanced three-band index for estimating 
#'   chlorophyll-a in turbid case-II waters: case studies of Lake Kasumigaura, Japan, and Lake 
#'   Dianchi, China[J]. IEEE Geoscience and Remote Sensing Letters, 2010, 7(4): 655-659.
FBA_Yang10 <- function(Rrs665, Rrs709, Rrs754){
  ind = (1/Rrs665-1/Rrs709)/(1/Rrs754-1/Rrs709)
  Chla = 161.24*ind+28.04
  result <- list(ind=ind,
                 Chla=Chla)
  return(result)
}

#' @title SCI_Shen10
#' @param Rrs560 Rrs560
#' @param Rrs620 Rrs620
#' @param Rrs665 Rrs665
#' @param Rrs681 Rrs681
#' @export
#' @family Algorithms: Chla concentration
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

#' @title Gons08
#' @param Rrs665 Rrs665
#' @param Rrs709 Rrs709
#' @param Rrs779 Rrs779
#' @export
#' @family Algorithms: Chla concentration
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

#' @title Chla_algorithms_name
#' @export
#' @family Algorithms: Chla concentration
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

#' @title run_all_Chla_algorithms
#' @param Rrs Dataframe that should with required colnames 443, 490, 510, 560, 620, 665, 681, 709, 754, 779
#' @param wv_range Number that used to define the range of wavelength to capture
#'   the center wavelength of required band
#' @details Please type \code{Chla_algorithm_name()} to see all Chla algorithms
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


