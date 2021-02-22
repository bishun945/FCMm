#Blend algorithms####


#' @name Blend_Smith18
#' @title Blending algorithm by Smith et al. (2018)
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
#' res = Blend_Smith18(WaterSpec35$`442.5`, WaterSpec35$`490`, WaterSpec35$`510`, 
#' WaterSpec35$`560`, WaterSpec35$`665`, WaterSpec35$`708.75`)
#' @references Smith M E, Lain L R, Bernard S. An optimized chlorophyll a switching 
#'   algorithm for MERIS and OLCI in phytoplankton-dominated waters[J]. 
#'   Remote Sensing of Environment, 2018, 215: 217-227.
Blend_Smith18 <- function(Rrs443, Rrs490, Rrs510, Rrs560, Rrs665, Rrs709) {
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
    Chla_blend = CHL_BLEND,
    PHI = PHI,
    a1  = a1,
    a2  = a2
  )
  return(result)
}

#' @name Blend_Jac17
#' @title Algorithm blending framework by Jackson et al. (2017)
#' @param Rrs Rrs data.frame input for matching centroids and calculating Chla
#' @param wv_range Number that used to define the range of wavelength to capture
#'   the center wavelength of required band
#' @param ... parameters of \link{apply_FCM_m}
#' @return A list includes \code{Chla_blend}, \code{Rrs}, \code{res_FCM} (the return of \link{apply_FCM_m}),
#'   and estimates from optimal algorithms.
#' @family Algorithms: Chla concentration
#' @examples 
#' library(FCMm)
#' data(WaterSpec35)
#' res = Blend_Jac17(WaterSpec35[, -c(1, 2)])
#' @references 
#' 
#' Jackson T, Sathyendranath S, Mélin F. An improved optical classification scheme for the Ocean 
#'   Colour Essential Climate Variable and its applications[J]. Remote Sensing of Environment, 2017, 
#'   203: 152-161.
#'   
#' Moore T S, Dowell M D, Bradt S, et al. An optical water type framework for 
#'   selecting and blending retrievals from bio-optical algorithms in lakes and coastal waters[J]. 
#'   Remote sensing of environment, 2014, 143: 97-111.
#' 
#' @export
#' 
Blend_Jac17 <- function(Rrs, wv_range = 5, ...) {
  
  # load Jackson centroids
  Jac17_cen <- system.file("Jac17_centroids_Rrs.csv", package = "FCMm") %>%
    read.csv(., comment.char = "#") %>%
    setNames(., gsub("X", "", colnames(.))) %>%
    magrittr::set_rownames(., .[,1]) %>%
    .[,-1]
  
  # select the required bands in the param Rrs
  wv <- stringr::str_subset(names(Rrs), '\\d') %>% as.numeric()
  wv_need <- colnames(Jac17_cen) %>% as.numeric()
  Rrs_need <- matrix(NA, ncol = ncol(Jac17_cen), nrow = nrow(Rrs)) %>% 
    as.data.frame() %>% setNames(., wv_need)
  
  for(i in 1:ncol(Rrs_need)){
    w = which.min(abs(wv-wv_need[i]))
    if(length(w) == 0){
      stop(sprintf("%s nm can not found from the input data.frame!", wv_need[i]))
    }else if(length(w) > 1){
      stop(paste0("Bands ", 
                  paste0(wv_need[w], collapse = ", "),
                  " duplicate in the input data.frame!"))
    }else {
      
      Rrs_need[, i] <- Rrs[, as.character(wv)][, w]
      
    }
  }
  
  if(!exists("m_used")) {
    
    res_FCM <- apply_FCM_m(
      Rrs = Rrs_need,
      wavelength = as.numeric(colnames(Rrs_need)),
      Rrs_clusters = Jac17_cen,
      do.stand = TRUE,
      default.cluster = FALSE,
      m_used = 1.5,
      ...
    )
    
  }else {
    
    res_FCM <- apply_FCM_m(
      Rrs = Rrs_need,
      wavelength = as.numeric(colnames(Rrs_need)),
      Rrs_clusters = Jac17_cen,
      do.stand = TRUE,
      default.cluster = FALSE,
      ...
    )
    
  }
  
  # OWT1-7 OCI
  Chl_OCI <- OCI_Hu12(
    Rrs443 = Rrs_need$`443`,
    Rrs490 = Rrs_need$`490`,
    Rrs510 = Rrs_need$`510`,
    Rrs560 = Rrs_need$`555`,
    Rrs665 = Rrs_need$`670`
  )$Chl_OCI
  
  # OWT13-14 OC3
  Chl_OC3 <- OC3_OLCI(
    Rrs443 = Rrs_need$`443`,
    Rrs490 = Rrs_need$`490`,
    Rrs560 = Rrs_need$`555`
  )$Chla
  
  # OWT8-12 OC5
  Chl_OC5 <- OC5_OLCI(
    Rrs412 = Rrs_need$`412`,
    Rrs443 = Rrs_need$`443`,
    Rrs490 = Rrs_need$`490`,
    Rrs510 = Rrs_need$`510`,
    Rrs560 = Rrs_need$`555`
  )$Chla
  
  # blend result
  Chl_opt <- matrix(NA, 
                    ncol = ncol(res_FCM$u),
                    nrow = nrow(res_FCM$u))
  
  Chl_opt[, 1:7] <- Chl_OCI
  Chl_opt[, 8:12] <- Chl_OC5
  Chl_opt[, 13:14] <- Chl_OC3
  
  Chl_opt[Chl_opt < 0] <- 0
  
  Chla_blend <- apply(Chl_opt * res_FCM$u, 1, sum)
  
  return(list(
    Chla_blend = Chla_blend,
    Rrs = Rrs_need,
    res_FCM = res_FCM,
    Chl_OC3 = Chl_OC3,
    Chl_OC5 = Chl_OC5,
    Chl_OCI = Chl_OCI
  ))
  
}


#' @name Blend_Moo14
#' @title Algorithm blending framework by Moore et al. (2014)
#' @param Rrs Rrs data.frame input for matching centroids and calculating Chla
#' @param wv_range Number that used to define the range of wavelength to capture
#'   the center wavelength of required band
#' @param ... parameters of \link{apply_FCM_m}
#' @return A list includes \code{Chla_blend}, \code{Rrs}, \code{res_FCM} (the return of \link{apply_FCM_m}),
#'   and estimates from optimal algorithms.
#' @importFrom magrittr %>% set_rownames
#' @importFrom utils read.csv
#' @importFrom stats setNames
#' @importFrom stringr str_subset
#' @family Algorithms: Chla concentration
#' @examples 
#' library(FCMm)
#' data(WaterSpec35)
#' res = Blend_Moo14(WaterSpec35[, -c(1, 2)])
#' @references Moore T S, Dowell M D, Bradt S, et al. An optical water type framework for 
#'   selecting and blending retrievals from bio-optical algorithms in lakes and coastal waters[J]. 
#'   Remote sensing of environment, 2014, 143: 97-111.
#' @export
#' 
Blend_Moo14 <- function(Rrs, wv_range = 3, ...) {
  
  # load Moore centroids
  Moo14_cen <- system.file("Moo14_centroids_Rrs0-.csv", package = "FCMm") %>%
    read.csv(., comment.char = "#") %>%
    setNames(., gsub("X", "", colnames(.))) %>%
    magrittr::set_rownames(., .[,1]) %>%
    .[,-1] %>%
    {0.52 / (1/. - 1.7)} # convert Rrs0- to Rrs0+
  
  # select the required bands in the param Rrs
  wv <- stringr::str_subset(names(Rrs), '\\d') %>% as.numeric()
  wv_need <- colnames(Moo14_cen) %>% as.numeric()
  Rrs_need <- matrix(NA, ncol = ncol(Moo14_cen), nrow = nrow(Rrs)) %>% 
    as.data.frame() %>% setNames(., wv_need)
  
  for(i in 1:ncol(Rrs_need)){
    w = which.min(abs(wv-wv_need[i]))
    if(length(w) == 0){
      stop(sprintf("%s nm can not found from the input data.frame!", wv_need[i]))
    }else if(length(w) > 1){
      stop(paste0("Bands ", 
                  paste0(wv_need[w], collapse = ", "),
                  " duplicate in the input data.frame!"))
    }else {
      
      Rrs_need[, i] <- Rrs[, as.character(wv)][, w]
      
    }
  }
  
  if(!exists("m_used")) {
    
    res_FCM <- apply_FCM_m(
      Rrs = Rrs_need,
      wavelength = as.numeric(colnames(Rrs_need)),
      Rrs_clusters = Moo14_cen,
      do.stand = TRUE,
      default.cluster = FALSE,
      m_used = 1.5,
      ...
    )
    
  }else {
    
    res_FCM <- apply_FCM_m(
      Rrs = Rrs_need,
      wavelength = as.numeric(colnames(Rrs_need)),
      Rrs_clusters = Moo14_cen,
      do.stand = TRUE,
      default.cluster = FALSE,
      ...
    )
    
  }
  
  
  # OWT 1 2 3 6
  Chl_OC4 <- OC4_OLCI(
    Rrs443 = Rrs_need$`443`,
    Rrs490 = Rrs_need$`490`,
    Rrs510 = Rrs_need$`510`,
    Rrs560 = Rrs_need$`560`
  )$Chla
  
  # OWT 4 5 7
  Chl_Git11 <- TBA_Git11(
    Rrs665 = Rrs_need$`665`,
    Rrs709 = Rrs_need$`709`,
    Rrs754 = Rrs_need$`753`
  )$Chla
  
  # blend result
  Chl_opt <- matrix(NA, 
                    ncol = ncol(res_FCM$u),
                    nrow = nrow(res_FCM$u))
  
  Chl_opt[, c(1, 2, 3, 6)] <- Chl_OC4
  Chl_opt[, c(4, 5, 7)] <- Chl_Git11
  
  Chl_opt[Chl_opt < 0] <- 0
  
  Chla_blend <- apply(Chl_opt * res_FCM$u, 1, sum)
  
  return(list(
    Chla_blend = Chla_blend,
    Rrs = Rrs_need,
    res_FCM = res_FCM,
    Chl_Git11 = Chl_Git11,
    Chl_OC4 = Chl_OC4
  ))
  
}



#' @name Blend_Bi21
#' @title Algorithm blending framework by Bi et al. (2021)
#' @param Rrs Rrs data.frame input for matching centroids and calculating Chla
#' @param wv_range Number that used to define the range of wavelength to capture
#'   the center wavelength of required band
#' @param ... parameters of \link{apply_FCM_m}
#' @return A list includes \code{Chla_blend}, \code{Rrs}, \code{res_FCM} (the return of \link{apply_FCM_m}),
#'   and estimates from optimal algorithms.
#' @importFrom magrittr %>% set_rownames
#' @importFrom utils read.csv
#' @importFrom stats setNames
#' @importFrom stringr str_subset
#' @family Algorithms: Chla concentration
#' @examples 
#' library(FCMm)
#' data(WaterSpec35)
#' res = Blend_Bi21(WaterSpec35[, -c(1,2)])
#' @export
#' 
Blend_Bi21 <- function(Rrs, wv_range = 3, ...) {
  
  # load Bi21 centroids
  Bi21_cen <- system.file("Bi21_centroids_Rrs_norm.csv", package = "FCMm") %>%
    read.csv(., comment.char = "#") %>%
    setNames(., gsub("X", "", colnames(.))) %>%
    magrittr::set_rownames(., .[,1]) %>%
    .[,-1]
  
  # select the required bands in the param Rrs
  wv <- stringr::str_subset(names(Rrs), '\\d') %>% as.numeric()
  wv_need <- colnames(Bi21_cen) %>% as.numeric()
  Rrs_need <- matrix(NA, ncol = ncol(Bi21_cen), nrow = nrow(Rrs)) %>% 
    as.data.frame() %>% setNames(., wv_need)
  
  for(i in 1:ncol(Rrs_need)){
    w = which.min(abs(wv-wv_need[i]))
    if(length(w) == 0){
      stop(sprintf("%s nm can not found from the input data.frame!", wv_need[i]))
    }else if(length(w) > 1){
      stop(paste0("Bands ", 
                  paste0(wv_need[w], collapse = ", "),
                  " duplicate in the input data.frame!"))
    }else {
      
      Rrs_need[, i] <- Rrs[, as.character(wv)][, w]
      
    }
  }
  
  if(!exists("m_used")) {
    
    res_FCM <- apply_FCM_m(
      Rrs = Rrs_need,
      wavelength = as.numeric(colnames(Rrs_need)),
      Rrs_clusters = Bi21_cen,
      do.stand = TRUE,
      default.cluster = FALSE,
      m_used = 1.5,
      ...
    )
    
  }else {
    
    res_FCM <- apply_FCM_m(
      Rrs = Rrs_need,
      wavelength = as.numeric(colnames(Rrs_need)),
      Rrs_clusters = Bi21_cen,
      do.stand = TRUE,
      default.cluster = FALSE,
      ...
    )
    
  }
  
  res_Chla <- run_all_Chla_algorithms(Rrs) %>% as.data.frame()
  res_Chla$Bloom <- res_Chla$C6
  res_Chla$C6 <- NULL
  
  source(system.file("LUT_OPT_BI2021.R", package = "FCMm"))
  
  # blending
  res_Chla_opt <- res_Chla[, LUT_OPT$Opt_algorithm]
  memb <- round(res_FCM$u, 5)
  res_blend <- Chla_algorithms_blend2(res_Chla_opt, memb, LUT_OPT$Opt_algorithm, LUT_OPT$Remove_algorithm)
  
  return(list(
    Chla_blend = res_blend$Chla_blend,
    Rrs = Rrs_need,
    res_FCM = res_FCM,
    dt_Chla = res_Chla,
    res_blend = res_blend
  ))
  
}



#' @name Blend_FCMm
#' @title Algorithm blending framework by Bi theses
#' @param Rrs Rrs data.frame input for matching centroids and calculating Chla
#' @param wv_range Number that used to define the range of wavelength to capture
#'   the center wavelength of required band
#' @param reparam Whether to apply reparameterization coefficients
#' @param ... parameters of \link{apply_FCM_m}
#' @return A list includes \code{Chla_blend}, \code{Rrs}, \code{res_FCM} (the return of \link{apply_FCM_m}),
#'   and estimates from optimal algorithms.
#' @importFrom magrittr %>% set_rownames
#' @importFrom utils read.csv
#' @importFrom stats setNames
#' @importFrom stringr str_subset
#' @family Algorithms: Chla concentration
#' @examples 
#' library(FCMm)
#' data(WaterSpec35)
#' # res = Blend_FCMm(WaterSpec35[, -c(1,2)])
#' @export
#' 
Blend_FCMm <- function(Rrs, wv_range = 3, reparam = TRUE, ...) {
  
  # load BIPHD centroids
  BIPHD_cen <- system.file("BiPHD_centroids_Rrs_norm.csv", package = "FCMm") %>%
    read.csv(., comment.char = "#") %>%
    setNames(., gsub("X", "", colnames(.))) %>%
    magrittr::set_rownames(., .[,1]) %>%
    .[,-1]
  
  # select the required bands in the param Rrs
  wv <- stringr::str_subset(names(Rrs), '\\d') %>% as.numeric()
  wv_need <- colnames(BIPHD_cen) %>% as.numeric()
  Rrs_need <- matrix(NA, ncol = ncol(BIPHD_cen), nrow = nrow(Rrs)) %>% 
    as.data.frame() %>% setNames(., wv_need)
  
  for(i in 1:ncol(Rrs_need)){
    w = which.min(abs(wv-wv_need[i]))
    if(length(w) == 0){
      stop(sprintf("%s nm can not found from the input data.frame!", wv_need[i]))
    }else if(length(w) > 1){
      stop(paste0("Bands ", 
                  paste0(wv_need[w], collapse = ", "),
                  " duplicate in the input data.frame!"))
    }else {
      
      Rrs_need[, i] <- Rrs[, as.character(wv)][, w]
      
    }
  }
  
  if(!exists("m_used")) {
    
    res_FCM <- apply_FCM_m(
      Rrs = Rrs_need,
      wavelength = as.numeric(colnames(Rrs_need)),
      Rrs_clusters = BIPHD_cen,
      do.stand = TRUE,
      default.cluster = FALSE,
      m_used = 1.5,
      ...
    )
    
  }else {
    
    res_FCM <- apply_FCM_m(
      Rrs = Rrs_need,
      wavelength = as.numeric(colnames(Rrs_need)),
      Rrs_clusters = BIPHD_cen,
      do.stand = TRUE,
      default.cluster = FALSE,
      ...
    )
    
  }
  
  res_Chla <- run_all_Chla_algorithms(Rrs) %>% as.data.frame()
  res_Chla$Bloom <- res_Chla$C6
  res_Chla$C6 <- NULL
  
  source(system.file("LUT_OPT_BIPHD.R", package = "FCMm"))
  
  # blending
  res_Chla_opt <- res_Chla[, LUT_OPT$Opt_algorithm]
  memb <- round(res_FCM$u, 5)
  res_blend <- Chla_algorithms_blend2(res_Chla_opt, memb, LUT_OPT$Opt_algorithm, LUT_OPT$Remove_algorithm)
  
  if(reparam) {
    
    blend_coef <- read.csv(system.file("blend_coef_BIPHD.csv", package = "FCMm"))
    
    tmp <- data.frame(
      lgP  = log10(res_blend$Chla_blend),
      clus = res_FCM$cluster
    )
    tmp$b <- tmp$a <- NA
    
    for(i in 1:nrow(blend_coef)) {
      w_cof <- tmp$clus == blend_coef$OWT[i]
      tmp[w_cof, c("a", "b")] <- blend_coef[i, c("a", "b")]
    }
    
    Chla_reparam <- 10^(tmp$a * tmp$lgP + tmp$b)
    
  } else {
    
    Chla_reparam <- NULL
    
  }
  
  
  return(list(
    Chla_blend = res_blend$Chla_blend,
    Rrs = Rrs_need,
    res_FCM = res_FCM,
    dt_Chla = res_Chla,
    res_blend = res_blend,
    Chla_reparam = Chla_reparam
  ))
  
}

