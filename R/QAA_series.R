#' @title TC2
#' @param Rrs443 Rrs443
#' @param Rrs560 Rrs560
#' @param Rrs665 Rrs665
#' @param Rrs709 Rrs709
#' @param Rrs754 Rrs754
#' @export
#' @return A list with names as Chla_final, Chla_clean, Chla_turbid, and flag
#' @family Algorithms: Chla concentration
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

#' @title TC2_clean
#' @param Rrs443 Rrs443
#' @param Rrs560 Rrs560
#' @param Rrs665 Rrs665
#' @param Rrs709 Rrs709
#' @export
#' @family Algorithms: Chla concentration
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


#' @title TC2_turbid
#' @param Rrs443 Rrs443
#' @param Rrs560 Rrs560
#' @param Rrs665 Rrs665
#' @param Rrs709 Rrs709
#' @param Rrs754 Rrs754
#' @export
#' @family Algorithms: Chla concentration 
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

#' @title QAA_v5
#' @param wv wv
#' @param Rrs Rrs
#' @param wv412 Wavelength index of 412 nm
#' @param wv443 Wavelength index of 443 nm
#' @param wv490 Wavelength index of 490 nm
#' @param wv560 Wavelength index of 560 nm
#' @param wv667 Wavelength index of 667 nm
#' @param verbose verbos (default as FASLE)
#' @export
#' @family Algorithms: Chla concentration 
#' @references Lee Z P, Carder K L, Arnone R A. Deriving inherent optical properties from water color: 
#'   a multiband quasi-analytical algorithm for optically deep waters[J]. Applied optics, 2002, 41(27): 
#'   5755-5772.
QAA_v5 <- function(wv, Rrs,
                   wv412=NULL, wv443=NULL, wv490=NULL, wv560=NULL, wv667=NULL,
                   verbose=F){
  
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
  
  # save results
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