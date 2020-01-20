#' @title Simulation of hyperspectral data by SRF
#' @name SRF_simulate
#' @description 
#' Simulate hyperspectral Rrs to multispectral bands
#'   via sensors SRF (Spectral response function).
#' @param Rrs A data.frame with colnames as Wavelength + SampleName, the first column is wavelength.
#' @param select_sensor Character. Select sensors. Use \code{show_sensor_names()} to print
#'   all supported sensors. Default as \code{All}
#' @param output_wavelength Character. \code{MED} (default) or \code{MAX}.
#'   Define the center wavelength. \code{MED} means the center wavelength is 
#'   middle position of half maximun of max peak. While \code{MAX} means the 
#'   poition at the maximun of SRF peak.
#' @param save_as_csv Logical. Choose to save the simulation results as single csv for each
#'   sensor. Default with \code{FALSE}
#' @param na.rm Logical. Should NA values be removed? Default as \code{TRUE}
#' @param wv_as_column Logical. If \code{TRUE} (default), the output result is a dataframe
#'   with wavelength as column names.
#' 
#' @export
#' @return 
#' A \code{list} with names as all selected sensors from parameters \code{select_sensor}.
#' For each list, including five elements:
#'   \itemize{
#'     \item \strong{sensor} Sensor name
#'     \item \strong{srf} Spectral response function of the sensor
#'     \item \strong{cw_med} Center wavelength by method \code{MED}
#'     \item \strong{cw_max} Center wavelength by method \code{MAX}
#'     \item \strong{Rrs_simu} The simulation of Rrs by supported SRF
#'   }
#' 
#' @examples 
#' \dontrun{
#' library(FCMm)
#' nm <- seq(400,900)
#' Rrs <- data.frame(nm=nm,Site1=exp(nm/1000)+runif(501))
#' result <- SRF_simulate(Rrs,select_sensor=c("OLCI","MODIS"))
#' }
#' 
#' @importFrom stats setNames
#' @family Utils

SRF_simulate <- function(Rrs,
                         select_sensor="All",
                         output_wavelength="MED",
                         save_as_csv=FALSE,
                         na.rm=TRUE,
                         wv_as_column=TRUE){
  
  if(select_sensor[1] == "All" & length(select_sensor) == 1){
    sensors <- show_sensor_names()
  }else{
    sensors <- select_sensor
    for(i in 1:length(sensors))
      sensors[i] <- match.arg(sensors[i], show_sensor_names())
  }
  print(sensors)
  result <- list()
  for(sensor in sensors){
    
    # load SRF
    SRF <- SRF_LIST[[sensor]]
    
    # creat Rrs_simu
    Rrs_simu <- matrix(nrow=length(SRF$cw_med),
                       ncol=ncol(Rrs)-1,
                       data=0) %>% as.data.frame
    names(Rrs_simu) <- names(Rrs)[-1]
    if(output_wavelength == "MED"){
      Rrs_simu <- cbind(nm=SRF$cw_med, Rrs_simu)
    }else{
      Rrs_simu <- cbind(nm=SRF$cw_max, Rrs_simu)
    }
    
    # subset the Rrs dataframe and SRF dataframe by wavelength intersect
    w <- intersect(Rrs[,1], SRF$srf[,1])
    
    # do simulation via SRF
    for(i in 2:ncol(Rrs_simu)){
      Rrs_simu[,i] <- cal_SRF(Rrs_single = Rrs[(Rrs[,1] %in% w),i],
                              srf = SRF$srf[(SRF$srf[,1] %in% w),])
    }
    if(na.rm==TRUE){
      SRF$Rrs_simu <- Rrs_simu %>% na.omit
    }else{
      SRF$Rrs_simu <- Rrs_simu
    }
    
    if(wv_as_column){
      SRF$Rrs_simu <- t(SRF$Rrs_simu)[-1,] %>% as.data.frame %>% setNames(., SRF$Rrs_simu[,1])
    }
    
    result[[sensor]] <- SRF
  }
  if(save_as_csv){
    for(sheet in names(result)){
      write.csv(result[[sheet]]$Rrs_simu,
                file=sprintf('./Simulated_Rrs_%s.csv',sheet),
                row.names=F)
    }
  }
  return(result)
}

cal_SRF <- function(Rrs_single, srf){
  result <- matrix(nrow=ncol(srf)-1, data=0) %>% as.data.frame()
  for(i in 2:(ncol(srf))){
    result[i-1,] <- sum(Rrs_single * srf[,i], na.rm=T) / sum(srf[,i], na.rm=T)
  }
  return(result)
}

read_srf_excel <- function(fn){
  SRF_LIST <- list()
  sheets <- excel_sheets(fn)
  for(sheet in sheets){
    dt <- read_excel(fn, sheet=sheet) %>% as.data.frame
    center_wavelength_med <- find_center_wavelength_med(dt)
    center_wavelength_max <- find_center_wavelength_max(dt)
    tmp <- list(sensor=sheet, 
                srf=dt,
                cw_med=center_wavelength_med,
                cw_max=center_wavelength_max)
    SRF_LIST[[sheet]] <- tmp
  }
  return(SRF_LIST)
}

find_center_wavelength_med <- function(dt){
  nm <- dt[,1]
  center_wv <- 2:dim(dt)[2]
  j <- 1
  for(i in 2:dim(dt)[2]){
    tmp <- dt[,i]
    center_wv[j] <- which(tmp >= max(tmp)/2) %>% 
      nm[.] %>% range %>% mean %>% round
    j <- j + 1
  }
  return(center_wv)
}

find_center_wavelength_max <- function(dt){
  nm <- dt[,1]
  center_wv <- 2:dim(dt)[2]
  j <- 1
  for(i in 2:dim(dt)[2]){
    tmp <- dt[,i]
    center_wv[j] <- which.max(tmp) %>% 
      nm[.]
    j <- j + 1
  }
  return(center_wv)
}

#' @title show_sensor_names
#' @name show_sensor_names
#' @description show_sensor_names
#' @usage show_sensor_names()
#' @export

show_sensor_names <- function(){
  return(names(SRF_LIST))
}
