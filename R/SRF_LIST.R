#' @title Collection of spectral responding functions
#' @description A list inlcuding spectral responding functions of several popular sensors
#' @docType data
#' @keywords datasets
#' @name SRF_LIST
#' @usage SRF_LIST
#' @format List with four elements: sensor (character), srf (data.frame), cw_med (num), and cw_max (num)
#' @details The cw_med and cw_max are different sicne 
#'   we use two functions \code{find_center_wavelength_med} and 
#'   \code{find_center_wavelength_max} to obtain thier center wavelength.
#'   \itemize{
#'     \item \code{find_center_wavelength_med} center wavelength is determined by the wavelength
#'       position of half-maximum width
#'     \item \code{find_center_wavelength_max} wavelength at the position of maximum srf
#'   }
#'   All included sensors are Sentinel-2A, Landsat-8, GF1-WFV1, GF1-WFV2, GF1-WFV3, GF1-WFV4, 
#'     GF2-PMS1, GF2-PMS2, GF4-PMI, HJ1A-CCD1, HJ1A-CCD2, HJ1B-CCD1, HJ1B-CCD2, MODISA, MERIS, 
#'     VIIRS, OLCI, GOCI, and GF5-VIMS.
#'   The data are collected from their corresponding official websits (see in references).
#'   
#' @references 
#'   \itemize{
#'     \item http://www.cresda.com/CN/
#'     \item https://oceancolor.gsfc.nasa.gov
#'     \item https://sentinels.copernicus.eu/web/sentinel/home
#'   }
#' @family Datasets
NULL
