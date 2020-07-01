#' @name Rrs_trans
#' @title To transform the Rrs values from one wavelength to another 
#'   based on the predefined relationship
#' 
#' @param Rrs490 Rrs490
#' @param Rrs560 Rrs560
#' @param Rrs674 Rrs674
#' @param Rrs754 Rrs754
#' @param Rrs865 Rrs865
#' @return A list including transformed Rrs value: Rrs4865, Rrs560, 
#'   Rrs674, Rrs754, and Rrs862
#' 
#' @note Maybe I will create a lookuptable for transformation across the whole range.
#' @details The trans list is:
#'   (OLCI) -> (VIIRS)
#'   Rrs490 -> Rrs486
#'   Rrs560 -> Rrs551
#'   Rrs674 -> Rrs671
#'   Rrs754 -> Rrs745
#'   Rrs865 -> Rrs862
#'
#' @examples Rrs_new = Rrs_trans(WaterSpec35$`490`, WaterSpec35$`560`, 
#' WaterSpec35$`673.75`, WaterSpec35$`753.75`, WaterSpec35$`865`) 
#' 
#' @noRd
#' 
Rrs_trans <- function(Rrs490, Rrs560, Rrs674, Rrs754, Rrs865){
  
  stopifnot({
    is.numeric(Rrs490)
    is.numeric(Rrs560)
    is.numeric(Rrs674)
    is.numeric(Rrs754)
    is.numeric(Rrs865)
  })
  
  return(list(
    Rrs486 = 0.985 * Rrs490,
    Rrs551 = 0.975 * Rrs560,
    Rrs671 = 0.999 * Rrs674,
    Rrs745 = 0.978 * Rrs754,
    Rrs862 = 1.001 * Rrs865
  ))
}


#' @name GAA_SPM
#' @title Globally applicable algorithm to seamlessly retrieve the 
#'   concentration of suspended particulate matter from remote sensing reflectance.
#' @param Rrs490 Rrs490
#' @param Rrs560 Rrs560
#' @param Rrs674 Rrs674
#' @param Rrs754 Rrs754
#' @param Rrs865 Rrs865
#' @return The result returned by \code{GAA_SPM} is a list including:
#' \itemize{
#'   \item \strong{param} Coefficients of C0, C1, C2, C3, a1, a2.
#'   \item \strong{W} The weighting factors (data.frame).
#'   \item \strong{GI} Generalized index for SPM.
#'   \item \strong{SPM} The predicted suspended particulate matter.
#' }
#' @noRd
#' 
#' @note 
#' The input wavelength are at OLCI bandset but this function will
#'   transform these wavelengthes to the GAA required bands.
#' 
#' 2020-06-27
#' 
#' The TSM functions are under test and will be exported in the next few versions.
#' 
#' @examples res = GAA_SPM(WaterSpec35$`490`, WaterSpec35$`560`, 
#' WaterSpec35$`673.75`, WaterSpec35$`753.75`, WaterSpec35$`865`) 
#' 
#' @references 
#' Yu X, Lee Z, Shen F, et al. An empirical algorithm to seamlessly retrieve the 
#'   concentration of suspended particulate matter from water color across ocean to 
#'   turbid river mouths[J]. Remote Sensing of Environment, 2019, 235: 111491.
#'   
GAA_SPM <- function(Rrs490, Rrs560, Rrs674, Rrs754, Rrs865){
  
  Rrs_ = Rrs_trans(Rrs490, Rrs560, Rrs674, Rrs754, Rrs865)
  Rrs486 = Rrs_$Rrs486
  Rrs551 = Rrs_$Rrs551
  Rrs671 = Rrs_$Rrs671 # lam1
  Rrs745 = Rrs_$Rrs745 # lam2
  Rrs862 = Rrs_$Rrs862 # lam3
  
  Rrs_W = data.frame(Rrs671, Rrs745, Rrs862)
  
  C0 = 0.04
  C1 = 1.17
  C2 = 0.4
  C3 = 14.86
  a1 = 20.43
  a2 = 2.15
  
  W = Rrs_W / rowSums(Rrs_W)
  
  GI = C0 * Rrs551 / Rrs486 + 
    C1 * W[,1] * Rrs671 / Rrs551 + 
    C2 * W[,2] * Rrs745 / Rrs551 + 
    C3 * W[,3] * Rrs865 / Rrs551
  
  SPM = a1 * GI ^ a2
  
  return(list(
    param = list(C0 = C0, C1 = C1, C2 = C2, C3 = C3, a1 = a1, a2 = a2),
    W  = W,
    GI = GI,
    SPM = SPM
  ))

}











