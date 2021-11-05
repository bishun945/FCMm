
#' @name Spectra.Centroids
#' @title Several built-in spectral centroids
#' @description A data frame contains the water optical clusters produced by Bi et al. (2019),
#'   Bi et al. (2021), Bi's PhD thesis paper in 2021, Jackson et al. (2017), and Moore et al. (2014)
#' @docType data
#' @keywords datasets
#' @usage 
#' data("B2019_norm")
#' data("B2021_norm")
#' data("BPHD_norm")
#' data("Jac17_raw")
#' data("Moo14_raw")
#' @family Datasets
#' @note
#' \code{Bi2019_norm}, \code{Bi2021_norm}, \code{BPHD_norm} are normalized centroids 
#'   based on spectral integrals.
#' 
#' \code{Jac17_raw} is the Rrs centroids based on the raw scale (using Mahalanobis distance)
#' 
#' \code{Moo14_raw} is converted from the Rrs(0-) by using Rrs(0+) = 0.52 / (1/Rrs(0-) - 1.7)
#'   
#' @references 
#' 
#' Bi, S., Li, Y., Liu, G., Song, K., Xu, J., Dong, X., Cai, X., Mu, M., Miao, S., & Lyu, H. (2021). 
#'   Assessment of Algorithms for Estimating Chlorophyll-a Concentration in Inland Waters: 
#'   A Round-Robin Scoring Method Based on the Optically Fuzzy Clustering. IEEE Transactions on 
#'   Geoscience and Remote Sensing.
#' 
#' Bi, S., Li, Y., Xu, J., Liu, G., Song, K., Mu, M., Lyu, H., Miao, S., & Xu, J. (2019). 
#'   Optical classification of inland waters based on an improved Fuzzy C-Means method. 
#'   Optics Express, 27(24), 34838–34856.
#' 
#' Jackson, T., Sathyendranath, S., & Mélin, F. (2017). An improved optical classification 
#'   scheme for the Ocean Colour Essential Climate Variable and its applications. Remote 
#'   Sensing of Environment, 203, 152–161.
#' 
#' Moore, T. S., Dowell, M. D., Bradt, S., & Verdu, A. R. (2014). An optical water type 
#'   framework for selecting and blending retrievals from bio-optical algorithms in
#'   lakes and coastal waters. Remote Sensing of Environment, 143, 97–111.
#' 
NULL

#' @name aowt
#' @title All owts in this package
#' @param plot Set \code{TRUE} to plot these centroids
#' @importFrom dplyr mutate
#' @importFrom cowplot plot_grid
#' @noRd
aowt <- function(plot = FALSE) {
  
  aowt_name <- c("B2019_norm", "B2021_norm", "BPHD_norm", "Jac17_raw", "Moo14_raw")
  
  if(plot) {
    
    owts <- list(
      B2019_norm, B2021_norm, BPHD_norm, Jac17_raw, Moo14_raw
    )
    
    lapply(owts, function(x) {
      names(x)[1] <- "OWT"
      x %>% 
        reshape2::melt(id = 1, variable.name = "wv") %>%
        dplyr::mutate(wv = as.numeric(as.vector(wv))) %>%
        ggplot(aes(x=wv, y=value, color = OWT)) + 
        geom_path()
    }) %>%
      cowplot::plot_grid(plotlist = . , labels = aowt_name)
    
  }
  
  return(aowt_name)
  
}

#' @name wavelength.default
#' @title Wavelength of built-in clusters
#' @description A vector contains the wavelength of water optical clusters produced by Bi et al. (2019)
#' @docType data
#' @keywords datasets
#' @usage wavelength.default
#' @format a vector with 15 elements
#' @family Datasets
NULL
