#' @import ppclust
#' @importFrom magrittr %>% %<>%
#' @import reshape2
#' @import ggplot2
#' @import raster
#' @import viridis
#' @import ggthemes
#' @import heatmaply
#' @import inaparc
#' @import stringr
#' @import fclust
#' @import tidyverse
#' @importFrom gridExtra arrangeGrob
#' @importFrom utils write.csv data
#' @importFrom readxl read_excel excel_sheets

utils::globalVariables(c('.','band','value','name','wv','variable',
                         'kmpp', 'imembrand','x','y','r','g','b',
                         'res.im.cluster','Chla','Rrs_clusters.default',
                         'wavelength.default','nm',
                         'SRF_LIST'))
