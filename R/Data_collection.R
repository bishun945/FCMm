#' @name SRF_LIST
#' @title Collection of spectral responding functions
#' @description A list including spectral responding functions of several popular sensors
#' @docType data
#' @keywords datasets
#' @usage SRF_LIST
#' @format List with four elements: sensor (character), srf (data.frame), cw_med (num), and cw_max (num)
#' @details The cw_med and cw_max are different since 
#'   we use two functions \code{find_center_wavelength_med} and 
#'   \code{find_center_wavelength_max} to obtain their center wavelength.
#'   \itemize{
#'     \item \code{find_center_wavelength_med} center wavelength is determined by the wavelength
#'       position of half-maximum width
#'     \item \code{find_center_wavelength_max} wavelength at the position of maximum srf
#'   }
#'   All included sensors are Sentinel-2A, Landsat-8, GF1-WFV1, GF1-WFV2, GF1-WFV3, GF1-WFV4, 
#'     GF2-PMS1, GF2-PMS2, GF4-PMI, HJ1A-CCD1, HJ1A-CCD2, HJ1B-CCD1, HJ1B-CCD2, MODISA, MERIS, 
#'     VIIRS, OLCI, GOCI, and GF5-VIMS.
#'   The data are collected from their corresponding official websites (see in references).
#'   
#' @references 
#'   \itemize{
#'     \item http://www.cresda.com/CN/
#'     \item https://oceancolor.gsfc.nasa.gov
#'     \item https://sentinels.copernicus.eu/web/sentinel/home
#'   }
#' @family Datasets
NULL



#' @name OLCI_TH
#' @title A test OLCI image raster in Lake Taihu
#' @description A dataframe contains OLCI_TH data
#' @docType data
#' @keywords datasets
#' @usage OLCI_TH
#' @format A RasterBrick with 791136 elements
#' @details This dataset is a OLCI raster in Lake Taihu for testing the package `FCMm`.  
#' @note The dataset should be only used in this package test. If any want 
#'   to use it in your own study, he or she should contact us via bishun1994 at foxmail or 
#'   liyunmei at njnu.edu.cn
#'   
#'   This OLCI image onboard Sentinel-3A (ESA) in Lake Taihu was scened on 17 August, 2018
#'   
#'   
#'   
#' @references 
#' * Bi S, Li Y, Xu J, et al. Optical classification of inland waters based on
#'   an improved Fuzzy C-Means method[J]. Optics Express, 2019, 27(24): 34838-34856.
#'   
#' * Liu G, Li Y, Lyu H, et al. An improved land target-based atmospheric correction method 
#'   for Lake Taihu[J]. IEEE Journal of Selected Topics in Applied Earth Observations and
#'   Remote Sensing, 2015, 9(2): 793-803.
#' @family Datasets
NULL



#' @name Nechad2015
#' @title Collection of in_situ data by Nechad et al. (2015)
#' @description A dataframe contains Nechad2015 data
#' @docType data
#' @keywords datasets
#' @usage Nechad2015
#' @format A dataframe with 336 rows by 18 cols
#' @details None
#' @note Colnames are "DataProvider, SAMPLE.ID, X412.5, X442.5, X490, X510, X560, X620,
#'   X665, X681.25, X708.75, X.Chl_a..ug.L., X.TSM..mg.l., CC_SITE, DATE.dd.mm.yyyy.,
#'  time, latitude, longitude"
#' @references  B. Nechad, K. Ruddick, T. Schroeder, K. Oubelkheir, D. Blondeau-Patissier, N. Cherukuru, 
#'   V. Brando, A. Dekker, L. Clementson, and A. C. Banks, "CoastColour Round Robin data sets: a database to
#'   evaluate the performance of algorithms for the retrieval of water quality parameters in coastal waters,"
#'   Earth system science data 7, 319-348 (2015). (https://doi.pangaea.de/10.1594/PANGAEA.841950)
#' @family Datasets
NULL



#' @name AeronetOC2019
#' @title Collection of coastal AeronetOC water spectra till 2019
#' @description A dataframe contains AeronetOC data
#' @docType data
#' @keywords datasets
#' @usage AeronetOC2019
#' @format A dataframe with 10667 rows by 14 cols
#' @details AeronetOC is collected from Level 2 products provided by Aeronet Ocean Color stations located
#'   in great inland lakes (i.e., Lake_Erie, Palgrunden) and coastal waters (i.e., COVE_SEAPRISM,
#'   Galata_Platform, Gloria, Gustav_Dalen_Tower, Helsinki_Lighthouse, LISCO, MVCO, and Zeebrugge-MOW1)
#'    with bands at 410, 440, 490, 530, 550, 667, 869 nm. 
#' @note Colnames are "SampleID, X410nm, X440nm, X490nm, X530nm, 
#'   X550nm, X667nm, X869nm, X1020nm, Pressure, Wind_Speed, Chlorophyll.a, Sea_Surface_Reflectance, Ozone"
#' @references  All these data are available at https://aeronet.gsfc.nasa.gov.
#' @family Datesets
NULL



#' @name WaterSpec35
#' @title A test-purpose collection of inland water data
#' @description A dataframe contains WaterSpec35 data
#' @docType data
#' @keywords datasets
#' @usage WaterSpec35
#' @format A dataframe with 35 rows by 17 cols
#' @details This dataset is a part of research of Bi et al. (2019) for
#'   testing the package `FCMm`. Colnames are "SampleID, Chla,
#'   400, 412.5, 442.5, 490, 510, 560, 620,
#'   665, 673.75, 681.25, 708.75, 753.75, 778.75, 865, 885"
#' @note The dataset should be only used in this package test. If any want 
#'   to use it in your own study, he or she should contact us via bishun1994 at foxmail or 
#'   liyunmei at njnu.edu.cn
#' @references Bi S, Li Y, Xu J, et al. Optical classification of inland waters based on
#'   an improved Fuzzy C-Means method[J]. Optics Express, 2019, 27(24): 34838-34856.
#' @family Datasets
NULL



#' @name Valente2019
#' @title Collection of in_situ data by Valente et al. (2019)
#' @description A dataframe contains Valente2019 data
#' @docType data
#' @keywords datasets
#' @usage Valente2019
#' @format A dataframe with 1205 rows by 14 cols
#' @details None
#' @note Colnames are "Date.Time, Lat, Lon, Depth.m, Chla.1, Chla.2,
#'   X412nm, X443nm, X490nm, X510nm, X560nm, X620nm, X665nm, X681nm"
#' @references A. Valente, S. Sathyendranath, V. Brotas, S. Groom, M. Grant, M. Taberner,
#'   D. Antoine, R. Arnone, W. M. Balch, and K. Barker, "A compilation of global bio-optical
#'   in situ data for ocean-colour satellite applications," Earth System Science Data 8, 235-252 (2016).
#'   (https://doi.org/10.1594/PANGAEA.898197)
#' @family Datasets
NULL



#' @name dt_water
#' @title Absorption and scattering coefficient of pure water
#' @description Absorption and scattering coefficient of pure water with unit m^-1
#' @docType data
#' @keywords datasets
#' @usage dt_water
#' @format A dataframe with 561 rows by 4 cols
#' @details The wavelength is from 340 nm to 900 nm.
#' @note Will add later.
#' @references Will add later.
#' @family Datasets
NULL



#' @name Rrs_clusters.default
#' @title Seven built-in clusters
#' @description A dataframe contains the water optical clusters produced by Bi et al. (2019)
#' @docType data
#' @keywords datasets
#' @usage Rrs_clusters.default
#' @format dataframe with 7 rows by 15 cols
#' @note You have to take a look at the description of this built-in cluster (water type) such as
#'   their corresponding water color parameters in Bi et al. (2019).
#' @family Datasets
NULL



#' @name wavelength.default
#' @title Wavelength of built-in clusters
#' @description A vector contains the wavelength of water optical clusters produced by Bi et al. (2019)
#' @docType data
#' @keywords datasets
#' @usage wavelength.default
#' @format a vector with 15 elements
#' @family Datasets
NULL



