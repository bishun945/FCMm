#' @name FCMm
#' 
#' @title Fuzzy Cluster Method Based on the Optimized m Value (Fuzzifier)
#' 
#' @description 
#' \code{FCMm} is a package for fuzzy clustering water spectra (or called water color). 
#'   Given that the most of water color spectra data sets are considered as the high dimensional 
#'   set, the advantage of this method is making FCM assign the membership (sum as 1) harder, 
#'   ensuring the desired water type are restricted to its belongings (not too soft). It is 
#'   possible to cluster the harm algal bloom water type which can not be produced by FCM with \code{m=2}.
#' 
#' \itemize{
#'   \item If you want to cluster your own data sets, it provides an improved Fuzzy Cluster 
#'     Method (FCM) by optimizing the fuzzifier value (default but not good being 2). 
#'   \item You can also use the built-in cluster of inland waters produced by Bi et al. (2019) 
#'     and can simply obtain the Chlorophyll-a concentration by blending three algorithms with 
#'     relatively low bias.  
#'   \item It supports raster (or called imagery) processing (see more details in help documents 
#'     or vignettes).
#'   \item It includes several data sets about water color spectra and corresponding water quality 
#'     parameters and a testing image raster (see help documents for details).
#'   \item Please see NEWS to get changes in each version.
#' }
#' 
#' @seealso Useful vignettes:
#'   \itemize{
#'     \item \code{vignette('Builtin_centrodis')}
#'     \item \code{vignette('Cluster_new_data')}
#'     \item \code{vignette('Assessment')}
#'   }
#' 
#' @docType package
#' 
#' @references 
#' \itemize{
#'   \item Bi S, Li Y, Xu J, et al. Optical classification of inland waters based on 
#'     an improved Fuzzy C-Means method[J]. Optics Express, 2019, 27(24): 34838-34856.
#'   \item Jackson T, Sathyendranath S, Mélin F. An improved optical classification 
#'     scheme for the Ocean Colour Essential Climate Variable and its applications[J]. 
#'     Remote Sensing of Environment, 2017, 203: 152-161.
#'   \item Moore T S, Dowell M D, Bradt S, et al. An optical water type framework for 
#'     selecting and blending retrievals from bio-optical algorithms in lakes and coastal 
#'     waters[J]. Remote sensing of environment, 2014, 143: 97-111.
#'   \item Spyrakos E, O'Donnell R, Hunter P D, et al. Optical types of inland and coastal 
#'     waters[J]. Limnology and Oceanography, 2018, 63(2): 846-870.
#'   \item Dembele D. Multi-objective optimization for clustering 3-way gene expression 
#'     data[J]. Advances in Data Analysis and Classification, 2008, 2(3): 211-225.
#' }
"_PACKAGE"

