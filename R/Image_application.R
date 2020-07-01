#' @title Apply FCM_m to raster data
#' @name apply_to_image
#' @description
#'   This function could apply the defined water cluster to corrected image files.
#'   Should run \link{generate_param} to generate a \code{res} list as an input
#'   of function \code{apply_to_image}
#'
#' @param input A \strong{raster} or a \strong{character} linking
#'   to the raster file on the disk.
#' @param res A required list that used for clustering the image data including:
#'   \itemize{
#'     \item \strong{FD}  The results from function \code{FuzzifierDetermination}
#'       providing wavelength
#'     \item \strong{K}  Cluster number
#'     \item \strong{res.FCM}  Cluster center and used fuzzifier value.
#'   }
#'   For the convenience, function \link{generate_param} supports to quickly
#'     generate this \code{list}. See more in examples.
#' @param color_palette The palette of cluster color. Default as \code{RdYlBu(res$K)}.
#'   In \code{FCMm}, it could be \link{RdYlBu}, \link{Spectral}, \link{HUE} or other color values
#'   with same length of cluster number.
#' @param output_image Logical, whether to save images as raster files to your disk. Default as \code{FALSE}.
#' @param output_resultpng Logical, whether to save png files to your disk. Default as \code{FALSE}.
#' @param output_imRrs.n Logical, whether to produce normalized Rrs files. Default as \code{FALSE}. 
#'   This parameter is used for inspection of FCM trainning. Will be \strong{deprecated} in the following version.
#' @param output_dir The directory of output files. Default as the current working directory (\code{"."})
#' @param output_format A string, the format of raster file, default as \code{"GTiff"}.
#'   See more in \code{\link[raster:writeRaster]{writeRaster}}
#' @param Chla_est Logical, whether to estimate Chla concentration using the function \link{FCM_m_Chla_estimation}.
#'   Default as \code{FALSE}.
#' @param title.name Character, the title name of ggplot for plotting. Default as \code{NULL}.
#' @param png_scale Numeric, scale of png. Recommended as \code{50}.
#' @param fn_memb A string, filename of membership raster file. Default as \code{"output_membership"}.
#' @param fn_cluster A string, filename of cluster raster file. Default as \code{"output_cluster"}.
#' @param fn_imRrs.n A string, filename of normalized Rrs raster file. Default as \code{"output_imRrs_normalized"}.
#' @param fn_truecolorpng A string, file name of truecolor png. Default as \code{"output_truecolor"}.
#' @param fn_Chla A string, file name of estimated Chla raster file. Default as \code{"output_Chla"}.
#'
#' @return A \code{list()} of all results and several inputs:
#'   \itemize{
#'     \item \strong{input}  What we input character link to image file or raster object
#'     \item \strong{res}  Condition input of apply_to_image
#'     
#'     \item \strong{res.FCM}  List of FCM.new result
#'     \item \strong{res.Chla}  Data.frame includes the template and final Chla estimation results
#'     
#'     \item \strong{coord}  Dataframe of xy coordinates
#'     \item \strong{imRrs.raw}  Dataframe of raw Rrs excluding NA pixels. Thus, it may have less pixels than 
#'       that of \code{input}
#'     \item \strong{imRrs.n}  Dataframe of normalized Rrs. Same as \code{imRrs.raw} but in the normalized scale
#'     \item \strong{im.cluster} Dataframe of cluster with coordinates.
#'     \item \strong{im.memb} Dataframe of membership returned by \link{apply_FCM_m}. It includes coordinates.
#'     
#'     \item \strong{raster.memb}  Raster object of membership value. \code{NULL} if \code{output_image} is \code{FALSE}.
#'     \item \strong{raster.cluster}  Raster object of cluster. \code{NULL} if \code{output_image} is \code{FALSE}.
#'     \item \strong{raster.Chla} Raster object of Chla concentration. \code{NULL} if \code{output_image} is \code{FALSE}.
#'     
#'     \item \strong{p.memb}  A ggplot list of membership value map
#'     \item \strong{p.cluster}  A ggplot list of cluster map
#'     \item \strong{p.truecolor} A ggplot list of truecolor map
#'     \item \strong{p.Chla}  A ggplot list of Chla map
#'   }
#'
#' @export
#' 
#' @rdname apply_to_image
#' 
#' @note 
#' 
#'   \strong{2019-12-12}: 
#' 
#'   The section FCM running used the subset of default Rrs clusters.
#'   Please see the section \strong{Part II: New coming raster data} by running \code{vignette('Builtin_centrodis')}
#'   if have not known how to do in that situation.
#'   
#'   Also, if it is your first time to get the image data into \strong{R}, you could
#'   load the raster data by typing \code{raster::brick(filename)} which filename
#'   is the path of your data such as '/data/Test.dat' or 'E://data//Test.tiff' or so.
#'   
#'   \strong{2020-06-27}: 
#'   
#'   Note that the inputs of \code{apply_to_image} support the file path of raster file.
#'
#' @examples 
#' \donttest{
#' library(FCMm)
#' data("OLCI_TH")
#' data("Bi_clusters")
#' res <- generate_param(c(413,443,490,510,560,620,665,674,709,754,865,885))
#' im_result <- apply_to_image(input=OLCI_TH, res=res, title.name="Test_image", Chla_est=TRUE)
#' }
#'   
#' @references 
#' Bi S, Li Y, Xu J, et al. Optical classification of inland waters based on
#'   an improved Fuzzy C-Means method[J]. Optics Express, 2019, 27(24): 34838-34856.
#' @family Fuzzy cluster functions
#' 
#' @import ggplot2
#' @importFrom raster raster brick as.data.frame rasterFromXYZ stretch rasterFromXYZ 
#'   crs crs<- writeRaster
#' @importFrom reshape2 melt
#' @importFrom magrittr %>% %<>%
#' @importFrom ggthemes theme_map
#' @importClassesFrom raster Raster RasterBrick RasterStack
#' 
apply_to_image <- function(input, res, 
                           color_palette = RdYlBu(res$K),
                           output_image=FALSE, output_resultpng=FALSE, output_imRrs.n=FALSE,
                           output_dir = ".", output_format="GTiff",
                           Chla_est=FALSE,
                           title.name = NULL, png_scale=50,
                           fn_memb="output_membership",
                           fn_cluster="output_cluster",
                           fn_imRrs.n="output_imRrs_normalized",
                           fn_truecolorpng="output_truecolor",
                           fn_Chla = "output_Chla"){

  if(missing(input))
    stop("The input of image file is missing!")
  if(missing(res))
    stop("FCM results is missing")
  
  if(length(color_palette) != res$K)
    stop("The length of color_palette shoud be same with cluster number!")
  
  if(!dir.exists(output_dir))
    stop("The specified output directory does not exist!")
  
  fn_memb <- file.path(output_dir, fn_memb)
  fn_cluster <- file.path(output_dir, fn_cluster)
  fn_imRrs.n <- file.path(output_dir, fn_imRrs.n)
  fn_truecolorpng <- file.path(output_dir, fn_truecolorpng)
  fn_Chla <- file.path(output_dir, fn_Chla)
  
  message("Since we have not check the input image file,")
  message("  please make sure the wavelength of image file and cluster dataframe match correspondingly.")
  message("The normalization of spectra is default in this version!")

  wv <- res$FD$wv
  if(class(input)[1] == "character"){
    im <- brick(input)
  }else if(class(input)[1] == "RasterBrick" | class(input)[1] == "RasterStack"){
    im <- input
  }else{
    stop('Unknow input of image. Make sure the input is charater or RasterBrick')
  }

  if(im@data@nlayers != length(res$FD$wv))
    stop("The band number of image file is different from wavelength length!")
  imdf <- raster::as.data.frame(im, na.rm=TRUE, xy=TRUE)
  
  x_name <- which(names(imdf) == "x")
  y_name <- which(names(imdf) == "y")
  imdf <- cbind(imdf[,c(x_name, y_name)], imdf[,-c(x_name, y_name)])
  
  coord <- data.frame(x = imdf$x, y = imdf$y)
  imRrs.raw <- imdf[,-c(1,2)]
  names(imRrs.raw) <- paste0("Rrs",wv)
  Rrs <- as.matrix(imRrs.raw)
  Area <- trapz(wv,Rrs)
  Area <- base::as.data.frame(Area)
  imRrs.n <- imRrs.raw
  for(i in 1:ncol(imRrs.raw))
    imRrs.n[,i] = imRrs.raw[,i] / Area
  
  raster.imRrs = rasterFromXYZ(cbind(coord, imRrs.raw), crs = im@crs)
  
  # save imRrs_normalized dataframe to ENVI rasterBrick
  if(output_imRrs.n){
    raster.imRrs.n <- rasterFromXYZ(data.frame(x=imdf$x, y=imdf$y, imRrs.n), crs = im@crs)
    crs(raster.imRrs.n) <- im@crs
    writeRaster(raster.imRrs.n, fn_imRrs.n, format='ENVI', overwrite=TRUE)
    message(paste0("Normalized Rrs map was generated, named: ", fn_imRrs.n))
  }

  # Apply FCM-m to image dataframe
  v <- res$res.FCM$v
  m <- res$res.FCM$m
  message('Apply fcm to image dataframe ......')
  res.FCM <- apply_FCM_m(Rrs=imRrs.raw, wavelength=wv, Rrs_clusters=v, m_used = m, 
                         stand = FALSE, default.cluster=FALSE)
  res.im <- list()
  res.im$u <- base::as.data.frame(res.FCM$u)
  res.im$cluster <- res.FCM$cluster

  # Save true color image
  rgb <- brick(im[[10]],im[[7]],im[[5]],im[[2]])
  rgb_stretch <- stretch(x=rgb, minv=0, maxv=255)
  rgb_df <- raster::as.data.frame(rgb_stretch, xy=TRUE)
  rgb_df <- data.frame(x=rgb_df$x, y=rgb_df$y,
                       n=rgb_df[,3],r=rgb_df[,4], g=rgb_df[,5],b=rgb_df[,6]) %>% na.omit
  p.truecolor = ggplot(data=rgb_df) +
    geom_raster(aes(x=x,y=y,fill=rgb(r,g,b, maxColorValue=255)), hjust=0.5, vjust=0.5) +
    scale_fill_identity() +
    labs(title=title.name) +
    coord_equal() + theme_map() +
    theme(plot.title=element_text(hjust=0.5),
          text = element_text(size=13))
  print(p.truecolor)

  # Save membership rasters
  sub.memb <- data.frame(imdf$x,imdf$y, res.im$u)
  names(sub.memb) <- c("x","y",paste('Cluster',seq(1,res$K)))
  if(output_image){
    raster.memb <- rasterFromXYZ(sub.memb, crs = im@crs)
    writeRaster(raster.memb, fn_memb, format=output_format, overwrite=TRUE)
    message(paste0("Membership map was generated, named: ", fn_memb))
  }else{
    raster.memb <- NULL
  }
  message("Plotting membership values ......")
  p.memb<- ggplot() +
    geom_raster(data=melt(sub.memb,id=c('x','y')), aes(x=x,y=y,fill=value),
                hjust=0.5, vjust=0.5) +
    scale_fill_viridis_c(na.value='white',begin=0,end=1,
                       breaks=c(0.0,0.25,0.5,0.75,1.0)) +
    facet_wrap(~variable, nrow=3) + 
    labs(fill='Membership', title=title.name) +
    coord_equal() + theme_map() +
    theme(text = element_text(size=13),
          plot.title=element_text(hjust=0.5),
          # legend.position=c(1,0),
          legend.position = "right",
          legend.key.width=unit(1,"lines"),
          legend.justification=c(0.5, 0.5),
          strip.background=element_rect(fill='white',color='white'))
  print(p.memb)

  # Save cluster raster
  im.cluster <- data.frame(imdf$x,imdf$y, cluster = res.im$cluster %>% as.character)
  names(im.cluster) <- c("x","y", "cluster")
  if(output_image){
    raster.cluster <- rasterFromXYZ(im.cluster, crs = im@crs)
    writeRaster(raster.cluster, fn_cluster, format=output_format, overwrite=TRUE)
    message(paste0("Cluster map was generated, named: ", fn_cluster))
  }else{
    raster.cluster <- NULL
  }
  message("Plotting clusters ......")
  im.cluster$cluster_f <- factor(im.cluster$cluster, levels=1:nrow(v), ordered = TRUE)
  cp <- color_palette
  names(cp) <- levels(im.cluster$cluster_f)
  p.cluster <- ggplot() +
    geom_raster(data=im.cluster,aes(x=x,y=y,fill=cluster_f),
                hjust=0.5, vjust=0.5) +
    scale_fill_manual(values=cp) +
    labs(fill='Cluster', title=title.name) +
    coord_equal() + theme_map() +
    theme(text=element_text(size=13),
          plot.title=element_text(hjust=0.5),
          legend.position='right',
          legend.justification='center',
          legend.key.width=unit(1.5,"lines"))
  print(p.cluster)

  # obtain width and height of ggsave
  width <- raster.imRrs@ncols/png_scale * 1.1
  height<- raster.imRrs@nrows/png_scale

  res.Chla <- NULL
  p.Chla <- NULL
  raster.Chla <- NULL
  # Chla estimation and plot Chla map
  if(Chla_est){
    X <- data.frame(Rrs665 = res.FCM$x$Rrs665,
                    Rrs709 = res.FCM$x$Rrs709,
                    Rrs754 = res.FCM$x$Rrs754,
                    M=res.FCM$u)
    res.Chla <- FCM_m_Chla_estimation(Rrs=X[,1:3],U=X[,4:10])
    sub.Chla <- data.frame(imdf$x,imdf$y, res.Chla$conc.Blend)
    names(sub.Chla) <- c("x","y","Chla")
    if(output_image){
      raster.Chla <- rasterFromXYZ(sub.Chla, crs = im@crs)
      writeRaster(raster.Chla, fn_Chla, format=output_format, overwrite=TRUE)
      message(paste0("Chla concentration map was generated, named: ", fn_Chla))
    }else{
      raster.Chla <- NULL
    }
    message("Plotting Chla concentration ......")
    oldoptions <- options(scipen=1000)
    on.exit(options(oldoptions))
    p.Chla<- ggplot() +
      geom_raster(data=sub.Chla, aes(x=x,y=y,fill=Chla),
                  hjust=0.5, vjust=0.5) +
      labs(fill='Chla [ug/L]', title=title.name) +
      coord_equal() + theme_map() +
      theme(plot.title=element_text(hjust=0.5),
            legend.position='right',
            legend.justification='center',
            legend.key.height=unit(2.5,"lines"),
            strip.background=element_rect(fill='white',color='white'))
    if(max(sub.Chla$Chla, na.rm=TRUE) >= 800){
      p.Chla <- p.Chla +
        scale_fill_viridis_c(na.value='gray', trans='log10',
                           limits=c(1,800),
                           breaks=c(10^seq(0,2.9,0.5)) %>% round(.,2))
    }else{
      p.Chla <- p.Chla +
        scale_fill_viridis_c(na.value='gray', trans='log10',
                           breaks=c(10^seq(0,2.9,0.5)) %>%
                             .[. < max(sub.Chla$Chla, na.rm=TRUE)] %>% round(.,2))
    }
    print(p.Chla)
  }

  # Save ggplot png for membership and cluster patterns
  if(output_resultpng){
    ggsave(paste0(fn_memb,'.png'), plot=p.memb, device='png')
    ggsave(paste0(fn_cluster,'.png'), plot=p.cluster, device='png', width=width, height=height, units='in')
    ggsave(paste0(fn_truecolorpng, '.png'), plot=p.truecolor, device='png', width=width, height=height, units='in')
    ggsave(paste0(fn_Chla, '.png'), plot=p.Chla, device='png', width=width, height=height, units='in')
  }

  message('Saving the results to the list ......')
  message(fn_memb)
  message(fn_cluster)

  # Save results to one list
  result = list()
  
  # input part
  result$input <- input
  result$res <- res
  
  # FCM and blending work
  result$res.FCM <- res.FCM
  result$res.Chla <- res.Chla
  
  # imdf results
  result$coord <- coord
  result$imRrs.raw <- imRrs.raw
  result$imRrs.n <- imRrs.n
  result$im.cluster <- im.cluster
  result$im.memb    <- sub.memb
  
  # raster outputs
  result$raster.memb <- raster.memb
  result$raster.cluster <- raster.cluster
  result$raster.Chla    <- raster.Chla
  
  # plots
  result$p.memb <- p.memb
  result$p.cluster <- p.cluster
  result$p.truecolor <- p.truecolor
  result$p.Chla <- p.Chla

  message('Done!')

  return(result)
}

#' @title Generate the input param \code{res} of \code{apply_to_image} based on 
#'   the built-in or user-defined centroids
#' @name generate_param
#' @param wl wavelength of subsetting
#' @export
#' @note \code{generate_param} only support the cluster centroids proposed by \code{Bi et al. (2019)}.
#'   However, the \code{generate_param_ex} could be used for user-defined centroids if you want.
#' @return The return of \code{generate_param} or \code{generate_param_ex} is a list used for 
#'   function \link{apply_to_image}.
#' @family Fuzzy cluster functions
#' @references 
#' \itemize{
#'   \item Bi S, Li Y, Xu J, et al. Optical classification of inland waters based on
#'     an improved Fuzzy C-Means method[J]. Optics Express, 2019, 27(24): 34838-34856.
#' }
#' @examples 
#' library(FCMm)
#' wl = c(413, 443, 490, 510, 560, 620, 665, 674, 681, 709, 754, 779, 865, 885)
#' res = generate_param(wl)
#' 
generate_param <- function(wl){
  w <- (wavelength.default %in% wl)
  wavelength <- wavelength.default[w]
  Rrs_clusters <- Rrs_clusters.default[,w]
  # generate the required res
  res <- list()
  res$FD$wv <- wavelength
  res$K <- nrow(Rrs_clusters)
  res$res.FCM <- list(v=Rrs_clusters,m=1.36)
  return(res)
}

#' @rdname generate_param
#' @param centroids The centorids values required for fuzzy clustering on images. Its colnames should
#'   be converted to numeric wavelength. The format should be like \code{Rrs_clusters.default}.
#' @param m The fuzzifier value of FCM. Default as \code{1.36}. 
#' @importFrom stringr str_detect
#' @export
generate_param_ex <- function(centroids, m = 1.36){
  
  wv = colnames(centroids)
  w  = str_detect(wv, "[:digit:]|\\.")
  if(sum(w) < ncol(centroids)){
    warning(sprintf(
      "Following columns are ignored as their names cannot be converted to numeric:\n%s",
      paste0(wv[!w], collapse = " ")
    ))
    centroids = centroids[, w]
  }
  
  wv = as.numeric(colnames(centroids))
  res <- list()
  res$FD$wv   = wv
  res$K       = nrow(centroids)
  res$res.FCM = list(v = centroids, m=m)
  
  return(res)
}


