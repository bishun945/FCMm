#' @title apply_to_image
#' @name apply_to_image
#' @description
#'   This function could apply the defined water cluster to corrected image files.
#'   Should run \code{generate_param()} to generate a \code{res} list as an input
#'   of function \code{apply_to_image}
#'
#' @usage apply_to_image(input, res,
#'   output_image=TRUE, output_resultpng=FALSE, output_imRrs.n=FALSE,
#'   Chla_est=FALSE,
#'   title.name = NULL, png_scale=50,
#'   fn_memb="output_membership",
#'   fn_cluster="output_cluster",
#'   fn_imRrs.n='output_imRrs_normalized',
#'   fn_truecolorpng='output_truecolor',
#'   fn_Chla = 'output_Chla',
#'   output_format='GTiff')
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
#'   For the convenience, function \code{generate_param()} supports to quickly
#'     generate this \code{list}. See more in examples.
#' @param output_image Logical, whether to produce image files
#' @param output_resultpng Logical, whether to produce png files
#' @param output_imRrs.n Logical, whether to produce normalized Rrs files
#' @param Chla_est Logical, whether to estimate Chla concentration
#' @param title.name Character, the title name of ggplot for plotting
#' @param png_scale Numeric, scale of png
#' @param fn_memb A string, filename of membership raster file
#' @param fn_cluster A string, filename of cluster raster file
#' @param fn_imRrs.n A string, filename of normalized Rrs raster file
#' @param fn_truecolorpng A string, file name of truecolor png
#' @param fn_Chla A string, file name of estimated Chla raster file
#' @param output_format A string, the format of raster file, default as \code{GTiff}.
#'   See more in \code{\link[raster:writeRaster]{writeRaster}}
#'
#' @return A \code{list()} of all results and several inputs:
#'   \itemize{
#'     \item \strong{input}  What we input character link to image file or raster object
#'     \item \strong{res}  Condition input of apply_to_image
#'     \item \strong{raster.memb}  Raster object of membership value
#'     \item \strong{raster.cluster}  Raster object of cluster
#'     \item \strong{p.memb}  A ggplot list of membership value map
#'     \item \strong{p.cluster}  A ggplot list of cluster map
#'     \item \strong{p.truecolor} A ggplot list of truecolor map
#'     \item \strong{p.Chla}  A ggplot list of Chla map
#'     \item \strong{imdf}  Data.frame of raw Rrs including xy coordinates
#'     \item \strong{imRrs.raw}  Data.frame of raw Rrs
#'     \item \strong{imRrs.n}  Data.frame of normalized Rrs
#'     \item \strong{res.FCM}  List of FCM.new result
#'     \item \strong{res.Chla}  Data.frame includes the template and final Chla estimation results
#'   }
#'
#' @export
#' 
#' @rdname apply_to_image
#' 
#' @details
#'   The section FCM running used the subset of default Rrs clusters.
#'   Please see the section \strong{One more thing} in vignettes \strong{New_Data_Running_FCMm}
#'   if have not known how to do in that situation.
#'   
#'   Also, if it is your first time to get the image data into \strong{R}, you could
#'   load the raster data by typing \code{raster::brick(filename)} which filename
#'   is the path of your data such as '/data/Test.dat' or 'E://data//Test.tiff' or so.
#'
#' @examples 
#' \dontrun{
#' library(FCMm)
#' data("OLCI_TH")
#' data("Bi_clusters")
#' res <- generate_param(c(413,443,490,510,560,620,665,674,709,754,865,885))
#' im_result <- apply_to_image(input=OLCI_TH, res=res,
#'   title.name="Test_image", Chla_est=T, output_image=F)
#' }
#'   
#' @references 
#' Bi S, Li Y, Xu J, et al. Optical classification of inland waters based on
#'   an improved Fuzzy C-Means method[J]. Optics Express, 2019, 27(24): 34838-34856.
#'  

apply_to_image <- function(input, res,
                           output_image=TRUE, output_resultpng=FALSE, output_imRrs.n=FALSE,
                           Chla_est=FALSE,
                           title.name = NULL, png_scale=50,
                           fn_memb="output_membership",
                           fn_cluster="output_cluster",
                           fn_imRrs.n='output_imRrs_normalized',
                           fn_truecolorpng='output_truecolor',
                           fn_Chla = 'output_Chla',
                           output_format='GTiff'){

  if(missing(input))
    stop("The input of image file is missing!")
  if(missing(res))
    stop("FCM results is missing")

  message("Since we have not check the input image file,")
  message("  please make sure the wavelength of image file and cluster dataframe match correspondingly.")
  message("The normalization of spectra is default in thie version!")

  wv <- res$FD$wv
  if(class(input)[1] == "character"){
    im <- brick(input)
  }else if(class(input)[1] == "RasterBrick" | class(input)[1] == "RasterStack"){
    im <- input
  }else{
    stop('Unknow input of image. Make sure the input is charater or RasterBrick')
  }

  if(im@data@nlayers!=length(res$FD$wv))
    stop("The band number of image file is different from wavelength length!")
  imdf <- as.data.frame(im,na.rm=T,xy=T)
  imRrs.raw <- imdf[,-c(1,2)]
  names(imRrs.raw) <- paste0("Rrs",wv)
  Rrs <- as.matrix(imRrs.raw)
  Area <- .trapz(wv,Rrs)
  imRrs.n <- imRrs.raw
  for(i in 1:ncol(imRrs.raw))
    imRrs.n[,i] = imRrs.raw[,i] / Area

  # save imRrs_normalized dataframe to ENVI rasterBrick
  if(output_imRrs.n){
    imRrs.n_df <- data.frame(imdf$x, imdf$y, imRrs.n)
    names(imRrs.n_df)[c(1,2)] <- c('x','y')
    raster.imRrs.n <- rasterFromXYZ(imRrs.n_df)
    crs(raster.imRrs.n) <- im@crs
    writeRaster(raster.imRrs.n, fn_imRrs.n, format='ENVI',overwrite=TRUE)
    message(paste0("Normalized Rrs map was generated, named: ", fn_imRrs.n))
  }

  # Apply FCM-m to image dataframe
  v <- res$res.FCM$v
  m <- res$res.FCM$m
  message('Apply fcm to image dataframe ......')
  res.FCM <- apply_FCM_m(Rrs=imRrs.n, wavelength=wv, Rrs_clusters=v,
                         default.cluster=F)
  res.im <- list()
  res.im$u <- res.FCM$u %>% as.data.frame
  res.im$cluster <- res.FCM$cluster

  # Save true color image
  rgb <- brick(input[[10]],input[[7]],input[[5]],input[[2]])
  rgb_stretch <- stretch(x=rgb, minv=0, maxv=255)
  rgb_df <- as.data.frame(rgb_stretch,xy=T)
  rgb_df <- data.frame(x=rgb_df$x, y=rgb_df$y,
                       n=rgb_df[,3],r=rgb_df[,4], g=rgb_df[,5],b=rgb_df[,6]) %>% na.omit
  p.truecolor=ggplot(data=rgb_df) +
    geom_raster(aes(x=x,y=y,fill=rgb(r,g,b, maxColorValue=255)), hjust=0.5, vjust=0.5) +
    scale_fill_identity() +
    labs(title=title.name) +
    coord_equal() + theme_map() +
    theme(plot.title=element_text(hjust=0.5))
  print(p.truecolor)

  # Save membership rasters
  sub.memb <- data.frame(imdf$x,imdf$y, res.im$u)
  names(sub.memb) <- c("x","y",paste0('Cluster_',seq(1,res$K)))
  raster.memb <- rasterFromXYZ(sub.memb)
  crs(raster.memb) <- im@crs
  if(output_image){
    writeRaster(raster.memb, fn_memb, format=output_format,overwrite=TRUE)
    message(paste0("Membership map was generated, named: ", fn_memb))
  }
  message("Plotting membership values ......")
  p.memb<- ggplot() +
    geom_raster(data=melt(sub.memb,id=c('x','y')), aes(x=x,y=y,fill=value),
                hjust=0.5, vjust=0.5) +
    scale_fill_viridis(na.value='white',begin=0,end=1,
                       breaks=c(0.0,0.25,0.5,0.75,1.0)) +
    labs(fill='Membership', title=title.name) +
    coord_equal() + theme_map() +
    theme(plot.title=element_text(hjust=0.5),
          legend.position=c(1,0),
          legend.key.width=unit(1,"lines"),
          legend.justification=c(1,0),
          strip.background=element_rect(fill='white',color='white')) +
    facet_wrap(~variable, nrow=2)
  print(p.memb)

  # Save cluster raster
  im.cluster <- data.frame(imdf$x,imdf$y, res.im$cluster)
  names(im.cluster)[c(1,2)] <- c("x","y")
  raster.cluster <- rasterFromXYZ(im.cluster)
  crs(raster.cluster) <- im@crs
  if(output_image){
    writeRaster(raster.cluster, fn_cluster, format=output_format,overwrite=TRUE)
    message(paste0("Cluster map was generated, named: ", fn_cluster))
  }
  message("Plotting clusters ......")
  cp <- RdYlBu(nrow(v))
  cp.sub <- cp[(unique(res.im$cluster)) %>% sort]
  p.cluster <- ggplot() +
    geom_raster(data=im.cluster,aes(x=x,y=y,fill=as.character(res.im.cluster)),
                hjust=0.5, vjust=0.5) +
    scale_fill_manual(values=cp.sub) +
    labs(fill='Cluster', title=title.name) +
    coord_equal() + theme_map() +
    theme(plot.title=element_text(hjust=0.5),
          legend.position='right',
          legend.justification='center',
          legend.key.width=unit(1.5,"lines"))
  print(p.cluster)

  # obtain width and height of ggsave
  width <- raster.cluster@ncols/png_scale+0.3
  height<- raster.cluster@nrows/png_scale

  res.Chla <- NULL
  p.Chla <- NULL
  # Chla estimation and plot Chla map
  if(Chla_est){
    X <- data.frame(Rrs665 = res.FCM$x$Rrs665,
                    Rrs709 = res.FCM$x$Rrs709,
                    Rrs754 = res.FCM$x$Rrs754,
                    M=res.FCM$u)
    res.Chla <- FCM_m_Chla_estimation(Rrs=X[,1:3],U=X[,4:10])
    sub.Chla <- data.frame(imdf$x,imdf$y, res.Chla$conc.Blend)
    names(sub.Chla) <- c("x","y","Chla")
    raster.Chla <- rasterFromXYZ(sub.Chla)
    crs(raster.Chla) <- im@crs
    if(output_image){
      writeRaster(raster.Chla, fn_Chla, format=output_format,overwrite=TRUE)
      message(paste0("Chla concentration map was generated, named: ", fn_Chla))
    }
    message("Plotting Chla concentration ......")
    options(scipen=1000)
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
    if(max(sub.Chla$Chla, na.rm=T) >= 800){
      p.Chla <- p.Chla +
        scale_fill_viridis(na.value='gray', trans='log10',
                           limits=c(1,800),
                           breaks=c(10^seq(0,2.9,0.5)) %>% round(.,2))
    }else{
      p.Chla <- p.Chla +
        scale_fill_viridis(na.value='gray', trans='log10',
                           breaks=c(10^seq(0,2.9,0.5)) %>%
                             .[. < max(sub.Chla$Chla, na.rm=T)] %>% round(.,2))
    }
    print(p.Chla)
    ggsave(paste0(fn_Chla, '.png'), plot=p.Chla, device='png', width=width, height=height, units='in')
  }

  # Save ggplot png for membership and cluster patterns
  if(output_resultpng){
    ggsave(paste0(fn_memb,'.png'), plot=p.memb, device='png')
    ggsave(paste0(fn_cluster,'.png'), plot=p.cluster, device='png', width=width, height=height, units='in')
    ggsave(paste0(fn_truecolorpng, '.png'), plot=p.truecolor, device='png', width=width, height=height, units='in')
  }

  message('Saving the results to the list ......')
  message(fn_memb)
  message(fn_cluster)

  # Save results to one list
  result = list()
  result$input <- input
  result$res <- res
  result$raster.memb <- raster.memb
  result$raster.cluster <- raster.cluster
  result$p.memb <- p.memb
  result$p.cluster <- p.cluster
  result$p.truecolor <- p.truecolor
  result$p.Chla <- p.Chla
  result$imdf <- imdf
  result$imRrs.raw <- imRrs.raw
  result$imRrs.n <- imRrs.n
  result$res.FCM <- res.FCM
  result$res.Chla <- res.Chla

  message('Done!')
  return(result)
}

#' @title generate_param
#' @name generate_param
#' @param wl wavelength of subsetting
#' @export
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
