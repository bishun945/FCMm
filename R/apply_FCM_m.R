#' @title Apply FCM_m to new input Rrs
#' @name apply_FCM_m
#' @description
#' Application of FCM_m method for new Rrs data based on default cluster settings
#'   or user-defined clusters (trained by FCM.new).
#'
#' @usage apply_FCM_m(Rrs, wavelength, Rrs_clusters,
#'   stand=FALSE, default.cluster=TRUE, m_used=1.36, option.plot=FALSE)
#'
#' @param Rrs Data.frame, the input Rrs of FCM.
#' @param wavelength Numeric vector, used for applying FCM.
#'   Default use the data from \code{Bi_clusters.rda}
#' @param Rrs_clusters Data.frame, used for applying FCM.
#'   Default use the data from \code{Bi_clusters.rda}
#' @param stand Logical, whether to normalized the Rrs data. Default as \code{FALSE} means do not.
#' @param default.cluster Logical, whether to use the default clusters.
#'   Default use the data from \code{Bi_clusters.rda}
#' @param m_used Number, Used fuzzifier value
#' @param option.plot Logical, whether to plot the result. Default as \code{FALSE}
#'
#' @export
#'
#' @return A \code{list} including several results of function \code{apply_FCM_m()}
#'   \itemize{
#'     \item \strong{x}  The raw input Rrs dataframe with unit sr^-1
#'     \item \strong{x.stand}  The standardized Rrs dataframe, if \code{stand=F}
#'     \item \strong{d}  Distance to each cluster
#'     \item \strong{u}  Membership values
#'     \item \strong{cluster}  Defined by the maximum of membership
#'     \item \strong{quality}  The quality of the cluster results.
#'     \item \strong{m.used}  The used value of fuzzifier(m)
#'     \item \strong{K}  Cluster number
#'     \item \strong{p.group}  A ggplot list for plotting the cluster result
#'     \item \strong{p.group.facet} \code{p.group} with facet to see each cluster resutls
#'       more clearly
#'     \item \strong{dt.plot}  Dataframe used for ggplot
#'   }
#' 
#' @examples
#' \dontrun{
#' library(FCMm)
#' data("WaterSpec35")
#' data("Bi_clusters")
#' Rrs <- WaterSpec35[,3:17]
#' result <- apply_FCM_m(Rrs=Rrs, option.plot=TRUE)
#' }
#' 
#' @references 
#' Bi S, Li Y, Xu J, et al. Optical classification of inland waters based on
#'   an improved Fuzzy C-Means method[J]. Optics Express, 2019, 27(24): 34838-34856.
#' 
#' @family Fuzzy cluster functions


apply_FCM_m <- function(Rrs, wavelength, Rrs_clusters,
                        stand=FALSE, default.cluster=TRUE, m_used=1.36,
                        option.plot=FALSE){

  if(default.cluster){
    v <- Rrs_clusters.default
    wavelength <- wavelength.default
  }else{
    v <- Rrs_clusters
  }

  if(!is.data.frame(Rrs))
    stop('The input param Rrs should be a data.frame!')
  if(!is.data.frame(v))
    stop('The input param Rrs_clusters should be a data.frame!')
  if(!is.numeric(wavelength))
    stop('The input param wavelength should be numeric!')
  
  if(length(wavelength) != ncol(v))
    stop("The number of wavelength is different from the col of Rrs_clusters, please check them!")
  if(length(wavelength) != ncol(Rrs))
    stop("The col of input Rrs dataframe is different from the number of wavelength, please check!")

  k <- nrow(v)
  x <- Rrs
  Area_x <- .trapz(wavelength, x)
  x.stand <- x
  for(i in 1:ncol(x)){x.stand[,i] <- x[,i]/Area_x}

  Area_v <- .trapz(wavelength, v)
  v.stand <- v
  for(i in 1:ncol(v)){v.stand[,i] <- v[,i]/Area_v}

  if(stand==FALSE){
    x_ <- x.stand
    v_ <- v.stand
  }else{
    x_ <- x
    v_ <- v
  }

  # Build distance and membership matrix
  d <- matrix(NA, ncol=nrow(v_), nrow=nrow(x_))
  for(j in 1:ncol(d)){
    v_tmp <- rep(as.matrix(v_[j,]),nrow(x_)) %>% matrix(.,ncol=ncol(v_), byrow=T)
    x_tmp <- x_ %>% as.matrix %>% as.data.frame
    row.names(x_tmp) <- row.names(v_tmp)
    d_tmp <- (x_tmp - v_tmp)^2 %>% apply(.,1,sum) %>% as.matrix
    d[,j] <- d_tmp
  }

  m <- m_used
  u <- 1/d^(2/(m-1))/(apply(1/d^(2/(m-1)),1,sum))
  cluster <- apply(u,1,which.max)

  quality <- rep('Dubious',nrow(x_))
  w <- which( (u %>% apply(.,1,max)) >= (k-1)/k)
  quality[w] <- 'Believable'

  result <- list()
  result$x <- Rrs
  result$x.stand <- x.stand
  result$d <- d
  result$u <- u
  result$cluster <- cluster
  result$quality <- quality
  result$m_used <- m_used
  result$K <- k

  if(option.plot){# plot job
    cp <- RdYlBu(k)
    cp.sub <- cp[(unique(cluster)) %>% sort]
    names(Rrs) <- as.character(wavelength)
    dt <- cbind(nm=as.character(seq(1,nrow(Rrs))),Rrs)
    dt$nm %<>% levels(.)[.]
    dt.melt <- melt(dt,id=c('nm'),variable.name='band',value.name='Rrs')
    dt.melt$band %<>% levels(.)[.] %>% as.numeric
    dt2.melt <- cbind(cluster=cluster, Rrs) %>% melt(., id=c('cluster'))
    dt3.melt <- cbind(quality=quality,Rrs) %>% melt(., id=c('quality'))
    dt3.melt$quality %<>% levels(.)[.]
    dt.plot <- cbind(dt.melt,cluster=as.character(dt2.melt$cluster),quality=dt3.melt$quality)
    dt.plot$cluster %<>% levels(.)[.]
    dt.plot$quality %<>% levels(.)[.]
    p.group <- ggplot(data=dt.plot)+
      geom_line(aes(x=band,y=Rrs,color=cluster,group=nm, linetype=quality),
                size=1, alpha=0.8) +
      scale_color_manual(values=cp.sub) +
      scale_linetype_manual(values=c('solid','longdash')) +
      labs(x='Wavelength (nm)', y=expression(Rrs~(sr^-1))) +
      theme_bw() +
      theme(text = element_text(face='bold',size=16), legend.position='right')
    
    result$p.group <- p.group
    result$p.group.facet <- p.group + facet_wrap(~cluster)
    result$dt.plot <- dt.plot
  }else{
    result$p.group <- NULL
    result$p.group.facet <- NULL
    result$dt.plot <- NULL
  }

  return(result)
}
