#' @name FuzzifierDetermination
#' @title Determine the optimized fuzzifier for FCM running
#' @description
#' To determine the optimized fuzzifier value for FCM running.
#' @usage FuzzifierDetermination(x, wv, max.m=10, stand=FALSE, dmetric='sqeuclidean')
#' @param x Data.frame. the input Rrs data
#' @param wv Wavelength of X
#' @param max.m Set max.m as for determination of m.mub. Default as 10
#' @param stand Whether the data set is standardized (T needs not standardization)
#' @param dmetric Distance method. Default as 'sqeuclidean'
#' @return \code{FD} list contains several result by \code{FuzzifierDetermination}:
#'   \itemize{
#'     \item \strong{x} The raw input Rrs dataframe with unit sr^-1
#'     \item \strong{x.stand} The standardized Rrs dataframe, if \code{stand=F}
#'     \item \strong{wv} Wavelength with unit nm
#'     \item \strong{max.m} The maximum fuzzifier of FCM as a restriction
#'     \item \strong{stand} A logic value for choosing whether to use standardization
#'     \item \strong{dmetric} A string value for choosing which distance metric to be used
#'     \item \strong{Area} A numeric vector for trapezoidal integral values
#'     \item \strong{m.ub} The upper boundary of fuzzifier(m) value
#'     \item \strong{m.used} The desired value of fuzzifier(m) value
#'   }
#' @export
#' @examples 
#' \dontrun{
#' library(FCMm)
#' library(tidyverse)
#' data("Nechad2015")
#' w <- Nechad2015 %>% names %>%
#'     str_extract(.,pattern="\\d") %>%
#'     is.na %>% {!.}
#' wv <- w %>% names(Nechad2015)[.] %>%
#'     gsub('X','',.) %>% as.numeric
#' x <- w %>% Nechad2015[,.]
#' names(x) <- wv
#' rm(w)
#' FD <- FuzzifierDetermination(x, wv, stand=F)
#' }
#' @references 
#' \itemize{
#'   \item Bi S, Li Y, Xu J, et al. Optical classification of inland waters based on
#'     an improved Fuzzy C-Means method[J]. Optics Express, 2019, 27(24): 34838-34856.
#'   \item Dembele D. Multi-objective optimization for clustering 3-way gene
#'     expression data[J]. Advances in Data Analysis and Classification, 2008, 2(3):
#'     211-225.
#' }
#' @family Fuzzy cluster functions
#' 
#' @importFrom stats sd

FuzzifierDetermination <- function(x, wv, max.m=10, stand=FALSE, dmetric="sqeuclidean"){
  
  # about x
  if(missing(x))
    stop("Missing input data set!")
  if(is.null(x))
    stop("Data set is null")
  if((is.data.frame(x))||(is.vector(x)))
    x <- as.matrix(x)
  if(!is.matrix(x))
    stop("Data set must be a vector, data frame or matrix")
  if(any(is.na(x)))
    stop("Data set should not contain NA values. Please remove NAs and try again")
  if(!is.numeric(x))
    stop("Data set must be a numeric vector, data frame or matrix")
  
  # about wv
  if(missing(wv))
    stop("Missing wavelength!")
  if(is.null(wv))
    stop("Wavelength is null")
  if(!is.vector(wv)){
    if(dim(wv)[2] != dim(x)[2]){
      stop("the dimension of wavelength and input data set is different! --A")}
    wv <- t(wv)
  }
  if(dim(t(wv))[2] != dim(x)[2]){
    stop("the dimension of wavelength and input data set is different! --B")}
  
  if(!is.numeric(wv))
    stop("Wavelength must be a numeric vector")
  
  # about max.m
  if(!is.matrix(x))
    stop("max.m must be a numeric")
  if(max.m > 50)
    stop("Are you kidding me? The max.m is too large!")
  
  # stand
  if(!is.logical(stand))
    stop("Stand must be TRUE or FALSE")
  
  # demteric
  if(dmetric != "sqeuclidean")
    stop("Sorry, this function could only support sqeuclidean distance metric!")
  
  Area <- NULL
  x.stand <- x
  if(stand==FALSE){
    Area <- .trapz(wv,x)
    x.stand <- x
    for(i in 1:ncol(x)){
      x.stand[,i] = x[,i] / Area
    }
  }
  
  m.ub <- .commub(x.stand, max.m)
  
  if(m.ub >= 10){
    warning("The upon boundary fuzzifier m is larger than 10!")
    stop("Maybe this Rrs data set is not good for clustering, Please check outliers and run again")
  }
  
  m.used <- 1+m.ub/10
  
  m <- m.used
  
  result = list()
  result$x <- x
  result$x.stand <- x.stand # x.stand is raw matrix if stand.flag is TRUE
  result$wv <- wv
  result$max.m <- max.m
  result$stand <- stand
  result$dmetric <- dmetric
  result$Area <- Area
  result$m.ub <- m.ub
  result$m.used <- m.used
  
  return(result)
  
}

.trapz <- function(wv,x){
  Area <-  as.matrix(x[,1])
  for(i in 1:dim(x)[1]){
    y <- x[i,]
    len <- length(y)
    Area[i] <- 0.5 * sum(diff(wv,1)*(y[1:(len-1)]+y[2:len]))
  }
  return(as.matrix(Area))
}

.commub <- function(x,max.m){
  dismt <- .compdismt(x)
  ii <- 1
  Ym <- NULL
  m.arr <- seq(from=1.1,to=max.m,by=0.1)
  for(m in m.arr){
    Ym[ii] <- abs(sd(dismt^(1/(m-1)))/mean(dismt^(1/(m-1))) - 0.03*dim(x)[2])/(0.03*dim(x)[2])
    ii <- ii + 1
  }
  if(which.min(Ym)==length(Ym))
    warning("The calculated m.ub is equal to max.m. Please set a larger max.m")
  m.ub <- m.arr[which.min(Ym)]
  return(m.ub)
}

.compdismt <-function(x){
  x <- as.matrix(x)
  n <- dim(x)[1]
  arr <- array(dim=c(n^2-n,1))
  ind <- 1
  for(i in 1:n){
    for(j in 1:n){
      if(i != j){
        a <- as.vector(x[i,])
        b <- as.vector(x[j,])
        arr[ind] <- sum(t(a-b) * (a-b))
        # arr[ind] <- .compdist(a,b)
        ind <- ind+1
      }
    }
  }
  return(arr)
}



#' @name FCM.new
#' @title Running the improved FCM by optimizing fuzzifier parameter
#' @description
#' An improved version of FCM for water spectra data sets.
#'
#' @usage FCM.new(FDlist, K, plot.jitter=TRUE, fast.mode=FALSE, stand=FALSE)
#'
#' @param FDlist A \code{list} from function \code{\link{FuzzifierDetermination}}
#' @param K Number, cluster number
#' @param plot.jitter Logical, choose to plot jitter plot using \code{ggplot2} functions
#' @param fast.mode Logical, \code{FALSE} (default)
#' @param stand Logical, choose to do normalization of input data
#'
#' @return A \code{list} of FCM:
#'   \itemize{
#'     \item \strong{FD} The return list by function \code{\link{FuzzifierDetermination}}
#'     \item \strong{res.FCM} The optimized FCM result generated by functions in package \code{ppclust}
#'     \item \strong{K} Cluster number
#'     \item \strong{plot.jitter} A logical value for the option of doing jitter plot by package \code{ggplot2}
#'     \item \strong{fast.mode} A logical value for choosing whether to use fast mode
#'   }
#'
#' @export
#' @examples 
#' \dontrun{
#' library(FCMm)
#' library(tidyverse)
#' data("Nechad2015")
#' w <- Nechad2015 %>% names %>%
#'     str_extract(.,pattern="\\d") %>%
#'     is.na %>% {!.}
#' wv <- w %>% names(Nechad2015)[.] %>%
#'     gsub('X','',.) %>% as.numeric
#' x <- w %>% Nechad2015[,.]
#' names(x) <- wv
#' rm(w)
#' FD <- FuzzifierDetermination(x, wv, stand=F)
#' nb <- 4
#' set.seed(54321)
#' result <- FCM.new(FD, nb)
#' print(result$p.jitter)
#' }
#' 
#' @references 
#' \itemize{
#'   \item Bi S, Li Y, Xu J, et al. Optical classification of inland waters based on
#'     an improved Fuzzy C-Means method[J]. Optics Express, 2019, 27(24): 34838-34856.
#'   \item Dembele D. Multi-objective optimization for clustering 3-way gene
#'     expression data[J]. Advances in Data Analysis and Classification, 2008, 2(3):
#'     211-225.
#' }
#'
#' @family Fuzzy cluster functions
#' 
#' @import ggplot2
#' @importFrom inaparc kmpp imembrand
#' @importFrom ppclust fcm
#' @importFrom reshape2 melt
#' 
FCM.new <- function(FDlist, K, plot.jitter=TRUE, fast.mode=FALSE,stand=FALSE){
  # K is the cluster number
  if(missing(FDlist)){
    warning("Have your run the function `FuzzifierDetermination` sucessfully?")
    stop("Missing FDlist!")
  }
  if(missing(K))
    stop("Missing cluster number K!")
  if(is.null(K))
    stop("Cluster number K is null")
  if(!is.numeric(K))
    stop("Cluster number must be a numeric value!")
  
  if(stand==FALSE){
    x <- FDlist$x.stand
  }else if(stand==TRUE){
    x <- FDlist$x
  }
  
  v <- kmpp(x, k=K)$v
  u <- imembrand(nrow(x), k=K)$u
  
  if(fast.mode==FALSE){
    res <- fcm(x, centers=v, memberships=u, m=FDlist$m.used, stand=F)
  }else if(fast.mode==TRUE){
    res <- fcm(x, centers=v, memberships=u, m=FDlist$m.used, stand=F,con.val=1e-4)
    message("The fast mode is used. Pay attention to the convergence of the obj. fun. and results!")
  }
  res$u <- .cal.new.membership(res)
  
  p.jitter <- NULL
  if(plot.jitter)
    p.jitter <- .plot.jitter(res)
  
  result = list()
  result$FD <- FDlist
  result$res.FCM <- res
  result$p.jitter <- p.jitter
  result$K <- K
  result$plot.jitter <- plot.jitter
  result$fast.mode <- fast.mode
  
  return(result)
}



#' @importFrom ppclust is.ppclust
.cal.new.membership <- function(res){
  if(!is.ppclust(res))
    stop("The input value is not a ppclust list!")
  return(1/res$d^(2/(res$m-1))/(apply(1/res$d^(2/(res$m-1)),1,sum)))
}



#' @importFrom ppclust is.ppclust
#' @importFrom reshape2 melt
.plot.jitter <- function(res){
  if(!is.ppclust(res))
    stop("The input value is not a ppclust list!")
  mem <- res$u
  m   <- res$m
  k   <- res$k
  SampleID <- rownames(mem)
  name<- res$algorithm
  tmp <- data.frame(mem)
  tmp <- cbind(SampleID,tmp)
  names(tmp)=c("SampleID",as.character(seq(from=1,to=k)))
  tmp.n <- tmp[,-1]
  ind<-sort(apply(tmp[,-1],2,sum),decreasing = T,index.return=T)$ix
  tmp.n2 <- tmp.n[,ind]
  tmp <- cbind(tmp$SampleID,tmp.n2)
  tmp.p <- melt(tmp,id.vars=1)
  
  p <- ggplot(data=tmp.p) +
    geom_jitter(aes(x=variable,y=value),alpha=I(1/10)) +
    ylim(c(0,1)) +
    theme(plot.title=element_text(hjust=0.5)) +
    labs(title=paste0(name," m=",as.character(m)),x="Cluster",y="Membership degree")
  return(p)
}



#' @name apply_FCM_m
#' @title Apply FCM_m to new input Rrs
#' @description
#' Application of FCM_m method for new Rrs data based on default cluster settings
#'   or user-defined clusters (trained by FCM.new).
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
#' @param quality_check Logical, quality chech option (default as \code{FALSE})
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
#'     \item \strong{Area} Spectral intergration of each sample
#'     \item \strong{cluster}  Defined by the maximum of membership
#'     \item \strong{quality}  The quality of the cluster results.
#'     \item \strong{m.used}  The used value of fuzzifier(m)
#'     \item \strong{K}  Cluster number
#'     \item \strong{p.group}  A ggplot list for plotting the cluster result
#'     \item \strong{p.group.facet} \code{p.group} with facet to see each cluster results
#'       more clearly
#'     \item \strong{dt.melt}  Dataframe used for ggplot
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
#' 
#' @import ggplot2
#' @importFrom magrittr %>% %<>%
#' @importFrom reshape2 melt
apply_FCM_m <- function(Rrs, wavelength = NULL, Rrs_clusters = NULL,
                        m_used = 1.36,
                        stand = FALSE,
                        default.cluster = TRUE,
                        quality_check = FALSE,
                        option.plot = FALSE){
  
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
  
  v <- as.matrix(v)
  
  k <- nrow(v)
  x <- as.matrix(Rrs)
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
  d <- matrix(ncol=nrow(v_), nrow=nrow(x_))
  for(j in 1:ncol(d)){
    v_tmp <- rep(as.matrix(v_[j,]),nrow(x_)) %>% matrix(.,ncol=ncol(v_), byrow=T)
    x_tmp <- x_ %>% as.matrix %>% as.data.frame
    row.names(x_tmp) <- row.names(v_tmp)
    d_tmp <- (x_tmp - v_tmp)^2 %>% apply(.,1,sum) %>% as.matrix
    d[,j] <- d_tmp
  }
  
  m <- m_used
  u <- 1/d^(2/(m-1))/(apply(1/d^(2/(m-1)),1,sum))
  cluster <- as.numeric(apply(u,1,which.max))
  
  if(quality_check == TRUE){
    quality <- rep('Dubious',nrow(x_))
    w <- which( (u %>% apply(.,1,max)) >= (k-1)/k)
    quality[w] <- 'Believable'
  }else{
    quality <- NULL
  }
  
  result <- list()
  result$x <- Rrs
  result$x.stand <- x.stand
  result$d <- d
  result$u <- u
  result$Area <- Area_x
  result$cluster <- cluster
  result$quality <- quality
  result$m_used <- m_used
  result$K <- k
  
  if(option.plot){
    cp <- RdYlBu(k)
    cp.sub <- cp[(unique(cluster)) %>% sort]
    names(Rrs) <- as.character(wavelength)
    dt <- cbind(nm=rownames(Rrs), Rrs, cluster = as.character(cluster))
    dt.melt <- melt(dt, id=c("nm","cluster"), variable.name='band', value.name='Rrs')
    dt.melt <- level_to_variable(dt.melt)
    dt.melt$band %<>% as.numeric
    p.group <- ggplot(data=dt.melt)+
      geom_line(aes(x=band,y=Rrs,color=cluster,group=nm),
                size=1, alpha=0.8) +
      scale_color_manual(values=cp.sub) +
      labs(x='Wavelength (nm)', y=expression(bold(Rrs~(sr^-1))), color = 'Cluster') +
      theme_bw() +
      theme(text = element_text(face='bold',size=16), legend.position='right')
    
    result$p.group <- p.group
    result$p.group.facet <- p.group + 
      facet_wrap(~cluster) + 
      theme(strip.text = element_blank())
    result$dt.melt <- dt.melt
  }else{
    result$p.group <- NULL
    result$p.group.facet <- NULL
    result$dt.melt <- NULL
  }
  
  return(result)
}



#' @name plot_spec
#' @title Plot the result of FCM_m
#' @description
#' Spectra plot for a FCM.new result
#'
#' @usage plot_spec(res, show.stand=FALSE, HABc=NULL)
#' @param res The result of function `FCM.new`
#' @param show.stand Logical option to show the standardized spectra (default as False)
#' @param HABc Numeric value to specify which cluster owning larger magnitude
#' 
#' @return A list of ggplot lists:
#'   \itemize{
#'     \item \strong{p.all.spec} Spectra of all data.
#'     \item \strong{p.cluster.spec} Spectra of clusters.
#'     \item \strong{p.group.spec.1} Spectra of groups in mode 1
#'     \item \strong{p.group.spec.2} Spectra of groups in mode 2
#'     \item \strong{res} Result list of \code{FCM.new}
#'     \item \strong{show.stand} Logical value whether use standardized y-axis
#'   }
#'
#' @export plot_spec
#' @method plot spec
#' 
#' @examples 
#' \dontrun{p.spec <- plot_spec(result, show.stand=FALSE, HABc=NULL)}
#' 
#' @references
#' Bi S, Li Y, Xu J, et al. Optical classification of inland waters based on
#'   an improved Fuzzy C-Means method[J]. Optics Express, 2019, 27(24): 34838-34856.
#' 
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom magrittr %>% %<>%
plot_spec <- function(res, show.stand=FALSE, HABc=NULL){
  
  p.all.spec <- .plot.all.spec(res,show.stand=show.stand)
  p.cluster.spec <- .plot.cluster.spec(res,show.stand=show.stand)
  p.group.spec.1 <- .plot.group.spec(res,show.stand=show.stand,mem.crisp=T, HABc=HABc)
  p.group.spec.2 <- .plot.group.spec(res,show.stand=show.stand,mem.crisp=F, HABc=HABc)
  
  result = list()
  result$p.all.spec <- p.all.spec
  result$p.cluster.spec <- p.cluster.spec
  result$p.group.spec.1 <- p.group.spec.1
  result$p.group.spec.2 <- p.group.spec.2
  result$res <- res
  result$show.stand = show.stand
  
  return(result)
}



#' @title plot_spec_from_df
#' @name plot_spec_from_df
#' @description
#' Spectra plot for a dataframe given numeric colnames
#'
#' @usage plot_spec_from_df(df)
#' @param df Data.frame for Rrs plotting
#' @note The colnames of input `df` should be numericalize (i.e. names(df) %>% as.numeric)
#' @return A ggplot list
#'
#' @export
#' @family Fuzzy cluster functions
#' 
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom magrittr %>% %<>%
plot_spec_from_df <- function(df){
  if(names(df) %>% as.numeric %>% anyNA)
    stop('The names of input dataframe cannot be converted to numeric wavelength!')
  tmp.spec <- cbind(nm=seq(1,nrow(df)), df) %>% as.data.frame
  names(tmp.spec)[1] <- "name"
  tmp.spec$name %<>% as.character
  spectra.p <- melt(tmp.spec, id=1,variable.name='band',value.name='value')
  spectra.p$band %<>% levels(.)[.] %>% as.numeric
  result <- ggplot(data=spectra.p,
                   aes(x=band,y=value,group=name,color=name)) +
    geom_path(alpha=0.5) +
    scale_colour_viridis_d() +
    theme_bw() +
    theme(legend.position='none')
  return(result)
}



.plot.all.spec <- function(res, show.stand){
  wavelength <- res$FD$wv %>% as.character
  if(show.stand==FALSE){
    tmp.spc <- res$FD$x
    y.lable <- "Raw data"
  }else if(show.stand==TRUE){
    tmp.spc <- res$FD$x.stand
    y.lable <- "Normalized data"
  }
  tmp.spec <- cbind(nm=seq(1,nrow(tmp.spc)), tmp.spc) %>% as.data.frame
  names(tmp.spec) <- c("name",wavelength)
  tmp.spec$name %<>% as.character
  spectra.p <- melt(tmp.spec,id='name',variable.name='band',value.name='value')
  result <- ggplot() +
    geom_path(data=spectra.p,
              aes(x=band,y=value,group=name,color=name),alpha=0.3) +
    scale_colour_viridis_d() +
    labs(title=paste0("All spectra of dataset N=",as.character(dim(res$FD$x)[1])),
         x="Wavelength", y=y.lable) +
    theme_bw() +
    theme(legend.position="none",plot.title=element_text(hjust=0.5))
  return(result)
}



.plot.cluster.spec <- function(res, show.stand){
  wavelength <- res$FD$wv
  p <- length(wavelength)
  if(show.stand==FALSE){
    Rrs.raw <- (as.matrix(res$FD$x))
    y.lable <- expression(R[rs](sr^-1))
  }else if(show.stand==TRUE){
    Rrs.raw <- (as.matrix(res$FD$x.stand))
    y.lable <- "Normalized Rrs"
  }
  
  Rrs.cluster <- aggregate(Rrs.raw, list(res$res.FCM$cluster), mean)
  names(Rrs.cluster)[1] <- 'name'
  tmp.spec <- melt(Rrs.cluster, id='name', variable.name='band', value.name='value')
  tmp.spec[,1] %<>% as.character
  tmp.spec[,2] %<>% levels(.)[.] %>% as.numeric
  result <- ggplot() +
    geom_path(data=tmp.spec,aes(x=band,y=value,group=name,color=name),size=1.1) +
    scale_color_manual(values=RdYlBu(res$K)) +
    labs(title=NULL,x="Wavelength (nm)", y=y.lable,color="Cluster") +
    theme_bw() +
    theme(legend.position='right',
          legend.key.height=unit(0.8,'lines'),
          text=element_text(size=15))
  result
  return(result)
}



.plot.group.spec <- function(res, show.stand, HABc=NULL, mem.crisp=TRUE, mem.threshold=0.8){
  # HABc could be user-defined if checked by function .plot.cluster.spec
  
  # mem.threshold is a number that defines the color for different membership
  # the default mem.threshold is set as 0.5
  if(is.numeric(mem.threshold)){
    if(mem.threshold <= 0 || mem.threshold >= 1)
      stop("The mem.threshold number should be set in (0,1)")
  }
  if(is.null(mem.threshold)){
    mem.threshold <- 0.8
  }
  
  wavelength <- res$FD$wv
  dimension <- length(wavelength)
  if(show.stand==FALSE){
    Rrs.raw <- (as.matrix(res$FD$x))
    y.lable <- expression(R[rs](sr^-1))
  }else if(show.stand==TRUE){
    Rrs.raw <- (as.matrix(res$FD$x.stand))
    y.lable <- "Normalized Rrs"
  }
  
  Rrs.cluster <- aggregate(Rrs.raw, list(res$res.FCM$cluster), mean)
  Rrs.cluster2 <- t(Rrs.cluster[,-c(1)])
  Rrs.cluster2 <- data.frame(wavelength,Rrs.cluster2)
  p.raw <- list()
  k <- res$K
  mem.matrix <- res$res.FCM$u
  data.raw=NULL
  data.raw$SampleID <- as.character(seq(1:(dim(Rrs.raw)[1])))
  u <- res$res.FCM$u
  
  # define ymax
  if(is.null(HABc)){
    HABc <- -1
    ymax.a <- max(Rrs.raw[which(res$res.FCM$cluster!=HABc),])*1.2
  }else if(!is.null(HABc)){
    ymax.a <- max(Rrs.raw[which(res$res.FCM$cluster!=HABc),])*1.2
    ymax.b <- max(Rrs.raw[which(res$res.FCM$cluster==HABc),])*1.2
  }
  
  ymin <- min(Rrs.raw)*0.8
  
  # start plot
  if(mem.crisp==TRUE){
    
    for(i in 1:k){
      w <- as.numeric(which(res$res.FCM$cluster == i & (u[,i] >= mem.threshold)))
      (length(w))
      N <- table(res$res.FCM$cluster)[[i]]
      ymax <- ymax.a
      if(i == HABc){ymax <- ymax.b}
      tmp.spec <- t(as.matrix(Rrs.raw[w,]))
      tmp <- data.frame(wavelength,tmp.spec)
      names(tmp) <- c("wv",as.character(data.raw$SampleID[w]))
      tmp.p <- melt(tmp,id.vars=1) # this is the important code for plot spectra
      cv1 <- rep("midnightblue",times=length(w))
      p <- ggplot() + geom_line(aes(x=wv,y=value,color=variable),data=tmp.p,alpha=0.8,size=0.8) +
        scale_color_manual(values=cv1) +
        theme(legend.position="none",plot.title=element_text(hjust=0.5),text=element_text(size=15)) +
        labs(title=paste0("Cluster ",as.character(i)," N=",as.character(N)),
             x="Wavelength (nm)", y=y.lable) +
        ylim(ymin, ymax)
      w <- as.numeric(which(res$res.FCM$cluster == i & u[,i] < mem.threshold))
      if(length(w)<=1){
        p1 <- p
      }else if(length(w)>1){
        tmp.spec <- t(as.matrix(Rrs.raw[w,]))
        tmp <- data.frame(wavelength,tmp.spec)
        names(tmp) <- c("wv",as.character(data.raw$SampleID[w]))
        tmp.p <- melt(tmp,id.vars=1) # this is the important code for plot spectra
        cv2 <- rep("seashell4",times=length(w))
        p1 <- p + geom_line(aes(x=wv,y=value,color=variable),data=tmp.p,alpha=0.5,size=0.8) +
          scale_color_manual(values = c(cv1,cv2))
      }
      g <- names(Rrs.cluster2)
      pp <- eval(substitute(p1 + geom_line(data=Rrs.cluster2,aes(x=wavelength,y=get(g[i+1])),color="Red",size=1.5),
                            list(i=i)))
      # print(pp)
      p.raw[[i]] <- pp
    }
    
  }else if(mem.crisp==FALSE){
    
    for(i in 1:k){
      w <- which(res$res.FCM$cluster == i)
      N <- table(res$res.FCM$cluster)[[i]]
      ymax <- ymax.a
      if(i == HABc){ymax <- ymax.b}
      tmp.spec <- t(as.matrix(Rrs.raw[w,]))
      tmp <- data.frame(wavelength,tmp.spec)
      names(tmp) <- c("wv",paste0("X",as.character(seq(1:length(w)))))
      tmp.p <- melt(tmp,id.vars=1) # this is the important code for plot spectra
      # cv1 <- rep("seashell4",times=N)
      # print(table(res$res.FCM$cluster)[[i]]==length(unique(tmp.p$variable)))
      p1 <- NULL
      p1 <- ggplot() + geom_line(aes(x=wv,y=value,color=variable),data=tmp.p,alpha=0.8,size=0.8) +
        scale_color_manual(values=rep("seashell4",times=N)) +
        theme(legend.position="none",plot.title=element_text(hjust=0.5),text=element_text(size=15)) +
        labs(title=paste0("Cluster ",as.character(i)," N=",as.character(N)),
             x="Wavelength (nm)", y=y.lable) + ylim(ymin, ymax)
      g <- names(Rrs.cluster2)
      pp <- eval(substitute(p1 + geom_line(data=Rrs.cluster2,aes(x=wavelength,y=get(g[i+1])),color="Black",size=1.5),
                            list(i=i)))
      p.raw[[i]] <- pp
      # print(pp)
    }
  }
  
  return(p.raw)
}
