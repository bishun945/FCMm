#' @name FuzzifierDetermination
#' @title Determine the optimized fuzzifier for Fuzzy Cluster Method (FCM) running
#' @description
#' To determine the optimized fuzzifier value for Fuzzy Cluster Method (FCM)running.
#' @usage FuzzifierDetermination(x, wv, max.m=10, do.stand=TRUE, stand=NULL, dmetric="sqeuclidean")
#' @param x Data.frame. the input Rrs data
#' @param wv Wavelength of X
#' @param max.m Set max.m as for determination of m.mub. Default as 10
#' @param do.stand Whether run standarization for the input data set. Default as \code{TRUE}.
#'   Note that the do.stand only be used for calculating the fuzzifier as both raw and standarzied 
#'   spectra will be restored in the return value. This is benifit to the spectra plotting.
#' @param stand Deprecated; Now \code{stand = !do.stand}
#' @param dmetric Distance method. Default as 'sqeuclidean'
#' @return \code{FD} list contains several result by \code{FuzzifierDetermination}:
#'   \itemize{
#'     \item \strong{x} The raw input Rrs dataframe with unit sr^-1
#'     \item \strong{x.stand} The standardized Rrs dataframe, if \code{do.stand=TRUE}
#'     \item \strong{wv} Wavelength with unit nm
#'     \item \strong{max.m} The maximum fuzzifier of FCM as a restriction
#'     \item \strong{do.stand} A logic value for whether we standardized the input data
#'     \item \strong{dmetric} A string value for choosing which distance metric to be used
#'     \item \strong{Area} A numeric vector for trapezoidal integral values
#'     \item \strong{m.ub} The upper boundary of fuzzifier(m) value
#'     \item \strong{m.used} The desired value of fuzzifier(m) value
#'   }
#' @export
#' @examples 
#' 
#' library(FCMm) 
#' library(magrittr)
#' data("Nechad2015")
#' x <- Nechad2015[,3:11]
#' wv <- gsub("X","",names(x)) %>% as.numeric
#' set.seed(1234)
#' w <- sample(1:nrow(x), 100) 
#' x <- x[w, ]
#' names(x) <- wv
#' set.seed(1234)
#' FD <- FuzzifierDetermination(x, wv, stand=FALSE)
#' 
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

FuzzifierDetermination <- function(x, wv, max.m=10, do.stand=TRUE,
                                   stand=NULL, dmetric="sqeuclidean"){
  
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
  
  # do.stand
  if(!is.logical(do.stand))
    stop("do.stand must be TRUE or FALSE")
  
  # demteric
  if(dmetric != "sqeuclidean")
    stop("Sorry, this function could only support sqeuclidean distance metric!")
  
  # convert tag stand to do.stand
  # do.stand == TRUE means we gonna do standarization for input data
  if(is.logical(stand)){
    do.stand = !stand
    warning("Parameter stand has been deprecated. Please use do.stand!")
  }
  
  x.stand <- x
  # Area <- trapz(wv,x)
  Area <- trapz2(x)
  for(i in 1:ncol(x)){
    x.stand[,i] = x[,i] / Area
  }

  if(do.stand==TRUE){
    m.ub <- .commub(x.stand, max.m)
  }else{
    m.ub <- .commub(x, max.m)
  }
  
  if(m.ub >= 10){
    warning("The upon boundary fuzzifier m is larger than 10!")
    stop("Maybe this Rrs data set is not good for clustering, Please check outliers and run again")
  }
  
  m.used <- 1+m.ub/10
  
  m <- m.used
  
  result = list()
  result$x <- x
  result$x.stand <- x.stand
  result$wv <- wv
  result$max.m <- max.m
  result$do.stand <- do.stand
  result$dmetric <- dmetric
  result$Area <- Area
  result$m.ub <- m.ub
  result$m.used <- m.used
  
  return(result)
  
}

#' @name fDN
#' @title Determining m value by the method of Schwammle and Jensen (2010)
#' @param x data.frame x
#' @return m value
#' @references 
#'   Schwammle V, Jensen O N. A simple and fast method to determine the parameters for 
#'     fuzzy c–means cluster analysis[J]. Bioinformatics, 2010, 26(22): 2841-2848.
#' @noRd
fDN <- function(x){
  stopifnot({
    is.data.frame(x)
    is.matrix(x)
  })
  N = nrow(x)
  D = ncol(x)
  1 + (1418 / N + 22.05)*D^(-2) + (12.33 / N + 0.243) * D ^ (-0.0406*log(N) - 0.1134)
}

#' @name FuzzifierDetermination2
#' @title The version 2 of FuzzifierDetermination
#' @param x data.frame x.
#' @param iter.max the max of iteration number. defualt as 300.
#' @param m_ini_method the initialization m value. 2 (default) or "Schwammle2010". 
#' @return m_ub and m_used
#' @export
#' @importFrom stats dist approx
#' @references 
#' \itemize{
#'   \item Bi S, Li Y, Xu J, et al. Optical classification of inland waters based on
#'     an improved Fuzzy C-Means method[J]. Optics Express, 2019, 27(24): 34838-34856.
#'   \item Dembele D. Multi-objective optimization for clustering 3-way gene
#'     expression data[J]. Advances in Data Analysis and Classification, 2008, 2(3):
#'     211-225.
#'   \item Schwammle V, Jensen O N. A simple and fast method to determine the parameters for 
#'     fuzzy c–means cluster analysis[J]. Bioinformatics, 2010, 26(22): 2841-2848.
#' }
FuzzifierDetermination2 <- function(x, iter.max = 300, m_ini_method = 2){
  
  find_minus <- function(m_ini, m_step){
    for(i in 1:iter.max){
      
      if(i == 1){
        m = round(m_ini, 2)
        Y = dist(x)^(2/(m-1))
      }else{
        Y = dist(x)^(2/(m-1))
      }
      jud <- sd(Y)/mean(Y) - 0.03*ncol(x)
      if(jud > 0){
        m = m + m_step
      }else{
        m_cen <- m
        break
      }
      if(i == iter.max){
        m_cen <- m
      }
    }
    return(list(m_cen=m_cen, m_step=m_step))
  }
  
  for(j in 1:iter.max){
    if(j == 1){
      if(m_ini_method == 2){
        m_cen_res <- find_minus(2, 0.5)
      }else if(m_ini_method == "Schwammle2010"){
        m_cen_res <- find_minus((fDN(x)-1) * 10, 0.2)
      }
      m_step = m_cen_res$m_step/2
      m_ini  = m_cen_res$m_cen - m_cen_res$m_step
    }else{
      m_cen_res <- find_minus(m_ini, m_step)
      m_step = m_cen_res$m_step/2
      m_ini  = m_cen_res$m_cen - m_cen_res$m_step
    }
    if(round(m_step, 2) == 0.01){
      break
    }
  }
  
  m_cen = m_cen_res$m_cen %>% round(2)
  
  m_range <- seq(m_cen-0.03, m_cen+0.03, 0.01)
  jud = rep(NA, length(m_range))
  for(i in 1:length(m_range)){
    Y = dist(x)^(2/(m_range[i]-1))
    jud[i] = (sd(Y)/mean(Y) - 0.03*ncol(x))
  }
  
  m_ub = approx(x=jud, y=m_range, xout=0)$y %>% round(2)
  
  if(m_ub >= 10){
    m_used <- 2
  }else{
    m_used <- 1 + m_ub / 10
  }
  
  return(list(m_ub=m_ub, m_used=m_used))
  
}


#' @name trapz
#' @title Trapezoid integral calculation
#' @param wv wavelength
#' @param x data.frame that need to be trapzed with ncol equal to length of \code{wv}.
#'   Note that the colnames of \code{x} should be converted to numbers for \code{trapz}.
#'   Do not name them as "Rrs560", "560nm" or so. Just set as "560".
#' @return A vector presenting the area result.
#' @details In 2020-11-02, I created an improved version of \code{trapz}, say \code{trapz}.
#'   So better to use \code{trapz} as it calculates much faster than before.
#' @export
#' @importFrom stringr str_extract_all str_c
#' @examples 
#' library(FCMm)
#' library(magrittr)
#' library(stringr)
#' data("Nechad2015")
#' w <- Nechad2015 %>% names %>%
#'     str_extract(.,pattern="\\d") %>%
#'     is.na %>% {!.}
#' wv <- w %>% names(Nechad2015)[.] %>%
#'     gsub('X','',.) %>% as.numeric
#' x <- w %>% Nechad2015[,.]
#' names(x) <- wv
#' rm(w)
#' Area <- trapz(wv, x)
#' Area2 <- trapz2(x) # the results are indentical to that from trapz
#' 
trapz <- function(wv,x){
  Area <-  as.matrix(x[,1])
  for(i in 1:dim(x)[1]){
    y <- x[i,]
    len <- length(y)
    Area[i] <- 0.5 * sum(diff(wv,1)*(y[1:(len-1)]+y[2:len]))
  }
  return(as.matrix(Area))
}

#' @rdname trapz
#' @export
trapz2 <- function(x){
  wv <- colnames(x)
  wv <- str_extract_all(wv, "[:digit:]|\\.", simplify = TRUE) %>% 
    apply(., 1, function(x) str_c(x, collapse = ""))
  if(any(is.na(as.numeric(wv)))) stop("Be sure the colname of x can be converted to number!")
  wv <- as.numeric(wv)
  wv_diff_arr <- matrix(data=rep(diff(wv), nrow(x)), ncol=ncol(x)-1, byrow=TRUE)
  wv_diff_arr <- as.data.frame(wv_diff_arr)
  tmp <- wv_diff_arr * (x[,1:(ncol(x)-1)] + x[,2:ncol(x)]) * 0.5
  Area <- apply(tmp, 1, sum)
  return(Area)
}

#' @name .commub
#' @title Compute the upper boundary of m values
#' @param x Input dataframe
#' @param max.m The defined maxmium of m values
#' @return The used m value.
#' @noRd
#' 
.commub <- function(x, max.m){
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

#' @name .compdismt
#' @title Compute distance matrix
#' @param x Input dataframe
#' @return Distance matrix.
#' @noRd
#' 
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
#' @title Running the improved Fuzzy Cluster Method (FCMm) by optimizing the fuzzifier m
#' @description An improved version of Fuzzy Cluster Method (FCM) for water spectra data sets.
#'
#' @param FDlist A \code{list} from function \code{\link{FuzzifierDetermination}}
#' @param K Number, cluster number
#' @param sort.pos The position to sort the cluster number. Default as the end position of the 
#'   input training matrix. It can be \code{NULL} if there is no need to re-sort the cluster index.
#' @param sort.decreasing Default as \code{FALSE} which means the sort are assigned from small values
#'   to large values.
#' @param plot.jitter Logical, choose to plot jitter plot using \code{ggplot2} functions
#' @param fast.mode Logical, \code{FALSE} (default)
#' @param do.stand Whether to use standarized data for FCM. Default as \code{TRUE}
#' @param stand Deprecated; Now \code{stand = !do.stand}
#'
#' @return A \code{list} of FCM:
#'   \itemize{
#'     \item \strong{FD} The return list by function \code{\link{FuzzifierDetermination}}
#'     \item \strong{res.FCM} The optimized FCM result generated by functions in package \code{ppclust}
#'     \item \strong{K} Cluster number
#'     \item \strong{centroids} The list of centroids (both at raw and normalized scale) of the cluster 
#'       result by aggregating each mean.
#'     \item \strong{plot.jitter} A logical value for the option of doing jitter plot by package \code{ggplot2}
#'     \item \strong{fast.mode} A logical value for choosing whether to use fast mode
#'   }
#'
#' @export
#' @examples 
#' library(FCMm) 
#' library(ggplot2) 
#' library(magrittr)
#' library(stringr)
#' data("Nechad2015")
#' x <- Nechad2015[,3:11]
#' wv <- gsub("X","",names(x)) %>% as.numeric
#' set.seed(1234)
#' w <- sample(1:nrow(x), 100) 
#' x <- x[w, ]
#' names(x) <- wv
#' nb = 4 # Obtained from the vignette "Cluster a new dataset by FCMm"
#' set.seed(1234)
#' FD <- FuzzifierDetermination(x, wv, stand=FALSE)
#' result <- FCM.new(FD, nb, fast.mode = TRUE)
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
#' @importFrom magrittr %>%
#' 
FCM.new <- function(FDlist, K, sort.pos = length(FDlist$wv), sort.decreasing = FALSE, 
                    plot.jitter=TRUE, fast.mode=FALSE, do.stand = TRUE, stand=NULL){
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
  
  if(!is.null(sort.pos)){
    if(sort.pos < 1 | sort.pos > length(FDlist$wv))
      stop("The wrong sort.pos! It should be a number between 1 and ncol of input matrix.")
  }
  
  if(is.logical(stand)){
    do.stand = !stand
    warning("Parameter stand has been deprecated. Please use do.stand!")
  }
  
  if(do.stand==TRUE){
    if(anyNA(FDlist$x.stand)){
      stop("NA found in FDlist$x.stand! Make sure FuzzifierDetermination run with do.stand == TRUE")
    }
    x <- FDlist$x.stand
  }else{
    x <- FDlist$x
  }
  
  v <- kmpp(x, k=K)$v
  u <- imembrand(nrow(x), k=K)$u
  
  if(fast.mode==FALSE){
    res <- fcm(x, centers=v, memberships=u, m=FDlist$m.used, stand=FALSE)
  }else if(fast.mode==TRUE){
    res <- fcm(x, centers=v, memberships=u, m=FDlist$m.used, stand=FALSE, con.val=1e-4)
    message("The fast mode is used. Pay attention to the convergence of the obj. fun. and results!")
  }
  res$u <- .cal.new.membership(res)
  
  # resort cluster
  if(!is.null(sort.pos)){
    x_end <- stats::aggregate(res$x[, ncol(res$x)], list(res$cluster), mean)[, 2]
    new_ind <- sort.int(x_end, index.return = TRUE, decreasing = sort.decreasing) %>% .$ix
    cluster_new <- res$cluster
    for(i in 1:length(new_ind)){
      w = which(res$cluster == new_ind[i])
      cluster_new[w] <- i
    }
    u_new <- res$u[, new_ind]
    colnames(u_new) = paste("Cluster", 1:K)
    
    res$cluster = cluster_new
    res$u       = u_new
  }
  
  # plot jitter
  p.jitter <- NULL
  if(plot.jitter)
    p.jitter <- .plot.jitter(res)
  
  centroids_stand <- stats::aggregate(FDlist$x.stand, list(res$cluster), mean)
  centroids_raw   <- stats::aggregate(FDlist$x, list(res$cluster), mean)
  names(centroids_stand)[1] <- names(centroids_raw)[1] <- "cluster"
  centroids = list(centroids_raw = centroids_raw, 
                   centroids_stand = centroids_stand)
  
  result = list()
  result$FD <- FDlist
  result$res.FCM <- res
  result$p.jitter <- p.jitter
  result$K <- K
  result$centroids <- centroids
  result$plot.jitter <- plot.jitter
  result$fast.mode <- fast.mode
  
  return(result)
}



#' @name .cal.new.membership
#' @title Calculate the membership of results from \code{ppclust}
#' @param x Result of \code{ppclust} object.
#' @return The corrected membership values.
#' @noRd
#' 
.cal.new.membership <- function(res){
  if(!is.ppclust(res))
    stop("The input value is not a ppclust list!")
  return(1/res$d^(2/(res$m-1))/(apply(1/res$d^(2/(res$m-1)),1,sum)))
}


#' @name .plot.jitter
#' @title .plot.jitter
#' @param res Result of \code{ppclust} object.
#' @return A ggplot list showing jitter.
#' @noRd
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
    theme_bw() + 
    theme(text=element_text(size=13),
          plot.title=element_text(hjust=0.5)) +
    labs(title=paste0(name," m=",as.character(m)),x="Cluster",y="Membership degree")
  return(p)
}



#' @name apply_FCM_m
#' @title Apply the improved Fuzzy Cluster Method (FCMm) to new data based on the defined centroids
#' @description
#' Application of the improved Fuzzy Cluster Method (FCMm) for new Rrs data based on 
#'   default cluster settings or user-defined clusters (trained by FCM.new).
#'
#' @param Rrs Data.frame, the input Rrs of FCM.
#' @param wavelength Numeric vector, used for applying FCM.
#'   Default use the data from \code{Bi_clusters.rda}
#' @param Rrs_clusters Data.frame, used for applying FCM.
#'   Default use the data from \code{Bi_clusters.rda}
#' @param stand Deprecated; same to \code{do.stand}.
#' @param do.stand Logical, whether to normalized the Rrs data (both for Input and Centroids).
#'   Default as \code{FALSE} means do not.
#' @param default.cluster Logical, whether to use the default clusters.
#'   Default use the data from \code{Bi_clusters.rda}
#' @param m_used Number, Used fuzzifier value
#' @param quality_check Logical, quality chech option (default as \code{FALSE})
#' @param option.plot Logical, whether to plot the result. Default as \code{FALSE}
#' @param color_palette The palette of cluster color. Default as \code{RdYlBu(res$K)} if \code{NULL}.
#'   It could be \code{RdYlBu(k)}, \code{Spectral(k)}, \code{HUE(k)} or other color values
#'   with same length of cluster number. \code{k} means the cluster number.
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
#'     \item \strong{p.group.facet} Recommened! \code{p.group} with facet to see each cluster results
#'       more clearly
#'     \item \strong{dt.melt} Dataframe used for ggplot including \code{dt_plot_x}, \code{dt_plot_v}, and \code{cp}.
#'   }
#' 
#' @examples
#' library(FCMm)
#' data("WaterSpec35")
#' data("Bi_clusters")
#' Rrs <- WaterSpec35[,3:17]
#' result <- apply_FCM_m(Rrs=Rrs, option.plot=TRUE)
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
                        do.stand = FALSE,
                        stand = NULL,
                        default.cluster = TRUE,
                        quality_check = FALSE,
                        option.plot = FALSE,
                        color_palette = NULL){
  
  
  if(is.logical(stand)){
    do.stand = stand
    warning("Parameter stand has been deprecated. Please use do.stand!")
  }
  
  if(is.null(stand)){
    stand = !do.stand
  }
  
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
  
  if(is.null(color_palette)){
    cp = RdYlBu(k)
    names(cp) <- paste("Cluster", 1:k)
  }else{
    if(length(color_palette) != k){
      stop("The length of color_palette shoud be same with cluster number!")
    }
    cp = color_palette
  }
  
  x <- as.matrix(Rrs)
  # Area_x <- trapz(wavelength, x)
  Area_x <- trapz2(x)
  x.stand <- x
  # for(i in 1:ncol(x)){x.stand[,i] <- x[,i]/Area_x}
  x.stand = x / Area_x
  
  # Area_v <- trapz(wavelength, v)
  Area_v <- trapz2(v)
  v.stand <- v
  # for(i in 1:ncol(v)){v.stand[,i] <- v[,i]/Area_v}
  v.stand <- v / Area_v
  
  if(stand==FALSE){
    x_ <- x.stand
    v_ <- v.stand
    y.lab <- "Normalized Rrs"
  }else{
    x_ <- x
    v_ <- v
    y.lab <- expression(Rrs~group("[",sr^-1,"]"))
  }
  
  # Build distance and membership matrix
  d <- matrix(ncol=nrow(v_), nrow=nrow(x_))
  for(j in 1:ncol(d)){
    v_tmp <- rep(as.matrix(v_[j,]),nrow(x_)) %>% matrix(.,ncol=ncol(v_), byrow=TRUE)
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

    # plot preparation of input data
    dt_plot_x <- cbind(stringsAsFactors = FALSE,
                            data.frame(x_) %>% setNames(., wavelength),
                            cluster = cluster %>% as.character,
                            ids = seq(1, nrow(x_)) %>% as.character) %>%
      melt(., id=c("ids","cluster"), factorsAsStrings = FALSE)
    dt_plot_x$variable %<>% levels(.)[.] %>% as.numeric
    dt_plot_x$cluster_f <- factor(paste("Cluster", dt_plot_x$cluster), 
                                  levels=paste("Cluster", 1:nrow(v_)), ordered = TRUE)
    dt_plot_x$cluster <- factor(dt_plot_x$cluster, levels=1:nrow(v_), ordered = TRUE)
    
    # plot preparation of centroids data
    dt_plot_v <- cbind(stringsAsFactors = FALSE,
                       data.frame(v_) %>% setNames(., wavelength),
                       cluster = 1:nrow(v_) %>% as.character,
                       ids = 1:nrow(v_) %>% as.character) %>%
      melt(., id=c("ids","cluster"), factorsAsStrings = FALSE)
    dt_plot_v$variable %<>% levels(.)[.] %>% as.numeric
    dt_plot_v$cluster_f <- factor(paste("Cluster", dt_plot_v$cluster), 
                                  levels=paste("Cluster", 1:nrow(v_)), ordered = TRUE)
    dt_plot_v$cluster <- factor(dt_plot_v$cluster, levels=1:nrow(v_), ordered = TRUE)
    
    # ggplot 
    p.group.facet <- ggplot() + 
      geom_path(data=dt_plot_x, aes(variable, value, group=ids), color="gray", alpha=0.3) + 
      geom_path(data=dt_plot_v, aes(variable, value, group=ids, color=cluster_f), alpha=1.0, size=1) +
      facet_wrap(~cluster_f) + 
      scale_color_manual(values=cp) + 
      labs(x='Wavelength [nm]', y=y.lab) + 
      theme_bw() + 
      theme(text=element_text(size=13),
            legend.position = "none")
    
    p.group <- ggplot() + 
      geom_path(data=dt_plot_x, aes(variable, value, group=ids, color=cluster_f), alpha=0.5) + 
      # geom_path(data=dt_plot_v, aes(variable, value, group=ids, color=cluster_f), alpha=1.0, size=1) +
      scale_color_manual(values=cp) + 
      # facet_wrap(~cluster_f) + 
      labs(x='Wavelength [nm]', y=y.lab, color = "Legend") + 
      theme_bw() + 
      theme(text=element_text(size=13))
    
    result$p.group <- p.group
    result$p.group.facet <- p.group.facet
    result$dt.melt <- list(dt_plot_x = dt_plot_x, dt_plot_v = dt_plot_v, cp = cp)
    
  }else{
    
    result$p.group <- NULL
    result$p.group.facet <- NULL
    result$dt.melt <- NULL
    
  }
  
  return(result)
}



#' @title Plot spectra from dataframe
#' @name plot_spec_from_df
#' @description
#' Spectra plot for a dataframe given numeric colnames
#'
#' @usage plot_spec_from_df(df)
#' @param df Data.frame for Rrs plotting
#' @note The colnames of input `df` should be numericalize.
#' @return The result of \code{plot_spec_from_df} is a \code{ggplot} object showing spectra curve.
#'
#' @export
#' @family Fuzzy cluster functions
#' 
#' @seealso plot_spec_group
#' 
#' @examples 
#' data(WaterSpec35)
#' plot_spec_from_df(WaterSpec35[, -c(1, 2)])
#' # which is equal to 
#' plot_spec_from_df(WaterSpec35)
#' 
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom magrittr %>% %<>%
#' @importFrom stringr str_detect str_extract_all str_c
plot_spec_from_df <- function(df){
  
  wv = colnames(df)
  w  = str_detect(wv, "[:digit:]|\\.")
  if(sum(w) < ncol(df)){
    warning(sprintf(
      "Following columns are ignored as their names cannot be converted to numeric:\n%s",
      paste0(wv[!w], collapse = " ")
    ))
    df = df[, w]
  }
  
  wv = str_extract_all(names(df),"[:digit:]|\\.",simplify = TRUE) %>% 
    apply(., 1, function(x) str_c(x, collapse = "")) %>% as.numeric
  colnames(df) <- wv
  
  tmp.spec <- cbind(nm=seq(1,nrow(df)) %>% as.character, df, stringsAsFactors = FALSE) %>% 
    as.data.frame
  names(tmp.spec)[1] <- "name"
  spectra.p <- melt(tmp.spec, id=1, variable.name='band', value.name='value')
  spectra.p$band %<>% levels(.)[.] %>% as.numeric
  
  result <- ggplot(data=spectra.p,
                   aes(x=band,y=value,group=name,color=name)) +
    geom_path(alpha=0.5) +
    scale_colour_viridis_d(option = "D") +
    theme_bw() +
    theme(legend.position='none')
  # result
  return(result)
}

#' @title Plot spectra from dataframe by group
#' @name plot_spec_group
#' @description 
#'   This function will help you quick plot the spectra if you have a data.frame
#'   that could be used for \link{plot_spec_from_df} and have a \code{group} parameter which indicates
#'   the different groups for spectra in your \code{x}.
#' @param x The input matrix with colnames could be converted to numeric values.
#' @param group The vector indicate the group-belongings for x
#' @param palette Color palette function. The default is \link{RdYlBu}
#' @param facet Whether to plot by \link{facet_wrap}. The default is \code{TRUE}.
#' @param group_num Whether to add sample numbers of each group in strip.text. Also support 
#'   character as the input.
#' @return A ggplot list.
#' @importFrom reshape2 melt
#' @importFrom magrittr %>% %<>%
#' @importFrom stats aggregate
#' @export
#' 
plot_spec_group <- function(x, group, palette = RdYlBu, facet = TRUE, group_num = TRUE){
  wv <- names(x)
  if(any(as.numeric(wv) %>% is.na)) stop("Be sure the colname of x can be converted to number!")
  if(length(group) != nrow(x)) stop("The length of group should be same as the row number of x!")
  x <- cbind(id=1:nrow(x), group=group, x)
  x_melt <- melt(x, id=c("id", "group"), variable.name="wv") %>% level_to_variable()
  x_melt$wv <- as.numeric(x_melt$wv)
  if(all(is.numeric(group))){
    x_melt$group <- as.numeric(x_melt$group)
  }else{
    x_melt$group <- as.character(x_melt$group)
  }
  # add number
  num <- aggregate(x[,1], list(group), length)[,2]
  x_melt$group_f <- factor(x_melt$group, 
                           levels = sort(unique(x_melt$group)), 
                           ordered = TRUE)
  p <- ggplot() + 
    geom_path(data = x_melt, alpha=0.5, 
              aes(x=wv, y=value, group=id, color=group_f)) + 
    scale_color_manual(values=palette(length(unique(x_melt$group_f))))
  if(facet == TRUE){
    if(group_num == TRUE){
      levels(x_melt$group_f) <- sprintf("%s\nN=%s", levels(x_melt$group_f), num)
    }else if(group_num == FALSE){
      levels(x_melt$group_f) <- sprintf("%s", levels(x_melt$group_f))
    }else if(is.character(group_num)){
      levels(x_melt$group_f) <- sprintf(group_num, levels(x_melt$group_f), num)
    }
    p <- ggplot() + 
      geom_path(data = x_melt, alpha=0.5, 
                aes(x=wv, y=value, group=id, color=group_f)) + 
      scale_color_manual(values=palette(length(unique(x_melt$group_f)))) + 
      facet_wrap(~group_f)
  }
  return(p)
}


#' @name plot_spec
#' @title Plot results of the improved Fuzzy Cluster Method (FCMm)
#' @description
#' Spectra plots for a \code{FCM.new} result
#'
#' @param res The result of function \link{FCM.new}
#' @param show.stand Whether to show the spectra on standarized scale (default as \code{NULL} 
#'   and assigned by \code{res$FD$stand}).
#' @param show.ribbon Whether to show the ribbon defined by sd values calculated from samples 
#'   in each cluster (default as \code{FALSE}).
#' @param color_palette The palette of cluster color. Default as \code{RdYlBu(res$K)}.
#'   In \code{FCMm}, it could be \link{RdYlBu}, \link{Spectral}, \link{HUE} or other color values
#'   with same length of cluster number.
#' 
#' @return The result of \code{plot_spec} is a list including:
#'   \itemize{
#'     \item \strong{p.all.spec} Spectra of all data.
#'     \item \strong{p.cluster.spec} Spectra of clusters.
#'     \item \strong{p.group.spec.1} Spectra of groups in mode 1
#'     \item \strong{p.group.spec.2} Spectra of groups in mode 2
#'     \item \strong{res} Result list of \link{FCM.new}
#'     \item \strong{show.stand} Logical value whether use standardized y-axis
#'   }
#'
#' @export plot_spec
#' 
#' @examples 
#' library(FCMm) 
#' library(ggplot2) 
#' library(magrittr)
#' data("Nechad2015")
#' x <- Nechad2015[,3:11]
#' wv <- gsub("X","",names(x)) %>% as.numeric
#' set.seed(1234) # Set this seed so that you can re-produce them
#' w <- sample(1:nrow(x), 100)
#' x <- x[w, ]
#' names(x) <- wv
#' nb = 4 # Obtained from the vignette "Cluster a new dataset by FCMm"
#' set.seed(1234)
#' FD <- FuzzifierDetermination(x, wv, stand=FALSE)
#' result <- FCM.new(FD, nb, fast.mode = TRUE)
#' p.spec <- plot_spec(result, show.stand = TRUE)
#' print(p.spec$p.cluster.spec)
#' 
#' @references
#' Bi S, Li Y, Xu J, et al. Optical classification of inland waters based on
#'   an improved Fuzzy C-Means method[J]. Optics Express, 2019, 27(24): 34838-34856.
#' 
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom magrittr %>% %<>%
plot_spec <- function(res, 
                      show.stand = NULL, show.ribbon = FALSE, 
                      color_palette = RdYlBu(res$K)){
  
  if(is.null(show.stand)) show.stand = !res$FD$do.stand
  
  if(length(color_palette) != res$K)
    stop("The length of color_palette shoud be same with cluster number!")
  
  p.all.spec <- .plot.all.spec(res, show.stand = show.stand)
  p.cluster.spec <- .plot.cluster.spec(res, 
                                       show.stand = show.stand, show.ribbon = show.ribbon,
                                       color_palette = color_palette)
  p.group.spec.1 <- .plot.group.spec(res, 
                                     show.stand = show.stand, mem.crisp = TRUE,
                                     color_palette = color_palette)
  p.group.spec.2 <- .plot.group.spec(res, 
                                     show.stand = show.stand, mem.crisp = FALSE,
                                     color_palette = color_palette)
  
  result = list()
  result$p.all.spec = p.all.spec
  result$p.cluster.spec = p.cluster.spec
  result$p.group.spec.1 = p.group.spec.1
  result$p.group.spec.2 = p.group.spec.2
  result$res = res
  result$show.stand = show.stand
  
  return(result)
}



#' @name .plot.all.spec
#' @title .plot.all.spec
#' @param res Result of \code{FCM.new()} object.
#' @param show.stand Whether to show the spectra on standarized scale (default as \code{NULL} 
#'   and assigned by \code{res$FD$stand}).
#' @return The result of \code{.plot.all.spec} is a \code{ggplot} object showing all spectra.
#' @noRd
#' 
.plot.all.spec <- function(res, show.stand = NULL){
  
  if(is.null(show.stand)) show.stand = !res$FD$stand
  
  if(show.stand==FALSE){
    Rrs_ <- res$FD$x
    y.lable <- expression(R[rs]~group("[",sr^-1,"]"))
  }else if(show.stand==TRUE){
    Rrs_ <- res$FD$x.stand
    y.lable <- "Normalized Rrs"
  }
  spectra.p <- data.frame(nm=seq(1,nrow(Rrs_)) %>% as.character, Rrs_, 
                         stringsAsFactors = FALSE) %>%
    setNames(., c("name", res$FD$wv %>% as.character)) %>%
    melt(., id='name', variable.name='band', value.name='value') %>%
    level_to_variable()
  spectra.p$band %<>% as.numeric
  
  result <- ggplot() +
    geom_path(data=spectra.p,
              aes(x=band, y=value, group=name, color=name), alpha=0.3) +
    scale_colour_viridis_d() +
    labs(title=paste0("All spectra of dataset N=",as.character(dim(res$FD$x)[1])),
         x="Wavelength [nm]", y=y.lable) +
    theme_bw() +
    theme(legend.position="none",plot.title=element_text(hjust=0.5))
  return(result)
}


#' @name .plot.cluster.spec
#' @title .plot.cluster.spec
#' @param res Result of \code{FCM.new()} object.
#' @param show.stand Whether to show the spectra on standarized scale (default as \code{NULL} 
#'   and assigned by \code{res$FD$stand}).
#' @param show.ribbon Whether to show the ribbon defined by sd values (default as \code{FALSE}).
#' @param color_palette The palette of cluster color. Default as \code{RdYlBu(res$K)}.
#'   In \code{FCMm}, it could be \link{RdYlBu}, \link{Spectral}, \link{HUE} or other color values
#'   with same length of cluster number.
#' @return The result of \code{.plot.cluster.spec} is a \code{ggplot} object showing the spectra of centroids.
#' @noRd
#' 
.plot.cluster.spec <- function(res, 
                               show.stand = NULL, show.ribbon = FALSE, 
                               color_palette = RdYlBu(res$K)){
  
  if(is.null(show.stand)) show.stand = !res$FD$stand
  
  if(length(color_palette) != res$K)
    stop("The length of color_palette shoud be same with cluster number!")
  
  if(show.stand==FALSE){
    Rrs_ <- (as.matrix(res$FD$x))
    y.lable <- expression(R[rs]~group("[",sr^-1,"]"))
  }else if(show.stand==TRUE){
    Rrs_ <- (as.matrix(res$FD$x.stand))
    y.lable <- "Normalized Rrs"
  }
  
  Rrs.cluster.mean <- stats::aggregate(Rrs_, list(res$res.FCM$cluster), mean) %>% melt(., id = 1)
  Rrs.cluster.sd <- stats::aggregate(Rrs_, list(res$res.FCM$cluster), sd) %>% melt(., id = 1) 
  
  Rrs.cluster <- cbind(Rrs.cluster.mean, Rrs.cluster.sd[,3]) %>% 
    setNames(., c("cluster", "band", "mean", "sd")) %>% level_to_variable()
  Rrs.cluster[, 1] %<>% as.character
  Rrs.cluster[, 2] %<>% as.numeric
  
  if(show.ribbon){
    
    result <- ggplot() +
      geom_ribbon(data = Rrs.cluster,
                  aes(x = band, ymin = mean-sd, ymax = mean+sd, 
                      group = cluster, fill = cluster), 
                  color = NA, show.legend = FALSE, alpha = 0.2) + 
      geom_path(data = Rrs.cluster,
                aes(x = band, y = mean, group = cluster, color = cluster), 
                size=1.1, alpha=1.5) +
      scale_color_manual(values = color_palette, breaks = seq(1,res$K)) +
      scale_fill_manual(values = color_palette, breaks = seq(1,res$K)) +
      labs(title = NULL, x = "Wavelength [nm]", y = y.lable, color = "Cluster") +
      theme_bw() +
      theme(legend.position = 'right',
            legend.key.height = unit(0.8,'lines'),
            text = element_text(size=13))
    
  }else{
    
    result <- ggplot() +
      geom_path(data = Rrs.cluster,
                aes(x = band, y = mean, group = cluster, color = cluster),
                size=1.1, alpha=1.5) +
      scale_color_manual(values = color_palette, breaks = seq(1,res$K)) +
      scale_fill_manual(values = color_palette, breaks = seq(1,res$K)) +
      labs(title = NULL, x = "Wavelength [nm]", y = y.lable, color = "Cluster") +
      theme_bw() +
      theme(legend.position = 'right',
            legend.key.height = unit(0.8,'lines'),
            text = element_text(size=13))
    
  }

  return(result)
}


#' @name .plot.group.spec
#' @title .plot.group.spec
#' @param res Result of \code{FCM.new()} object.
#' @param show.stand Whether to show the spectra on standarized scale (default as \code{NULL} 
#'   and assigned by \code{res$FD$stand}).
#' @param mem.crisp Default as \code{TRUE} to decide to present spectra in two colors.
#' @param mem.threshold Default as \code{(K-1)/K} as the threshold to show the different colors.
#' @param color_palette The palette of cluster color. Default as \code{RdYlBu(res$K)}.
#'   In \code{FCMm}, it could be \link{RdYlBu}, \link{Spectral}, \link{HUE} or other color values
#'   with same length of cluster number.
#' @return The result of \code{.plot.group.spec} is a \code{ggplot} object showing the spectra by facet_wrap
#' @noRd
#' 
.plot.group.spec <- function(res, show.stand = NULL,
                             mem.crisp=TRUE, mem.threshold=NULL, 
                             color_palette = RdYlBu(res$K)){
  
  if(is.null(show.stand)) show.stand = !res$FD$stand
  
  if(is.numeric(mem.threshold)){
    if(mem.threshold <= 0 || mem.threshold >= 1)
      stop("The mem.threshold number should be set in (0, 1)")
  }
  if(is.null(mem.threshold)){
    mem.threshold <- (res$K - 1) / res$K
  }
  
  if(length(color_palette) != res$K)
    stop("The length of color_palette shoud be same with cluster number!")
  
  if(show.stand==FALSE){
    Rrs_ <- res$FD$x
    y.lable <- expression(R[rs]~group("[",sr^-1,"]"))
  }else if(show.stand==TRUE){
    Rrs_ <- res$FD$x.stand
    y.lable <- "Normalized Rrs"
  }
  
  Rrs.number <- stats::aggregate(Rrs_[,1], list(res$res.FCM$cluster), length) %>%
    setNames(., c("cluster","number"))
  Rrs.number$cluster_f <- paste("Cluster",Rrs.number$cluster) %>%
    factor(., levels=paste("Cluster", 1:res$K), ordered = TRUE)
  
  Rrs.cluster <- stats::aggregate(Rrs_, list(res$res.FCM$cluster), mean) %>% 
    melt(., id = 1) %>% level_to_variable() %>% 
    setNames(., c("cluster", "band", "value"))
  Rrs.cluster[, 1] %<>% as.character
  Rrs.cluster[, 2] %<>% as.numeric

  Rrs.group <- data.frame(stringsAsFactors = FALSE,
                          cluster = as.character(res$res.FCM$cluster),
                          x = paste0("ID",1:nrow(Rrs_)),
                          Rrs_) %>% setNames(., c("cluster", "ids", colnames(Rrs_))) %>%
    melt(., id=c("cluster","ids")) %>% level_to_variable() %>%
    setNames(., c("cluster", "ids", "band", "value"))
  Rrs.group[,3] %<>% as.numeric
  
  # re-factor
  Rrs.cluster$cluster_f <- paste("Cluster", Rrs.cluster$cluster) %>% 
    factor(., levels=paste("Cluster", 1:res$K), ordered = TRUE)
  Rrs.group$cluster_f <- paste("Cluster", Rrs.group$cluster) %>%
    factor(., levels=paste("Cluster", 1:res$K), ordered = TRUE)
  
  
  # Rrs.ribbon <- data.frame(stringsAsFactors = FALSE,
  #                         cluster = as.character(res$res.FCM$cluster),
  #                         x = paste0("ID",1:nrow(Rrs_)),
  #                         Rrs_) %>% setNames(., c("cluster", "ids", colnames(Rrs_)))
  # Rrs.ribbon_ <- aggregate(Rrs_, list(res$res.FCM$cluster), sd)
  
  if(mem.crisp==TRUE){
    
    # memb threshold
    w_confident = which(apply(res$res.FCM$u, 1, max) >= mem.threshold)
    Rrs.group_ <- data.frame(stringsAsFactors = FALSE,
                             cluster = as.character(res$res.FCM$cluster),
                             x = paste0("ID",1:nrow(Rrs_)),
                             Rrs_)[w_confident,] %>% 
      setNames(., c("cluster", "ids", colnames(Rrs_))) %>%
      melt(., id=c("cluster","ids")) %>% level_to_variable() %>%
      setNames(., c("cluster", "ids", "band", "value"))
    Rrs.group_[,3] %<>% as.numeric
    Rrs.group_$cluster_f <- paste("Cluster", Rrs.group_$cluster) %>%
      factor(., levels=paste("Cluster", 1:res$K), ordered = TRUE)
    
    p <- ggplot() + 
      geom_path(data = Rrs.group, 
                aes(x = band, y = value, group = ids) , color = "gray50", alpha=1, size = 1.3)+
      geom_path(data = Rrs.group_, 
                aes(x = band, y = value, group = ids) , color = "seagreen", alpha=0.3, size = 1.3)+
      geom_path(data = Rrs.cluster, 
                aes(x=band, y=value, group= cluster, color = cluster), alpha = 1.5, size=1.5) + 
      geom_text(data = Rrs.number,
                aes(x=-Inf, y=Inf, label = sprintf("N=%s",number)), vjust=1.8, hjust = -0.2) + 
      scale_color_manual(values = color_palette, breaks = seq(1,res$K)) +
      facet_wrap(~cluster_f) + 
      labs(x = "Wavelength [nm]", y = y.lable, color = "Cluster") + 
      theme_bw() + 
      theme(text = element_text(size=13),
            panel.spacing.x = unit(1.5, "lines"),
            legend.position = "none",
            strip.text = element_text(face="bold"))
    
  }else{
    
    # no memb threshold
    p <- ggplot() + 
      geom_path(data = Rrs.group, 
                aes(x = band, y = value, group = ids),
                color = "gray50", alpha=1, size=1.3)+
      geom_path(data = Rrs.cluster,
                aes(x=band, y=value, group= cluster, color = cluster),
                size = 1.5) + 
      geom_text(data = Rrs.number,
                aes(x=-Inf, y=Inf, label = sprintf("N=%s",number)), vjust=1.8, hjust = -0.2) + 
      scale_color_manual(values = color_palette, breaks = seq(1,res$K)) +
      labs(x = "Wavelength [nm]", y = y.lable, color = "Cluster") + 
      facet_wrap(~cluster_f) + 
      theme_bw() + 
      theme(text = element_text(size=13),
            panel.spacing.x = unit(1.5, "lines"),
            legend.position = "none",
            strip.text = element_text(face="bold"))
    
  }
  
  return(p)
  
}
