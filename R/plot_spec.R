#' @title Plot the result of FCM_m
#' @name plot_spec
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

plot_spec_from_df <- function(df){
  if(names(df) %>% as.numeric %>% anyNA)
    stop('The names of input dataframe cannot be converted to numeric wavelength!')
  tmp.spec <- cbind(nm=seq(1,nrow(df)), df) %>% as.data.frame
  names(tmp.spec)[1] <- "name"
  tmp.spec$name %<>% as.character
  spectra.p <- melt(tmp.spec,id='name',variable.name='band',value.name='value')
  spectra.p$band %<>% levels(.)[.] %>% as.numeric
  result <- ggplot(data=spectra.p,
                   aes(x=band,y=value,group=name,color=name)) +
    geom_path(alpha=0.5) +
    scale_colour_viridis(discrete=TRUE) +
    theme_bw() +
    theme(legend.position='none')
  return(result)
}

#' @export
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
    scale_colour_viridis(discrete=TRUE) +
    labs(title=paste0("All spectra of dataset N=",as.character(dim(res$FD$x)[1])),
         x="Wavelength", y=y.lable) +
    theme_bw() +
    theme(legend.position="none",plot.title=element_text(hjust=0.5))
  return(result)
}

#' @export
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

#' @export
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
