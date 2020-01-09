## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>",
                      fig.align='center')

## ----message=FALSE, warning=FALSE---------------------------------------------
rm(list=ls())
library(FCMm)
library(tidyverse)
data("OLCI_TH")
data("Bi_clusters")

class(OLCI_TH)

# convert RasterBrick to imdf with NA value removed and xy coordinates recorded
imdf <- raster::as.data.frame(OLCI_TH, na.rm=T, xy=TRUE)
names(imdf)[c(-1,-2)] <- c(412.5,442.5,490,510,
                           560,620,665,673.75,
                           708.75,753.75,865,885)

## ----fig.height=4, fig.width=6------------------------------------------------
set.seed(54321)
sample_n(imdf[,c(-1,-2)],50) %>% 
  plot_spec_from_df(.) + 
  labs(x='Wavelength (nm)',y=expression(Rrs~(sr^-1))) + 
  theme_bw() + 
  theme(legend.position='none', text=element_text(size=18))

## ----message=FALSE, warning=FALSE, fig.height=4, fig.width=6------------------
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

res <- generate_param(c(413,443,490,510,560,620,665,674,709,754,865,885))


im_result <- apply_to_image(input=OLCI_TH, res=res, title.name="Test_image",
                            Chla_est=T, output_image=F)
summary(im_result)

