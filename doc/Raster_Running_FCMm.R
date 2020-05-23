## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
	fig.align = "center",
	message = FALSE,
	warning = FALSE,
	collapse = TRUE,
	comment = "#>"
)

## ----message=FALSE, warning=FALSE---------------------------------------------
rm(list=ls())
library(FCMm)
library(ggplot2)
library(magrittr)
library(dplyr)
data("OLCI_TH")
data("Bi_clusters")

class(OLCI_TH)

# convert RasterBrick to imdf with NA value removed and xy coordinates recorded
imdf <- raster::as.data.frame(OLCI_TH, na.rm=T, xy=TRUE)

x_name <- which(names(imdf) == "x")
y_name <- which(names(imdf) == "y")
imdf <- cbind(imdf[,c(x_name, y_name)], imdf[,-c(x_name, y_name)])

names(imdf)[-c(1,2)] <- c(412.5,442.5,490,510,
                           560,620,665,673.75,
                           708.75,753.75,865,885)

## ----fig.height=4, fig.width=6------------------------------------------------
set.seed(54321)
sample_n(imdf[,-c(1,2)],50) %>% 
  plot_spec_from_df(.) + 
  labs(x='Wavelength (nm)',y=expression(Rrs~(sr^-1))) + 
  theme_bw() + 
  theme(legend.position='none', text=element_text(size=18))

## ----message=TRUE, warning=FALSE, fig.height=4, fig.width=6-------------------
res <- generate_param(c(413,443,490,510,560,620,665,674,709,754,865,885))

im_result <- apply_to_image(input=OLCI_TH, res=res, title.name="Test_image",
                            Chla_est=T, output_image=F)
summary(im_result)

