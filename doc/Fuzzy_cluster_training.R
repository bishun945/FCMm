## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
	fig.align = "center",
	message = FALSE,
	warning = FALSE,
	collapse = TRUE,
	comment = "#>"
)

## ----message=F, warning=F-----------------------------------------------------
rm(list=ls())
library(FCMm)
library(ggplot2)
library(magrittr)
library(stringr)
library(dplyr)
data("Nechad2015")
w <- Nechad2015 %>% names %>%
  str_extract(.,pattern="\\d") %>%
  is.na %>% {!.}
wv <- w %>% names(Nechad2015)[.] %>%
  gsub('X','',.) %>% as.numeric
x <- w %>% Nechad2015[,.]
set.seed(1234)
w_sample <- sample_n(seq(nrow(x)) %>% data.frame, 100) %>% as.matrix %>% c
x <- x[w_sample,]
names(x) <- wv
rm(w)

## ----fig.height=4, fig.width=6------------------------------------------------
p.spec <- plot_spec_from_df(x) + 
  labs(x='Wavelength (nm)',y=expression(Rrs~(sr^-1))) + 
  theme_bw() + 
  theme(legend.position='none', text=element_text(size=18))
print(p.spec)

## ----message=FALSE, warning=FALSE, fig.height=4, fig.width=6------------------
library(ppclust)
FD <- FuzzifierDetermination(x, wv, stand=F)
(FD$m.used)
nb <- 4
set.seed(54321) # I just set this seed so that you can re-produce them
result <- FCM.new(FD, nb, fast.mode = T)
summary(result)
result$p.jitter + theme(text = element_text(size=13))

## ----message=FALSE, warning=FALSE, include=TRUE-------------------------------
if(library('cowplot', logical.return = T)){
  library('cowplot')
}else{
  install.packages("cowplot")
  library('cowplot')
}

p.spec <- plot_spec(result, show.stand=T, HABc=NULL)
print(p.spec$p.cluster.spec)
plot_grid(plotlist = p.spec$p.group.spec.2)

