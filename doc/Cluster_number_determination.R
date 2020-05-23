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
x <- sample_n(x, 100)
names(x) <- wv
rm(w)

## ----fig.height=5, fig.width=6------------------------------------------------
p.spec <- plot_spec_from_df(x) + 
  labs(x='Wavelength (nm)',y=expression(Rrs~(sr^-1))) + 
  theme_bw() + 
  theme(legend.position='none', text=element_text(size=13))
print(p.spec)

## -----------------------------------------------------------------------------
FD <- FuzzifierDetermination(x, wv, stand=F)
summary(FD)

## ----fig.height=3, fig.width=6, message=FALSE, warning=FALSE------------------
library(ppclust)
library(fclust)
set.seed(1234)
FD <- x %>% 
  # sample_n(., size=300) %>%
  FuzzifierDetermination(., wv, stand=F)
nb_min <- 3
nb_max <- 10
idxsf <- idxpe <- idxpc <- idxmpc <- seq(nb_min,nb_max,1)
i <- 1
for(nb in nb_min:nb_max){
  res <- FCM.new(FD, nb, fast.mode=T) # open the fast mode
  # print(res$p.jitter)
  tmp <- ppclust2(res$res.FCM, otype="fclust")
  idxsf[i] <- SIL.F(tmp$Xca, tmp$U, alpha=1) # optimal with maximum value
  idxpe[i] <- PE(tmp$U) # optimal with minimum value
  idxpc[i] <- PC(tmp$U) # optimal with maximum value
  idxmpc[i] <- MPC(tmp$U) # optimal with maximum value
  i <- i + 1
}

dt <- data.frame(nb=seq(nb_min,nb_max,1),idxsf,idxpe,idxpc,idxmpc)
opt.num <- c(apply(dt,2,which.max)[-c(1,3)]+1,apply(dt,2,which.min)[3]+1)

dt %>% reshape2::melt(., id='nb') %>% 
  ggplot(data=.,aes(x=nb,y=value,group=variable,color=variable)) + 
  geom_path() + 
  labs(x = "Cluster number", y = 'Metric values') +
  facet_wrap(~variable, scales='free_y', nrow=2) + 
  theme_bw() + 
  theme(text=element_text(size=13), legend.position='none')

