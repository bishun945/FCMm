## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>",
                      fig.align='center')

## ----message=F, warning=F-----------------------------------------------------
rm(list=ls())
library(FCMm)
library(tidyverse)
data("Nechad2015")
w <- Nechad2015 %>% names %>%
  str_extract(.,pattern="\\d") %>%
  is.na %>% {!.}
wv <- w %>% names(Nechad2015)[.] %>%
  gsub('X','',.) %>% as.numeric
x <- w %>% Nechad2015[,.]
names(x) <- wv
rm(w)

## ----fig.height=4, fig.width=6------------------------------------------------
p.spec <- plot_spec_from_df(x) + 
  labs(x='Wavelength (nm)',y=expression(Rrs~(sr^-1))) + 
  theme_bw() + 
  theme(legend.position='none', text=element_text(size=18))
print(p.spec)

## ----message=FALSE, warning=FALSE, fig.height=4, fig.width=6------------------
FD <- FuzzifierDetermination(x, wv, stand=F)
(FD$m.used)
nb <- 4
set.seed(54321) # I just set this seed so that you can re-produce them
result <- FCM.new(FD, nb, fast.mode = T)
summary(result)
print(result$p.jitter)

## ----message=FALSE, warning=FALSE, include=TRUE-------------------------------
p.spec <- plot_spec(result, show.stand=F, HABc=NULL)
# print(p.spec$p.all.spec)
# print(p.spec$p.cluster.spec)
# print(p.spec$p.group.spec.2)

## ----message=FALSE, warning=FALSE, fig.height=5, fig.width=6------------------
library(magrittr)
# melt data
tmp <- data.frame(cluster = result$res.FCM$cluster %>% as.character,
                  X = FD$x.stand)
names(tmp)[2:ncol(tmp)] <- as.character(wv)
tmp.p <- reshape2::melt(tmp, id='cluster', variable.name='band', value.name='value')
tmp.p$cluster %<>% levels(.)[.]
tmp.p$band %<>% levels(.)[.] %>% as.numeric

tmp2 <- data.frame(name = seq(1,nrow(tmp)) %>% as.character,
                   tmp[,2:ncol(tmp)])
names(tmp2)[2:ncol(tmp2)] <- as.character(wv)
tmp2.p <- reshape2::melt(tmp2, id='name', variable.name='band', value.name='value')
tmp2.p$name %<>% levels(.)[.]
tmp2.p$band %<>% levels(.)[.] %>% as.numeric
tmp.p$name <- tmp2.p$name  
tmp.p$cluster <- paste0('Cluster ',tmp.p$cluster)

# plot
cp2 <- cp <- heatmaply::RdYlBu(result$K)
names(cp) <- seq(1,result$K)
names(cp2) <- paste0('Cluster ',seq(1,result$K))
p.group.facet <- ggplot(data=tmp.p) + 
  geom_path(aes(x=band, y=value, group=name, color=cluster), alpha=0.9) + 
  scale_color_manual(values=cp2,
                     labels=paste0('Cluster ',seq(result$K))) + 
  labs(color=paste0('Cluster ',seq(result$K)),y='Normalized_scale') + 
  facet_wrap(~cluster, nrow=1) + 
  theme_bw() + 
  theme(legend.position='none')
# plot(p.group.facet)

Nechad2015$cluster <- result$res.FCM$cluster %>% as.character

tmp <- data.frame(cluster=Nechad2015$cluster,Chla=Nechad2015$`X.Chl_a..ug.L.`) %>% na.omit
tmp$cluster <- reorder(tmp$cluster, tmp$Chla, mean)

p.boxplot <- ggplot(data=tmp, aes(x=cluster,y=Chla,fill=cluster)) + 
  geom_boxplot(outlier.color=NA) + scale_y_log10() + 
  scale_fill_manual(values=cp) + 
  theme_bw() + theme(legend.position='none')

p.cluster.spec.n <- plot_spec_from_df(result$res.FCM$v) + 
  geom_path(size=1.5) + labs(y="Normalized_scale") + 
  scale_color_manual(values=cp) +
  theme(legend.position='right')

library(gridExtra)
pm <- arrangeGrob(grobs=list(p.cluster.spec.n, p.boxplot),nrow=1) %>%
  arrangeGrob(grobs=list(., p.group.facet),
                          nrow=2,top='Nechad2015: Normalized cluster result')
plot(pm)


