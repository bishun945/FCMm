## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>",
                      fig.align='center')

## ----message=FALSE, warning=FALSE---------------------------------------------
rm(list=ls())
library(FCMm)
library(tidyverse)
data("WaterSpec35")
data("Bi_clusters")
Rrs <- WaterSpec35[,3:17]

## ----fig.height=4, fig.width=6------------------------------------------------
qplot(data=WaterSpec35, x=Chla, geom='density', xlab='The raw Chla')
qplot(data=WaterSpec35, x=Chla, geom='density', log='x', xlab='Log10 transformed Chla')
p.spec <- plot_spec_from_df(Rrs) + 
  labs(x='Wavelength (nm)',y=expression(Rrs~(sr^-1))) + 
  theme_bw() + 
  theme(legend.position='none', text=element_text(size=18))
print(p.spec)

## ----message=FALSE, warning=FALSE, fig.height=4, fig.width=6------------------
result <- apply_FCM_m(Rrs=Rrs, option.plot=TRUE)
summary(result)
result$p.group

## ----fig.height=4, fig.width=6------------------------------------------------
dt_Chla <- FCM_m_Chla_estimation(Rrs=data.frame(Rrs665=Rrs$`665`,
                                                Rrs709=Rrs$`708.75`,
                                                Rrs754=Rrs$`753.75`),
                                 U=result$u)
dt_Chla$cluster <- result$cluster %>% as.character
dt_Chla$Chla_true <- WaterSpec35$Chla

options(scipen=10000)

subset(dt_Chla, select=c('cluster','Chla_true','BR','TBA','C6','conc.Blend')) %>%
  reshape2::melt(., id=c('cluster','Chla_true')) %>%
  ggplot(data=.) + 
  geom_point(aes(x=Chla_true,y=value,group=cluster,color=cluster),
             alpha=0.8, size=4) +
  scale_x_log10(limits=c(1,800)) + 
  scale_y_log10(limits=c(1,800)) +
  scale_color_manual(values=heatmaply::RdYlBu(result$K)) + 
  labs(x='True value of Chla concentration (ug/L)',
       y='Estimated value of Chla concentration (ug/L)',
       color='Cluster') + 
  geom_abline(intercept=0, slope=1, linetype=2) + 
  facet_wrap(~variable, nrow=2) + 
  theme_bw() + 
  theme(axis.text.x.bottom = element_text(hjust=1))

MAPEs <- NULL
i <- 1
for(model in c('BR','TBA','C6','conc.Blend')){
  MAPEs[i] <- cal.metrics(x=dt_Chla$Chla_true %>% log10,
                          y=dt_Chla[,model]  %>% log10,
                          name='MAPE')
  names(MAPEs)[i] <- model
  i <- i + 1
}
print(MAPEs)

MAEs <- NULL
i <- 1
for(model in c('BR','TBA','C6','conc.Blend')){
  MAEs[i] <- cal.metrics(x=dt_Chla$Chla_true,
                         y=dt_Chla[,model],
                         name='MAE2')
  names(MAEs)[i] <- model
  i <- i + 1
}
print(MAEs)

## ----fig.height=4, fig.width=6------------------------------------------------
Rrs_sub <- subset(Rrs, select=c(`412.5`,`442.5`,`490`,`510`,
                                `560`,`620`,`665`,`673.75`,
                                `708.75`,`753.75`,`865`,`885`))
wavelength.sub <- c(412.5,442.5,490,510,
                    560,620,665,673.75,
                    708.75,753.75,865,885)
Rrs_clusters.sub <- which(names(Rrs_clusters.default) != 'X400' &
                            names(Rrs_clusters.default) != 'X681' & 
                            names(Rrs_clusters.default) != 'X779') %>%
  Rrs_clusters.default[,.]

# Note the parameter settings in this function `default.cluster=F`
result_sub <- apply_FCM_m(Rrs=Rrs_sub, wavelength=wavelength.sub,
                          Rrs_clusters=Rrs_clusters.sub,
                          stand=F, default.cluster=F, option.plot=T)
result_sub$p.group

