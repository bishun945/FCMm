# Why sort-based scoring method could not be so NB?

rm(list=ls())

library(magrittr)
library(ggplot2)
library(FCMm)

#################################################### Generate toy data

x_true <- 50 + 10^seq(0, 3, 0.01)
# normal error added
set.seed(1234)
dt_normal_error <- data.frame(
  x  = x_true,
  y  = x_true * (1+runif(length(x_true), min=-0.05, max=0.05)),
  y1 = x_true * (1+runif(length(x_true), min=0.2, max=0.3)),
  y2 = x_true * (1-runif(length(x_true), min=0.2, max=0.3)),
  y3 = x_true * (1+runif(length(x_true), min=0.4, max=0.4)),
  y4 = x_true * (1-runif(length(x_true), min=0.4, max=0.4)),
  y5 = x_true * (1+runif(length(x_true), min=0.6, max=0.6)),
  y6 = x_true * (1-runif(length(x_true), min=0.6, max=0.6))
)


################################################## Begin assessment

dt_test <- dt_normal_error
summary(dt_test)

dt_p <- reshape2::melt(dt_test, id="x")

ggplot(dt_p) + 
  geom_point(aes(x=x, y=(x-value), color=variable)) + 
  geom_hline(yintercept=0) + 
  scale_x_log10()

ggplot(dt_p) + 
  geom_point(aes(x=x, y=value, color=variable)) + 
  geom_abline(slope=1, intercept=0) +
  scale_x_log10() + scale_y_log10()

memb = rep(1, length(x_true)) %>% as.matrix

Asses <- 
  Assessment_via_cluster(pred=dt_test[,-1],
                         meas=dt_test$x,
                         memb=memb,
                         metrics = c('MAE', 'BIAS',
                                     'CMAPE','CMRPE'),
                         total=FALSE,
                         na.process=TRUE,
                         hard.mode=TRUE,
                         cal.precision = TRUE,
                         log10=TRUE,
                         plot.col=TRUE)

er <- data.frame(
  A.MAE  = Asses$MAE     %>% round(3) %>% t %>% as.numeric, 
  A.BIAS = Asses$BIAS    %>% round(3) %>% t %>% as.numeric, 
  A.CAPE = Asses$CMAPE   %>% round(2) %>% t %>% as.numeric, 
  A.CRPE = Asses$CMRPE   %>% round(2) %>% t %>% as.numeric,
  P.MAE  = Asses$MAE_p   %>% round(3) %>% t %>% as.numeric, 
  P.BIAS = Asses$BIAS_p  %>% round(3) %>% t %>% as.numeric, 
  P.CAPE = Asses$CMAPE_p %>% round(2) %>% t %>% as.numeric, 
  P.CRPE = Asses$CMRPE_p %>% round(2) %>% t %>% as.numeric
)

rownames(er) <- names(dt_test)[-1]

er

score_sort <- apply(er, 2, Score_algorithms_sort) %>% rowSums()

score_inte <- apply(er, 2, function(x) Score_algorithms_interval(x, reward.punishment=FALSE)$score) %>% rowSums()


score_sort
score_inte


plot(score_sort, score_inte)


######################################## Simulate Er to Score

rm(list=ls())

er <- data.frame(
  A = c(1,  10, 50, 100),
  B = c(10,  1, 100, 10),
  C = c(1,  10, 100, 50),
  D = c(100, 1, 10,  50),
  E = c(10,  1, 50, 100)
)

(sort_sc <- apply(er, 2, Score_algorithms_sort) %>% rowSums())

(scale_sc <- apply(er, 2, function(x) scales::rescale(x)) %>% rowSums() %>% round(2))





