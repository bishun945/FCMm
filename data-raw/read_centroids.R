rm(list=ls())

B2021_norm <- read.csv("./inst/Bi21_centroids_Rrs_norm.csv", 
                        comment.char = "#", check.names = FALSE,
                        stringsAsFactors = FALSE)

save(B2021_norm, file = "./data/B2021_norm.rda", compression_level = 9)



BPHD_norm <- read.csv("./inst/BiPHD_centroids_Rrs_norm.csv",
                      comment.char = "#", check.names = FALSE,
                      stringsAsFactors = FALSE)

save(BPHD_norm, file = "./data/BPHD_norm.rda", compression_level = 9)




Jac17_raw <- read.csv("./inst/Jac17_centroids_Rrs.csv",
                      comment.char = "#", check.names = FALSE,
                      stringsAsFactors = FALSE)

save(Jac17_raw, file = "./data/Jac17_raw.rda", compression_level = 9)




Moo14_raw <- read.csv("./inst/Moo14_centroids_Rrs0-.csv",
                      comment.char = "#", check.names = FALSE,
                      stringsAsFactors = FALSE)

Moo14_raw[, -1] <- 0.52 / (1 / Moo14_raw[, -1] - 1.7)

save(Moo14_raw, file = "./data/Moo14_raw.rda", compression_level = 9)


B2019_raw <- read.csv("./inst/Bi19_centroids_Rrs_raw.csv",
                       comment.char = "#", check.names = FALSE,
                       stringsAsFactors = FALSE)
B2019_norm[, -1] <- B2019_raw[,-1] / trapz2(B2019_raw[,-1])

save(B2019_norm, file = "./data/B2019_norm.rda", compression_level = 9)


H2017_lgnorm <- read.csv("./inst/ONNS_lg_centroids.csv",
                         comment.char = "#", check.names = FALSE,
                         stringsAsFactors = FALSE)

save(H2017_lgnorm, file = "./data/H2017_lgnorm.rda", compression_level = 9)



# ----------------     Import Moore hyperspectral mean and covariance data

owtMean <- stars::read_stars("./inst/Moore2014/rrs_lake_owt_means_131213.hdf")[[1]]
owtWavelen <- as.numeric(readr::read_lines("./inst/Moore2014/wl_hyper.txt"))
colnames(owtMean) <- owtWavelen
owtCov  <- stars::read_stars("./inst/Moore2014/rrs_lake_owt_cov_131213.hdf")[[1]]
note <- "The mean value is Rrs with unit sr-1, wavelength is with nm"

owtMoore2014 <- list(owtMean = owtMean, 
                     owtCov = owtCov, 
                     owtWavelen = owtWavelen,
                     note = note)
save(owtMoore2014, file = "./data/owtMoore2014.rda", compression_level = 9)

# -----------------    Import Martin's data

nc <- ncdf4::nc_open("./inst/MH2017/class_means.nc")
owtWavelen <- nc$dim$wavelength$vals
owtWavelen <- owtWavelen[-length(owtWavelen)]
owtMean <- ncdf4::ncvar_get(nc, "cluster_means")
owtArea <- owtMean[, ncol(owtMean)]
owtMean <- owtMean[, -ncol(owtMean)]
colnames(owtMean) <- owtWavelen
ncdf4::nc_close(nc)

nc <- ncdf4::nc_open("./inst/MH2017/class_cov.nc")
owtCov  <- ncdf4::ncvar_get(nc, "covariance")
owtCov <- owtCov[-12, -12, ]
ncdf4::nc_close(nc)

note <- "The mean value is sum-normalized based on log10(Rrs + 1), wavelength is with nm. The normalization is sum of the 11 bands."

owtMH2017 <- list(owtMean = owtMean,
                  owtCov  = owtCov,
                  owtWavelen = owtWavelen,
                  note = note)
save(owtMH2017, file = "./data/owtMH2017.rda", compression_level= 9)

# ------------------   Import Jackson's data (for Mahalanobis distance calculation)

Jac17_raw <- read.csv("./inst/Jac17_centroids_Rrs.csv",
                      comment.char = "#", check.names = FALSE,
                      stringsAsFactors = FALSE) %>%
  .[,-1] %>% as.matrix()

Jac17_sd <- read.csv("./inst/Jac17_sd_Rrs.csv",
                      comment.char = "#", check.names = FALSE,
                      stringsAsFactors = FALSE) %>%
  .[,-1] %>% as.matrix() %>% {.^2}

owtWavelen <- as.numeric(colnames(Jac17_raw))

owtJackson2017 <- list(owtMean = Jac17_raw,
                       owtCov  = Jac17_sd,
                       owtWavelen = owtWavelen,
                       note = "Only diagonals of variance. Not covariance matrix.")

save(owtJackson2017, file = "./data/owtJack2017.rda", compression_level = 9)




# ------------------- Import Spyrakos's normalized centroids

Spy18_norm <- read.csv("./inst/Spyrakos_OWT_inland.csv", check.names = FALSE)[,-1] %>% as.matrix()
plot_spec_from_df(Spy18_norm)
owtSyprakos2018 <- list(owtmean = Spy18_norm,
                        owtWavelen = as.numeric(colnames(Spy18_norm)),
                        note = "Normalization are based on trapz area")
save(owtSyprakos2018, file = "./data/owtSpyrakos2018.rda", compression_level = 9)
