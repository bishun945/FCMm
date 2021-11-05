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






