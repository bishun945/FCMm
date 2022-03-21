apply_FCM_m_maha <- function(Rrs, global_cov = FALSE) {
  
  v <- BPHD_raw
  k <- nrow(v) # cluster number
  
  BandName <- colnames(v)
  
  v <- v[, BandName]
  x <- Rrs[, BandName]
  DF <- ncol(x)
  
  DistMaha <- MembMaha <- matrix(NA, ncol = k, nrow = nrow(x))
  
  if(global_cov) {
    
    s <- BPHD_cov$cov_glob
    s <- s[BandName, BandName]
    
    for(i in 1:k) {
      
      d <- mahalanobis(x, 
                       center = as.numeric(v[i,]),
                       cov = limSolve::Solve(s),
                       inverted = TRUE)
      DistMaha[,i] <- d
      MembMaha[,i] <- 1 - pchisq(d, df = DF)
    }
      
  } else {
    
    s <- BPHD_cov$cov_clus
    s <- s[BandName, BandName, ]
    
    for(i in 1:k) {
      d <- mahalanobis(x, 
                       center = as.numeric(v[i,]),
                       cov = limSolve::Solve(s[,,i]),
                       inverted = TRUE)
      DistMaha[,i] <- d
      MembMaha[,i] <- 1 - pchisq(d, df = DF)
    }
    
  }
  
  MembMaha[MembMaha < 1e-4] <- 0
  
  ClusMaha <- apply(MembMaha, 1, which.max)
  
  MembTotal <- apply(MembMaha, 1, sum)  
  MembRenor <- MembMaha / MembTotal
  
  result <- list(
    DistMaha = DistMaha,
    MembMaha = MembMaha,
    ClusMaha = ClusMaha,
    MembTotal = MembTotal,
    MembRenor = MembRenor
  )
  
  return(result)
  
}