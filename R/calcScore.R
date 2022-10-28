### Functions used to help run simulations with RCD
## Called from the python driver
calcScore <- function(Y, trueDirect, trueBidirect, estDirect, estBidirect){
  p <- ncol(Y)
  n <- nrow(Y)

  trueSig = BCD::ricf(trueDirect, trueBidirect, t(Y), tol = 1e-8)
  estSig = BCD::ricf(estDirect, estBidirect, t(Y), tol = 1e-8)

  trueScore = mean(mvtnorm::dmvnorm(Y, sigma = trueSig$SigmaHat, log = T)) - (sum(trueDirect) + (sum(trueBidirect) - p)/2) * log(n) / n
  estScore = mean(mvtnorm::dmvnorm(Y, sigma = estSig$SigmaHat, log = T)) - (sum(estDirect) + (sum(estBidirect) - p)/2) * log(n) / n
  return(c(trueScore, estScore))
}


### Functions used to help run simulations with RCD
## Called from the python driver
rcdHelper_rBAP <- function(runInd){
  p <- 6
  ### param.grid size: 288 ###
  param.grid <- expand.grid(n = c(500, 1000, 1500), dist = c("gamma", "unif", "t", "lognormal"),
                            regime = c("sparse", "medium", "dense"), signs = c(T, F))
  param.grid <- rbind(param.grid, param.grid)




  n <- param.grid[runInd, 1]
  dist <- as.character(param.grid[runInd, 2])
  regime <- as.character(param.grid[runInd, 3])
  signs <- param.grid[runInd, 4]



  if(regime == "sparse"){

    d <- p/2
    b <- p/2

  } else if (regime == "medium") {

    d <- p
    b <- p

  } else {

    d <- p * 1.5
    b <- p

  }
  dat <- ngBap::rBAP(n = n, p = p, dist = dist, d = d, b = b, ancestral = F, shuffle = T, signs = signs)

  return(dat)

}


