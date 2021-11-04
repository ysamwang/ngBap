runInd <- 1
args <- commandArgs(TRUE)
for(i in 1:length(args)){
  eval(parse(text = args[[i]]))
}
print(runInd)

library("ngBap")
library("pcalg")
library("greedyBAPs")
library("BCD")
library("dHSIC")

#######
# 7/5/21: Include dhsic test
#

one.run.baps <- function(p, n,  dist, K = 4, d, b, restrict = 1, signs = T){
  dat <- ngBap::rBAP(n = n, p = p, dist = dist, d = d, b = b, ancestral = F, shuffle = T, signs = signs)

  Y <- dat$Y

  out05EL <- ngBap::bang(Y, K = K, level = .05, verbose = F, restrict = restrict, testType = "empLik")
  out01EL <- ngBap::bang(Y, K = K, level = .01, verbose = F, restrict = restrict, testType = "empLik")
  out001EL <- ngBap::bang(Y, K = K, level = .001, verbose = F, restrict = restrict, testType = "empLik")

  bangMod05EL <- BCD::ricf(out05EL$dEdge, out05EL$bEdge, t(dat$Y), tol = 1e-8)
  bangMod01EL <- BCD::ricf(out01EL$dEdge, out01EL$bEdge, t(dat$Y), tol = 1e-8)
  bangMod001EL <- BCD::ricf(out001EL$dEdge, out001EL$bEdge, t(dat$Y), tol = 1e-8)



  gbap_out <- greedyBAPs::greedySearch(cov(Y), n = n, n.restarts = 100)

  gbapB <- t(ifelse(gbap_out$final.bap == 1, 1, 0))
  gbapOm <- ifelse(gbap_out$final.bap == 100, 1, 0) + diag(rep(1, p))


  trueMod <- BCD::ricf(dat$directEdges,   dat$bidirectEdges, t(dat$Y), tol = 1e-8)
  gbapMod <- BCD::ricf(gbapB, gbapOm, t(dat$Y), tol = 1e-8)




  if(n < 5000){

    out05DH <- ngBap::bang(Y, K = K, level = .05, verbose = F, restrict = restrict, testType = "dhsic")
    out01DH <- ngBap::bang(Y, K = K, level = .01, verbose = F, restrict = restrict, testType = "dhsic")
    out001DH <- ngBap::bang(Y, K = K, level = .001, verbose = F, restrict = restrict, testType = "dhsic")


    bangMod01DH <- BCD::ricf(out01DH$dEdge, out01DH$bEdge, t(dat$Y), tol = 1e-8)
    bangMod001DH <- BCD::ricf(out001DH$dEdge, out001DH$bEdge, t(dat$Y), tol = 1e-8)
    bangMod05DH <- BCD::ricf(out05DH$dEdge, out05DH$bEdge, t(dat$Y), tol = 1e-8)

    bangScore01DH = mean(mvtnorm::dmvnorm(dat$Y, sigma = bangMod01DH$SigmaHat, log = T)) - (sum(out01DH$dEdge) + (sum(out01DH$bEdge) - p)/2) * log(n) / n
    bangScore001DH = mean(mvtnorm::dmvnorm(dat$Y, sigma = bangMod001DH$SigmaHat, log = T)) - (sum(out001DH$dEdge) + (sum(out001DH$bEdge) - p)/2) * log(n) / n
    bangScore05DH = mean(mvtnorm::dmvnorm(dat$Y, sigma = bangMod05DH$SigmaHat, log = T)) - (sum(out05DH$dEdge) + (sum(out05DH$bEdge) - p)/2) * log(n) / n

    bangExact01DH = sum(abs(out01DH$bEdge - dat$bidirectEdges) + abs(out01DH$dEdge - dat$directEdges)) == 0
    bangExact001DH = sum(abs(out001DH$bEdge - dat$bidirectEdges) + abs(out001DH$dEdge - dat$directEdges)) == 0
    bangExact05DH = sum(abs(out05DH$bEdge - dat$bidirectEdges) + abs(out05DH$dEdge - dat$directEdges)) == 0

  } else {

    bangScore01DH = NA
    bangScore001DH = NA
    bangScore05DH = NA

    bangExact01DH = NA
    bangExact001DH = NA
    bangExact05DH = NA
  }






  ret <- list(trueScore = mean(mvtnorm::dmvnorm(dat$Y, sigma = trueMod$SigmaHat, log = T)) - (sum(dat$directEdges) + (sum(dat$bidirectEdges) - p)/2) * log(n) / n,
              gbsScore = mean(mvtnorm::dmvnorm(dat$Y, sigma = gbapMod$SigmaHat, log = T)) - (sum(gbapB) + (sum(gbapOm) - p)/2) * log(n) / n,

              bangScore05EL = mean(mvtnorm::dmvnorm(dat$Y, sigma = bangMod05EL$SigmaHat, log = T)) - (sum(out05EL$dEdge) + (sum(out05EL$bEdge) - p)/2) * log(n) / n,
              bangScore01EL = mean(mvtnorm::dmvnorm(dat$Y, sigma = bangMod01EL$SigmaHat, log = T)) - (sum(out01EL$dEdge) + (sum(out01EL$bEdge) - p)/2) * log(n) / n,
              bangScore001EL = mean(mvtnorm::dmvnorm(dat$Y, sigma = bangMod001EL$SigmaHat, log = T)) - (sum(out001EL$dEdge) + (sum(out001EL$bEdge) - p)/2) * log(n) / n,

              bangExact05EL = sum(abs(out05EL$bEdge - dat$bidirectEdges) + abs(out05EL$dEdge - dat$directEdges)) == 0,
              bangExact01EL = sum(abs(out01EL$bEdge - dat$bidirectEdges) + abs(out01EL$dEdge - dat$directEdges)) == 0,
              bangExact001EL = sum(abs(out001EL$bEdge - dat$bidirectEdges) + abs(out001EL$dEdge - dat$directEdges)) == 0,


              bangScore05DH,
              bangScore01DH,
              bangScore001DH,


              bangExact05DH,
              bangExact01DH,
              bangExact001DH
              )

  return(ret)

}



##########################################################
library(parallel)
sim.size <- 250
clust.size <- 7
p <- 6


### param.grid size: 288 ###
param.grid <- expand.grid(n = c(500, 1000, 1500, 10000, 25000, 50000), dist = c("gamma", "unif", "t", "lognormal"),
                          regime = c("sparse", "medium", "dense"), signs = c(T, F))
param.grid <- rbind(param.grid, param.grid)




n <- param.grid[runInd, 1]
dist <- as.character(param.grid[runInd, 2])
regime <- as.character(param.grid[runInd, 3])
signs <- param.grid[runInd, 4]

cover <- data.frame(p = rep(p, sim.size), n = rep(n, each = sim.size), dist = rep(dist, each = sim.size),
                    regime = rep(regime, sim.size), signs = rep(signs, sim.size))


if (dist == "gauss"){
  K <- 4
} else if (dist == "gamma") {
  K <- 3
} else if (dist == "unif") {
  K <- 4
} else if (dist == "t") {
  K <- 6
} else if (dist == "lognormal") {
  K <- 3
} else if (dist == "laplace") {
  K <- 4
}



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

cl <- makeCluster(clust.size)
clusterExport(cl, ls())


out <- parSapply(cl, 1:sim.size, function(i){one.run.baps(p = p, n = n, dist = dist, K = K, d = d, b = b, restrict = 1, signs = signs)})

out <- t(out)


one.run.baps(p = p, n = n, dist = dist, K = K, d = d, b = b, restrict = 1, signs = signs)

end <- proc.time()
stopCluster(cl)

cover <- cbind(cover, data.matrix(out))
cover <- data.frame(apply(cover, MAR = 2, unlist))

write.table(cover, paste("../resultsNew/bap_",runInd,".csv", sep = ""), row.names = F, sep = ",")
