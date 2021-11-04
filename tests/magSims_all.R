runInd <- 1
args <- commandArgs(TRUE)
for(i in 1:length(args)){
  eval(parse(text = args[[i]]))
}
print(runInd)

library("ngBap")
library("pcalg")
library("greedyBAPs")


one.run.mags <- function(p, n,  dist, K = 4, d, b, restrict = 1, signs = T){

  dat <- ngBap::rBAP(n = n, p = p, dist = dist, d = d, b = b, ancestral = T, shuffle = T, signs = signs)

  Y <- dat$Y


  truePAG <- ngBap::mag2pag(dat$direc, dat$bidirectEdges)

  out05EL <- ngBap::bang(Y, K = K, level = .05, verbose = F, restrict = restrict, testType = "empLik")
  out01EL <- ngBap::bang(Y, K = K, level = .01, verbose = F, restrict = restrict, testType = "empLik")


  bangPAG05EL <-ngBap::mag2pag(out05EL$dEdge, out05EL$bEdge)
  bangPAG01EL <-ngBap::mag2pag(out01EL$dEdge, out01EL$bEdge)




  gbap_out <- greedyBAPs::greedySearch(cov(Y), n = n, n.restarts = 100)

  gbapB <- t(ifelse(gbap_out$final.bap == 1, 1, 0))
  gbapOm <- ifelse(gbap_out$final.bap == 100, 1, 0) + diag(rep(1, p))
  gbapPAG <- ngBap::mag2pag(gbapB, gbapOm)

  fciPAG01 <- pcalg::fciPlus(suffStat = list(C = cor(Y), n = n), indepTest = pcalg::gaussCItest, alpha = .01, p = p)
  fciPAG05 <- pcalg::fciPlus(suffStat = list(C = cor(Y), n = n), indepTest = pcalg::gaussCItest, alpha = .05, p = p)




  if(n < 5000){

    out05DH <- ngBap::bang(Y, K = K, level = .05, verbose = F, restrict = restrict, testType = "dhsic")
    out01DH <- ngBap::bang(Y, K = K, level = .01, verbose = F, restrict = restrict, testType = "dhsic")


    bangPAG05DH <-ngBap::mag2pag(out05DH$dEdge, out05DH$bEdge)
    bangPAG01DH <-ngBap::mag2pag(out01DH$dEdge, out01DH$bEdge)

    bangMAG05DH = all(bangPAG05DH@amat == truePAG@amat)
    bangMAG01DH = all(bangPAG01DH@amat == truePAG@amat)

    bangExact01DH = sum(abs(out01DH$bEdge - dat$bidirectEdges) + abs(out01DH$dEdge - dat$directEdges)) == 0
    bangExact05DH = sum(abs(out05DH$bEdge - dat$bidirectEdges) + abs(out05DH$dEdge - dat$directEdges)) == 0



  } else {

    bangScore01DH = NA
    bangScore05DH = NA

    bangExact01DH = NA
    bangExact05DH = NA

    bangMAG05DH = NA
    bangMAG01DH = NA
  }


  ret <- list(bangMAG05EL = all(bangPAG05EL@amat == truePAG@amat),
              bangMAG01EL = all(bangPAG01EL@amat == truePAG@amat),

              bangMAG05DH = bangMAG05DH,
              bangMAG01DH = bangMAG01DH,

              gbsMAG =  all(gbapPAG@amat == truePAG@amat),

              fciMAG05 = all(fciPAG05@amat == truePAG@amat),
              fciMAG01 = all(fciPAG01@amat == truePAG@amat),


              bangExact05EL = sum(abs(out05EL$bEdge - dat$bidirectEdges) + abs(out05EL$dEdge - dat$directEdges)) == 0,
              bangExact01EL = sum(abs(out01EL$bEdge - dat$bidirectEdges) + abs(out01EL$dEdge - dat$directEdges)) == 0,

              bangExact05DH = bangExact05DH,
              bangExact01DH = bangExact01DH
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

out <- parSapply(cl, 1:sim.size, function(i){one.run.mags(p = p, n = n, dist = dist, K = K, d = d, b = b, restrict = 1, signs = signs)})

out <- t(out)



stopCluster(cl)



cover <- cbind(cover, data.matrix(out))
cover <- data.frame(apply(cover, MAR = 2, unlist))

write.table(cover, paste("../resultsNew/anc_",runInd,".csv", sep = ""), row.names = F, sep = ",")
