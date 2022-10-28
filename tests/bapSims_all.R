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
# 10/21/22: Update to include RCD call
#

calc_shd <- function(B, D, trueB, trueD){
  p <- ncol(B)
  correctEdges <- (B == trueB) & (D == trueD)

  return(p * (p-1)/2 - correctEdges)
}

one.run.baps <- function(p, n,  dist, K = 4, d, b, restrict = 1, signs = T, simInd, runInd){
  dat <- ngBap::rBAP(n = n, p = p, dist = dist, d = d, b = b, ancestral = F, shuffle = T, signs = signs)

  Y <- dat$Y

  timeEL05 <- system.time(out05EL <- ngBap::bang(Y, K = K, level = .05, verbose = F, restrict = restrict, testType = "empLik"))[3]
  timeEL01 <- system.time(out01EL <- ngBap::bang(Y, K = K, level = .01, verbose = F, restrict = restrict, testType = "empLik"))[3]
  timeEL001 <- system.time(out001EL <- ngBap::bang(Y, K = K, level = .001, verbose = F, restrict = restrict, testType = "empLik"))[3]

  bangMod05EL <- BCD::ricf(out05EL$dEdge, out05EL$bEdge, t(dat$Y), tol = 1e-8)
  bangMod01EL <- BCD::ricf(out01EL$dEdge, out01EL$bEdge, t(dat$Y), tol = 1e-8)
  bangMod001EL <- BCD::ricf(out001EL$dEdge, out001EL$bEdge, t(dat$Y), tol = 1e-8)

  bangScore05EL = mean(mvtnorm::dmvnorm(dat$Y, sigma = bangMod05EL$SigmaHat, log = T)) - (sum(out05EL$dEdge) + (sum(out05EL$bEdge) - p)/2) * log(n) / n
  bangScore01EL = mean(mvtnorm::dmvnorm(dat$Y, sigma = bangMod01EL$SigmaHat, log = T)) - (sum(out01EL$dEdge) + (sum(out01EL$bEdge) - p)/2) * log(n) / n
  bangScore001EL = mean(mvtnorm::dmvnorm(dat$Y, sigma = bangMod001EL$SigmaHat, log = T)) - (sum(out001EL$dEdge) + (sum(out001EL$bEdge) - p)/2) * log(n) / n

  bangExact05EL = calc_shd(out05EL$bEdge, out05EL$dEdge, dat$bidirectEdges, dat$directEdges)
  bangExact01EL = calc_shd(out01EL$bEdge, out01EL$dEdge, dat$bidirectEdges, dat$directEdges)
  bangExact001EL = calc_shd(out001EL$bEdge, out001EL$dEdge, dat$bidirectEdges, dat$directEdges)

  timeGBS <- system.time(gbap_out <- greedyBAPs::greedySearch(cov(Y), n = n, n.restarts = 100))[3]

  gbapB <- t(ifelse(gbap_out$final.bap == 1, 1, 0))
  gbapOm <- ifelse(gbap_out$final.bap == 100, 1, 0) + diag(rep(1, p))


  trueMod <- BCD::ricf(dat$directEdges, dat$bidirectEdges, t(dat$Y), tol = 1e-8)
  gbapMod <- BCD::ricf(gbapB, gbapOm, t(dat$Y), tol = 1e-8)




  if(n < 5000) {

    timeDH05 <- system.time(out05DH <- ngBap::bang(Y, K = K, level = .05, verbose = F, restrict = restrict, testType = "dhsic"))[3]
    timeDH01 <- system.time(out01DH <- ngBap::bang(Y, K = K, level = .01, verbose = F, restrict = restrict, testType = "dhsic"))[3]
    timeDH001 <- system.time(out001DH <- ngBap::bang(Y, K = K, level = .001, verbose = F, restrict = restrict, testType = "dhsic"))[3]

    bangMod05DH <- BCD::ricf(out05DH$dEdge, out05DH$bEdge, t(dat$Y), tol = 1e-8)
    bangMod01DH <- BCD::ricf(out01DH$dEdge, out01DH$bEdge, t(dat$Y), tol = 1e-8)
    bangMod001DH <- BCD::ricf(out001DH$dEdge, out001DH$bEdge, t(dat$Y), tol = 1e-8)

    bangScore05DH = mean(mvtnorm::dmvnorm(dat$Y, sigma = bangMod05DH$SigmaHat, log = T)) - (sum(out05DH$dEdge) + (sum(out05DH$bEdge) - p)/2) * log(n) / n
    bangScore01DH = mean(mvtnorm::dmvnorm(dat$Y, sigma = bangMod01DH$SigmaHat, log = T)) - (sum(out01DH$dEdge) + (sum(out01DH$bEdge) - p)/2) * log(n) / n
    bangScore001DH = mean(mvtnorm::dmvnorm(dat$Y, sigma = bangMod001DH$SigmaHat, log = T)) - (sum(out001DH$dEdge) + (sum(out001DH$bEdge) - p)/2) * log(n) / n



    bangExact05DH = calc_shd(out05DH$bEdge, out05DH$dEdge, dat$bidirectEdges, dat$directEdges)
    bangExact01DH = calc_shd(out01DH$bEdge, out01DH$dEdge, dat$bidirectEdges, dat$directEdges)
    bangExact001DH = calc_shd(out001DH$bEdge, out001DH$dEdge, dat$bidirectEdges, dat$directEdges)




    ### Write data to disk, and then run RCD and pass adjacency matrix back to R
    write.csv(round(Y, 8), paste("temp/bap_obs_",runInd, "_", simInd, ".csv", sep = ""), row.names = F)

    system(paste('python3 rcd_sim3_call.py ', runInd, ' ', simInd, sep = ""), wait = T)

    rcd_adj_05 <- read.csv(paste("temp/rcd_bap05_", runInd, "_", simInd, ".csv", sep = ""), header=T)[, -1]
    rcd05_dEdge <- ifelse(is.na(rcd_adj_05), 0, rcd_adj_05 != 0)
    rcd05_bEdge <- is.na(rcd_adj_05) + diag(p)

    rcd_adj_01 <- read.csv(paste("temp/rcd_bap01_", runInd, "_", simInd, ".csv", sep = ""), header=T)[,-1]
    rcd01_dEdge <- ifelse(is.na(rcd_adj_01), 0, rcd_adj_01 != 0)
    rcd01_bEdge <- is.na(rcd_adj_01) + diag(p)

    rcd_adj_001 <- read.csv(paste("temp/rcd_bap001_", runInd, "_", simInd, ".csv", sep = ""), header=T)[,-1]
    rcd001_dEdge <- ifelse(is.na(rcd_adj_001), 0, rcd_adj_001 != 0)
    rcd001_bEdge <- is.na(rcd_adj_001) + diag(p)

    rcd_time <- c(unlist(read.csv(paste("temp/rcd_bap_time_", runInd, "_", simInd, ".csv", sep = ""), header=T)[,-1]))




    rcdMod05 <- BCD::ricf(rcd05_dEdge, rcd05_bEdge, t(dat$Y), tol = 1e-8)
    rcdScore05 = mean(mvtnorm::dmvnorm(dat$Y, sigma = rcdMod05$SigmaHat, log = T)) - (sum(rcd05_dEdge) + (sum(rcd05_bEdge) - p)/2) * log(n) / n
    rcdExact05 = calc_shd(rcd05_bEdge, rcd05_dEdge, dat$bidirectEdges, dat$directEdges)

    rcdMod01 <- BCD::ricf(rcd01_dEdge, rcd01_bEdge, t(dat$Y), tol = 1e-8)
    rcdScore01 = mean(mvtnorm::dmvnorm(dat$Y, sigma = rcdMod01$SigmaHat, log = T)) - (sum(rcd01_dEdge) + (sum(rcd01_bEdge) - p)/2) * log(n) / n
    rcdExact01 = calc_shd(rcd01_bEdge, rcd01_dEdge, dat$bidirectEdges, dat$directEdges)

    rcdMod001 <- BCD::ricf(rcd001_dEdge, rcd001_bEdge, t(dat$Y), tol = 1e-8)
    rcdScore001 = mean(mvtnorm::dmvnorm(dat$Y, sigma = rcdMod001$SigmaHat, log = T)) - (sum(rcd001_dEdge) + (sum(rcd001_bEdge) - p)/2) * log(n) / n
    rcdExact001 = calc_shd(rcd001_bEdge, rcd001_dEdge, dat$bidirectEdges, dat$directEdges)


    system(paste("rm temp/*bap*_", runInd, "_", simInd, ".csv", sep = ""))



  } else {

    bangScore01DH = NA
    bangScore001DH = NA
    bangScore05DH = NA

    bangExact01DH = NA
    bangExact001DH = NA
    bangExact05DH = NA
    timeDH05 <- NA; timeDH01 <- NA; timeDH001 <- NA;

    rcdScore05 = NA
    rcdScore01 = NA
    rcdScore001 = NA

    rcdExact05 = NA
    rcdExact01 = NA
    rcdExact001 = NA
    rcd_time <- rep(NA, 3)

  }






  ret <- c(trueScore = mean(mvtnorm::dmvnorm(dat$Y, sigma = trueMod$SigmaHat, log = T)) - (sum(dat$directEdges) + (sum(dat$bidirectEdges) - p)/2) * log(n) / n,
           gbsScore = mean(mvtnorm::dmvnorm(dat$Y, sigma = gbapMod$SigmaHat, log = T)) - (sum(gbapB) + (sum(gbapOm) - p)/2) * log(n) / n,

           bangScore05EL = bangScore05EL,
           bangScore01EL = bangScore01EL,
           bangScore001EL = bangScore001EL,

           bangExact05EL = bangExact05EL,
           bangExact01EL = bangExact01EL,
           bangExact001EL = bangExact001EL,

          bangScore05DH = bangScore05DH,
          bangScore01DH = bangScore01DH,
          bangScore001DH = bangScore001DH,

          bangExact05DH = bangExact05DH,
          bangExact01DH = bangExact01DH,
          bangExact001DH = bangExact001DH,

          rcdScore05 = rcdScore05,
          rcdScore01 = rcdScore05,
          rcdScore001 = rcdScore001,

          rcdExact05 = rcdExact05,
          rcdExact01 = rcdExact01,
          rcdExact001 = rcdExact001,
          timeGBS= timeGBS,
          timeEL05 = timeEL05, timeEL01 = timeEL01, timeEL001 = timeEL001,
          timeDH05 = timeDH05, timeDH01 = timeDH01, timeDH001 = timeDH001,
          timeRCD05 = rcd_time[1], timeRCD01 =  rcd_time[2], timeRCD001 =  rcd_time[3])

  return(ret)

}



##########################################################
library(parallel)
sim.size <- 250
clust.size <- 7

p <- 6
### param.grid size: 288 ###
param.grid <- expand.grid(n = c(500, 1000, 2500, 10000, 25000, 50000), dist = c("gamma", "unif", "t", "lognormal"),
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


out <- parSapply(cl, 1:sim.size, function(i){one.run.baps(p = p, n = n, dist = dist, K = K, d = d, b = b, restrict = 1, signs = signs, simInd = i, runInd = runInd)})

out <- t(out)



stopCluster(cl)

cover <- cbind(cover, data.matrix(out))
cover <- data.frame(apply(cover, MAR = 2, unlist))

write.table(cover, paste("../resultsNew/bap_",runInd,".csv", sep = ""), row.names = F, sep = ",")

####
