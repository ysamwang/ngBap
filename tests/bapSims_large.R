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
# 10/15/22: Setting with p = 16
#

one.run.baps <- function(p, n,  dist, K = 4, d, b, restrict = 1, signs = T){
  dat <- ngBap::rBAP(n = n, p = p, dist = dist, d = d, b = b, ancestral = F, shuffle = T, signs = signs)

  Y <- dat$Y

  t1 <- system.time(out05EL <- ngBap::bang(Y, K = K, level = .05, verbose = F, restrict = restrict, testType = "empLik"))
  t2 <- system.time(out01EL <- ngBap::bang(Y, K = K, level = .01, verbose = F, restrict = restrict, testType = "empLik"))
  t3 <- system.time(out001EL <- ngBap::bang(Y, K = K, level = .001, verbose = F, restrict = restrict, testType = "empLik"))


  ret <- matrix(
    c(sum(abs(out05EL$bEdge - dat$bidirectEdges))/2, sum(abs(out05EL$dEdge - dat$directEdges)), t1[3],
      sum(abs(out01EL$bEdge - dat$bidirectEdges))/2, sum(abs(out01EL$dEdge - dat$directEdges)), t2[3],
      sum(abs(out001EL$bEdge - dat$bidirectEdges))/2, sum(abs(out001EL$dEdge - dat$directEdges)), t3[3]), nrow = 1)

  colnames(ret) <- paste(rep(c("b", "d", "time"), times = 3), rep(c("05", "01", "001"), each = 3), sep = "_")

  # ret <- list(bangExact05EL =
  #             # bangExact01EL = sum(abs(out01EL$bEdge - dat$bidirectEdges) + abs(out01EL$dEdge - dat$directEdges)) == 0,
  #             # bangExact001EL = sum(abs(out001EL$bEdge - dat$bidirectEdges) + abs(out001EL$dEdge - dat$directEdges)) == 0
  #             )

  return(ret)

}



##########################################################
library(parallel)
sim.size <- 125
clust.size <- 15
p <- 16


### param.grid size: 384 ###
param.grid <- expand.grid(n = c(1500, 10000, 25000, 50000), dist = c("gamma", "unif", "t", "lognormal"),
                          regime = c("sparse", "medium", "dense"), signs = c(T, F))
param.grid <- rbind(param.grid, param.grid, param.grid, param.grid)




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
  b <- p * 1.5

}

cl <- makeCluster(clust.size)
clusterExport(cl, ls())


out <- t(parSapply(cl, 1:sim.size, function(i){one.run.baps(p = p, n = n, dist = dist, K = K, d = d, b = b, restrict = 1, signs = signs)}))
colnames(out) <- paste(rep(c("b", "d", "time"), times = 3), rep(c("05", "01", "001"), each = 3), sep = "_")

stopCluster(cl)

cover <- cbind(cover, out)
cover <- data.frame(apply(cover, MAR = 2, unlist))

write.table(cover, paste("../resultsRev/bapLarge/bapLarge_",runInd,".csv", sep = ""), row.names = F, sep = ",")



### param.grid size: 288 ###
# param.grid <- expand.grid(n = c(1500, 10000, 25000, 50000), dist = c("gamma", "unif", "t", "lognormal"),
#                           regime = c("sparse", "medium", "dense"), signs = c(T, F))
# param.grid <- rbind(param.grid, param.grid, param.grid, param.grid)
#
# missing <- c()
# runInd <- 1
# dat <- read.csv(paste("~/Dropbox/mixedGraphID/testOutputs/resultsRev/bapLarge/bapLarge_",runInd,".csv", sep = ""),sep = ",")
# for(runInd in 2:nrow(param.grid)){
#   if(file.exists(paste("~/Dropbox/mixedGraphID/testOutputs/resultsRev/bapLarge/bapLarge_",runInd,".csv", sep = ""))){
#     dat <- rbind(dat,
#                  read.csv(paste("~/Dropbox/mixedGraphID/testOutputs/resultsRev/bapLarge/bapLarge_",runInd,".csv", sep = ""),sep = ","))
#
#   } else{
#     missing <- c(missing, runInd)
#   }
# }

