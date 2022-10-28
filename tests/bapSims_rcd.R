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
library("reticulate")





#######
# 10/14/22: Run simulations which use RCD (from python lingam module)


one.run.baps <- function(p, n,  dist, K = 4, d, b, restrict = 1, signs = T){
  dat <- ngBap::rBAP(n = n, p = p, dist = dist, d = d, b = b, ancestral = F, shuffle = T, signs = signs)

  Y <- dat$Y





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

sim.size <- 250
p <- 6


### param.grid size: 288 ###
param.grid <- expand.grid(n = c(500, 1000, 1500, 10000, 25000, 50000), dist = c("gamma", "unif", "t", "lognormal"),
                          regime = c("sparse", "medium", "dense"), signs = c(T, F))
param.grid <- rbind(param.grid, param.grid)



i <- 1
n <- param.grid[i, 1]
dist <- as.character(param.grid[i, 2])
regime <- as.character(param.grid[i, 3])
signs <- param.grid[i, 4]
cover <- data.frame(p = rep(p, sim.size), n = rep(n, each = sim.size), dist = rep(dist, each = sim.size),
                    regime = rep(regime, sim.size), signs = rep(signs, sim.size))
dat <- cbind(cover, read.csv(paste("~/Dropbox/mixedGraphID/testOutputs/resultsNew/rcd_",i,".csv", sep = ""), sep = ",")[,-1])
colnames(dat) <- c("p", "n", "dist", "regime", "signs", "trueScore", "estScore", "missEdges", "time")

for(i in 2:nrow(param.grid)){
  n <- param.grid[i, 1]
  dist <- as.character(param.grid[i, 2])
  regime <- as.character(param.grid[i, 3])
  signs <- param.grid[i, 4]
  cover <- data.frame(p = rep(p, sim.size), n = rep(n, each = sim.size), dist = rep(dist, each = sim.size),
                      regime = rep(regime, sim.size), signs = rep(signs, sim.size))

  if(file.exists(paste("~/Dropbox/mixedGraphID/testOutputs/resultsNew/rcd_",i,".csv", sep = ""))){
    temp <- cbind(cover, read.csv(paste("~/Dropbox/mixedGraphID/testOutputs/resultsNew/rcd_",i,".csv", sep = ""), sep = ",")[,-1])
    colnames(temp) <- c("p", "n", "dist", "regime", "signs", "trueScore", "estScore", "missEdges", "time")
    dat <- rbind(dat, temp)
  }


}

aggregate(missEdges == 0 ~ n + dist + regime + signs, data = dat, FUN = mean)




