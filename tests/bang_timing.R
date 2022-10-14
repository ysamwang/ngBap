runInd <- 1
args <- commandArgs(TRUE)
for(i in 1:length(args)){
  eval(parse(text = args[[i]]))
}
print(runInd)

library("ngBap")
library("pcalg")
library("greedyBAPs")
library("dHSIC")

n.list <- c(500, 1000, 2000, 10000, 25000, 50000)
p.list <- c(6, 8, 10)

param.grid <- expand.grid(n.list, p.list)
dist <- "gamma"
signs <- T
K <- 3

final.rec <- matrix(0, nrow = nrow(param.grid), ncol = 2)

run.once <- function(n, p, dist, d, b){
  rec <- rep(0,2)
  dat <- ngBap::rBAP(n = n, p = p, dist = dist, d = d, b = b, ancestral = F, shuffle = T, signs = T)
  Y <- dat$Y

  rec[1] <- system.time(out001EL <- ngBap::bang(Y, K = K, level = .01,
                                                   verbose = F, restrict = 1, testType = "empLik"))[3]

  rec[2] <- system.time(gbap_out <- greedyBAPs::greedySearch(cov(Y), n = n, n.restarts = 100))[3]
  return(rec)
}

sim.size <- 200
clust.size <- 15
cl <- makeCluster(clust.size)
clusterExport(cl, ls())


n <- param.grid[runInd, 1]
p <- param.grid[runInd, 2]
d <- p
b <- p

rec <- parSapply(cl, 1:sim.size, function(x){run.once(n, p, dist, d, b)})

finalTab <- cbind(param.grid, final.rec)
stopCluster(cl)

write.table(finalTab, paste("../resultsNew/timingTab_",runInd,".csv", sep = ""), row.names = F, sep = ",")



