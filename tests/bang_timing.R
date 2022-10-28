runInd <- 1
args <- commandArgs(TRUE)
for(i in 1:length(args)){
  eval(parse(text = args[[i]]))
}
print(runInd)

library("parallel")
library("ngBap")
library("pcalg")
library("greedyBAPs")
library("dHSIC")

sim.size <- 200
clust.size <- 15

n.list <- c(500, 1000, 2000, 10000, 25000, 50000)
p.list <- c(6, 8, 10)

param.grid <- expand.grid(n.list, p.list)
dim(param.grid)

dist <- "gamma"
signs <- T
K <- 3


run.once <- function(n, p, dist, d, b){
  rec <- rep(0,2)
  dat <- ngBap::rBAP(n = n, p = p, dist = dist, d = d, b = b, ancestral = F, shuffle = T, signs = T)
  Y <- dat$Y

  rec[1] <- system.time(out001EL <- ngBap::bang(Y, K = K, level = .01,
                                                   verbose = F, restrict = 1, testType = "empLik"))[3]

  rec[2] <- system.time(gbap_out <- greedyBAPs::greedySearch(cov(Y), n = n, n.restarts = 100))[3]
  return(rec)
}





n <- param.grid[runInd, 1]
p <- param.grid[runInd, 2]
d <- p
b <- p

cl <- makeCluster(3)
clusterExport(cl, ls())

rec <- parSapply(cl, 1:sim.size, function(x){run.once(n, p, dist, d, b)})

finalTab <- cbind(param.grid, rec)
stopCluster(cl)

write.table(finalTab, paste("../resultsNew/timingTab_",runInd,".csv", sep = ""), row.names = F, sep = ",")


###################
n.list <- c(500, 1000, 2000, 10000, 25000, 50000)
p.list <- c(6, 8, 10)

param.grid <- expand.grid(n.list, p.list)
dim(param.grid)

runInd <- 1
dat <- read.csv(paste("~/Dropbox/mixedGraphID/testOutputs/resultsRev/timing/timingTab_",runInd,".csv", sep = ""))

for(runInd in 2:nrow(param.grid)){
  dat <- rbind(dat,read.csv(paste("~/Dropbox/mixedGraphID/testOutputs/resultsRev/timing/timingTab_",runInd,".csv", sep = "")))
}

out <- cbind(aggregate(BANG_time ~ n + p, data = dat, FUN = mean),
      aggregate(GBS_time ~ n + p, data = dat, FUN = mean)[,3])
colnames(out) <- c("n", "p", "BANG", "GBS")

bang_mat <- matrix(out$BANG, nrow = 6)
gbs_mat <- matrix(out$GBS, nrow = 6)
comb_mat <- matrix(0, nrow = 6, ncol = ncol(bang_mat) * 2)

for(i in 1:ncol(bang_mat)){
  comb_mat[, (2 * i - 1):(2*i)] <- cbind(bang_mat[, i], gbs_mat[,i])
}
rownames(comb_mat) <- n.list
colnames(comb_mat) <- rep(c("BANG", "GBS"), 3)
print(xtable::xtable(comb_mat, digits = 1))



