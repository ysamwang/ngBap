runInd <- 1
args <- commandArgs(TRUE)
for(i in 1:length(args)){
  eval(parse(text = args[[i]]))
}
print(runInd)



sim.size <- 50
p <- 10

rec <- matrix(0, nrow = sim.size, ncol = 2)
for(i in 1:sim.size){
  cat("Sim ")
  cat(i)
  cat(": ")
  bap <- greedyBAPs::generateUniformBAP(p = p, N = 1)
  directedEdges <- (bap[[1]] == 1) + 0
  bidirectedEdges <- (bap[[1]] == 100) + diag(p)

  B <- sample(c(-1, 1), size = p^2, replace = T) * runif(p^2, 0, 1) * directedEdges

  Omega <- sample(c(-1, 1), size = p^2, replace = T) * runif(p^2, 0, .5) * bidirectedEdges
  Omega <- Omega + t(Omega)
  Omega <- Omega - diag(diag(Omega))
  # ensure that Omega.edge is PD
  diag(Omega) <- rowSums(abs(Omega)) + 1






  Sigma <- solve(diag(p) - B, Omega) %*% t(solve(diag(p) - B))
  cat("Fitting Model... ")

  Y <- ngBap::generate_data(B)
  errMomentList <- ngBap::generateMoments(p = p, K = 3, Omega, noiseMoments = c(3), .4, .8)

  oracle_est <- ngBap::bang_symb(Y, Sigma, errMomentList, K = 3, cutoff = 1e-8, verbose = F)
  oracle_est2 <- ngBap::bang_symb2(Y, Sigma, errMomentList, K = 3, cutoff = 1e-8, verbose = F)

  rec[i, 1 ] <- sum(abs(oracle_est$dEdge - directedEdges)) + sum(abs(oracle_est$bEdge - bidirectedEdges))
  rec[i, 2 ] <- sum(abs(oracle_est2$dEdge - directedEdges)) + sum(abs(oracle_est2$bEdge - bidirectedEdges))
  cat("Discrepancy: ")
  cat(rec[i, ])
  cat("\n")

  if(any(rec[i, ] > 0)){

    write.csv(B, paste("temp/oracle_B_", runInd, "_", i, ".csv", sep = ""), row.names = F)
    write.csv(Omega, paste("temp/oracle_Om_", runInd, "_", i, ".csv", sep = ""), row.names = F)
    write.csv(errMomentList, paste("temp/oracle_err_", runInd, "_", i, ".csv", sep = ""), row.names = F)

    cat("\n")
    print(directedEdges)
    print(oracle_est$dEdge)
    cat("\n")
    print(bidirectedEdges)
    print(oracle_est$bEdge)

  }

}

write.csv(rec, paste("../resultsRev/oracle/oracle_", runInd, ".csv", sep = ","), row.names = F)






###############
p <- 10
runInd <- 3
i <- 18
B <- as.matrix(read.csv(paste("temp/oracle_B_", runInd, "_", i, ".csv", sep = "")))
Omega <- as.matrix(read.csv(paste("temp/oracle_Om_", runInd, "_", i, ".csv", sep = "")))
errMomentList <- as.matrix(read.csv(paste("temp/oracle_err_", runInd, "_", i, ".csv", sep = "")))

dEdge <- (round(B, 12) != 0) + 0
bEdge <- (round(Omega, 12) != 0) + 0
ngBap::getOrdering(solve(diag(p) - B))

Sigma <- solve(diag(p) - B, Omega) %*% t(solve(diag(p) - B))
Y <- ngBap::generate_data(B)



oracle_est <- ngBap::bang_symb(Y, Sigma, errMomentList, K = 3, cutoff = 1e-12, verbose = F)
(oracle_est$dEdge - dEdge)
(oracle_est$bEdge - bEdge)

oracle_est2 <- ngBap::bang_symb2(Y, Sigma, errMomentList, K = 3, cutoff = 1e-12, verbose = F)

(oracle_est2$dEdge - dEdge)
(oracle_est2$bEdge - bEdge)
