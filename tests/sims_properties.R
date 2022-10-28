one.run.baps <- function(regime = "sparse"){

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

  dat <- ngBap::rBAP(n = 2, p = 6, dist = "gamma", d = d, b = b, ancestral = F, shuffle = T, signs = T)

  ancest <- (solve(diag(p) - dat$B) != 0 ) - diag(p)
  sibs <- dat$bidirectEdges - diag(p)


  return(sum(ancest * sibs))

}

###
### param.grid size: 288 ###
regime.list <- c("sparse", "medium", "dense")
i <- 1
rec <- matrix(0, ncol = 2, nrow = 3)
for(i in 1:3){
  regime <- regime.list[i]
  out <- replicate(10000, one.run.baps(regime = regime))
  rec[i,] <- c(mean(out == 0), mean(out))
}

rec


