### Evaluate a the value of a polynomial of error moments
# given a polynomial of "eps" terms created using mpoly
# evaluates the polynomial using errMomentList
# where the first p columns of errMomentList
# correspond to the exponent of each term
# and the last column corresponds to the value

evalMoment <- function(term, errMomentList){
  p <- ncol(errMomentList)-1
  if(length(term) > 1){
    return(sum(sapply(unclass(term), evalHelper, errMomentList)))
  } else {
    return(evalHelper(unclass(term)[[1]], errMomentList))
  }

}

### Helper to evaluate a the value of a polynomial of error moments
# given a monomial of "eps" terms created using mpoly
# evaluates the polynomial using errMomentList
# where the first p columns of errMomentList
# correspond to the exponent of each term
# and the last column corresponds to the value
evalHelper <- function(monomial, errMomentList){
  p <- ncol(errMomentList) - 1

  # number of terms in monomial
  numVars <- length(monomial) - 1
  # indices of those terms
  ind <- as.numeric(gsub("eps", "", names(monomial)[1:numVars]))
  deg <- rep(0, p)
  deg[ind] <- monomial[1:numVars]
  # pull moment evaluation from table
  momentEval <- errMomentList[which(apply(errMomentList[,1:p], MAR = 1, function(x){all(deg == x)} )), p+1]
  # multiply by coefficient and return
  return(monomial[numVars + 1] * momentEval)
}

# create an mp with eps
# wrapper for mpoly::mp
eps <- function(i){
  return(mpoly::mp(paste("eps", i, sep = "")))
}

# Generate Data given some linear SEM
# Given a matrix of linear coefficients, creates a list of symbolic Y variables
# from the "eps" terms
generate_data <- function(B){

  totalEffect <- round(solve(diag(p) - B), 12)
  ordering <- rev(ngBap::getOrdering(totalEffect))

  p <- ncol(B)
  varList <- vector("list", p)


  varList[[ ordering[1] ]] <- ngBap::eps( ordering[1] )

  for(i in 2:p){
    varList[[ ordering[i] ]] <- ngBap::eps( ordering[i] )

    for(j in 1:i){
     varList[[ ordering[i] ]] <- varList[[ ordering[i] ]] + as.numeric(B[ordering[i], ordering[j] ]) * varList[[ ordering[j] ]]
    }
  }

  return(varList)
}



## Generate up to K-degree moments for a p-variate distribution
# Given a covariance matrix omega, it starts with the moments of a multivariate
# Gaussian, then it peturbs the marginal moments in the noiseMoments vector
# by +/- U(noiseLow, noiseHigh); i.e., if noiseMoments = c(3), it adds uniform
# noise to each of the eps1^3, eps2^3, etc
generateMoments <- function(p, K, omega, noiseMoments, noiseLow, noiseHigh){

  z <- list()
  for(i in 1:p){
    z[[i]] <- c(0:K)
  }

  momentInd <- expand.grid(z)
  momentInd <- momentInd[-1, ]
  momentInd <- as.matrix(momentInd[which(rowSums(momentInd) <= K), ])

  colnames(momentInd) <- paste("eps", 1:p, sep = "")

  omegaVec <-  c(omega[lower.tri(omega, diag = T)])

  evalMoments <- rep(0, nrow(momentInd))
  for(i in 1:nrow(momentInd)){
    if(sum(momentInd[i,]) %%2 == 1){
      evalMoments[i] <- 0
    } else {
      evalMoments[i] <- symmoments::evaluate(symmoments::callmultmoments(momentInd[i,]),omegaVec)
    }
  }


  for(j in noiseMoments){
    ind <- which(apply(momentInd, MAR = 1, function(x){sum(x==0) == p-1 & max(x) == j}))
    evalMoments[ind] <- evalMoments[ind] + sample(c(-1, 1), size = length(ind), replace = T) * runif(length(ind), noiseLow, noiseHigh)

  }



  momentInd <- cbind(momentInd, evalMoments)
  return(momentInd)

}


