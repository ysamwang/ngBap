
### Symbolic evaluation
# Prune siblings
# Returns updated siblings list
#
pruneSib_symb <- function(v, Y.current, sibSet, errMomentList, cutoff = 1e-8){
    sibs <- which(sibSet[v,] != 0)

    for(s in sibs){
      temp <- Y.current[[v]] * Y.current[[s]]

      if(abs(evalMoment(temp, errMomentList)) < cutoff){
        sibSet[v,s] <- sibSet[s,v] <- 0
      }

    }
  return(sibSet)
}




# Certify if C are parents of v
# Returns true if C can be certified
# returns false if C cannot be certified
certifyAncest_symb <- function(C, v, directEffect, Sigma, Y.orig, Y.current, K, errMomentList, cutoff = 1e-8, verbose = F){

  p <- ncol(directEffect)

  totalEffect <- round(solve(diag(p) - directEffect), 12)
  A <- unique(c(unlist(apply(totalEffect[C, , drop = F], MAR = 1, function(x){which(x != 0)}))))
  l1 <- (diag(p) - directEffect)[C,A]

  delta <- round(solve(l1 %*% Sigma[A, C]) %*% l1 %*% Sigma[A,v], 12)

  gamma_v <- Y.orig[[v]]
  for(j in 1:length(C)){
    gamma_v <- gamma_v - delta[j] * Y.orig[[ C[j] ]]
  }

  momentsToTest <- rep(0, length(C))
  for(j in 1:length(C)){
    temp <- Y.current[[ C[j] ]]^(K-1) * gamma_v
    momentsToTest[j] <- evalMoment(temp, errMomentList)
  }

  if(verbose){

    cat(formatC(max(abs(momentsToTest)), format = "e", digits = 2))

  }
  return(all(abs(momentsToTest) < cutoff))

}



updateDebiased_symb <- function(v, C, Y.orig, Y.current, Sigma, directEffect, verbose = F){

  p <- ncol(directEffect)
  totalEffect <- round(solve(diag(p) - directEffect), 12)


  if(length(C) > 0){



    A <- unique(c(unlist(apply(totalEffect[C, , drop = F], MAR = 1, function(x){which(x != 0)}))))
    l1 <- (diag(p) - directEffect)[C,A]

    beta.debiased <-  round(solve(l1 %*% Sigma[A, C]) %*% l1 %*% Sigma[A,v], 12)

    directEffect[v, C] <- beta.debiased
    directEffect[v, -C] <- 0
    totalEffect <- round(solve(diag(rep(1, p)) - directEffect), 12)

    gamma_v <- Y.orig[[v]]
    for(j in 1:length(C)){
      gamma_v <- gamma_v - beta.debiased[j] * Y.orig[[ C[j] ]]
    }

    return(list(totalEffect = totalEffect, directEffect = directEffect, res.v = gamma_v))

  } else {

    return(list(totalEffect = totalEffect, directEffect = directEffect, res.v = Y.current[[v]] ))

  }
}

# Certify if C are parents of v
# Returns T if u can be pruned
pruneAncest_symb <- function(u, C, v, directEffect, Sigma, Y.orig, Y.current, K, errMomentList, cutoff = 1e-8){

  p <- ncol(directEffect)
  totalEffect <- round(solve(diag(p) - directEffect), 12)

  # Set of regressors to test
  # parents with u removed
  C.test <- setdiff(C, u)

  if(length(C.test) > 0){
    A <- unique(c(unlist(apply(totalEffect[C.test, , drop = F], MAR = 1, function(x){which(x != 0)}))))
    l1 <- (diag(p) - directEffect)[C.test,A]

    delta <- round(solve(l1 %*% Sigma[A, C.test]) %*% l1 %*% Sigma[A,v], 12)

    gamma_v <- Y.orig[[v]]
    for(j in 1:length(C.test)){
      gamma_v <- gamma_v - delta[j] * Y.orig[[ C.test[j] ]]
    }
  } else {

    gamma_v <- Y.orig[[v]]

  }



  momentsToTest <- rep(0, length(C))
  for(j in 1:length(C)){
    temp <- Y.current[[ C[j] ]]^(K-1) * gamma_v
    momentsToTest[j] <- evalMoment(temp, errMomentList)
  }

  return(all(abs(momentsToTest) < cutoff))

}



bang_symb <- function(Y, Sigma, errMomentList, K, cutoff = 1e-8, verbose = F){

  if(verbose){
    cat("======= Estimating BAP ========\n")
  }


  ### directEffect is B matrix
  ### totalEffect is (I-B)^{-1}
  totalEffect <- directEffect <- matrix(0, nrow = p, ncol = p)
  diag(totalEffect) <- 1

  siblings <- matrix(1, nrow = p, ncol = p)
  diag(siblings) <- 0

  Y.orig <- Y


  # counter for set size
  l <- 1

  while(any(rowSums(siblings) >= l)){
    old.directEffect <- directEffect




    # cycle through variables
    for (v in 1:p) {
      poss.pa <- which(siblings[v, ] != 0)

      if (verbose) {
        cat(paste("\n\n=== Checking ", v, " ===\n"))
        cat("Remaining siblings: ")
        cat(paste(poss.pa, collapse = " "))
        cat(paste("; sets of size: ", l, "\n"))
      }

      # prune anything that is independent
      if(length(poss.pa) > 0){

        ## Call to symbolic eval
        # takes in current Y, siblings
        # returns updated siblings matrix
        #
        siblings <- pruneSib_symb(v, Y, siblings, errMomentList, cutoff)
        if(verbose){
          cat("Unpruned siblings: ")
          cat(which(siblings[v,] != 0))
          cat("\n")
        }


      } # End independence testing


      # Set the possible parents to be any nodes which are not certified and are not ancestors of v
      poss.pa <- setdiff(which(siblings[v, ] != 0), which(totalEffect[, v] != 0))

      if(length(poss.pa) >= l){



        # See if we've already checked all the remaining sets
        if(length(poss.pa) == 1){
          # if only 1 possible parent
          # combn gives unwanted result, so manually make matrix
          condSets <- matrix(c(poss.pa,
                               which(directEffect[v, ] != 0)), ncol = 1)

        } else {

          ## Get all sets of possible parents of size l
          ## Add current parents
          condSets <- combn(poss.pa, l)
          condSets <- rbind(condSets, matrix(which(directEffect[v, ] != 0),
                                             nrow = length(which(directEffect[v, ] != 0)),
                                             ncol = dim(condSets)[2]))

        }

          ### Call to Symbolic Evaluation
          # certifyAncest takes in C, v, D, Sigma,

        j <- 0
        foundCertified <- FALSE


        while(j < ncol(condSets) & !foundCertified){
          j <- j + 1




          if(verbose){
            cat(paste(condSets[, j], collapse = "-"))
            cat(">")
          }

          foundCertified <- certifyAncest_symb(condSets[, j], v, directEffect, Sigma,
                                               Y.orig, Y, K, errMomentList, cutoff, verbose = verbose)
          if(verbose){
            cat(" ")

          }

        }

        if(foundCertified){
          certifiedParents <- condSets[, j]
          if(verbose){
            cat("Certified: ")
            cat(paste(certifiedParents, collapse = "-"))
            cat("\n")

          }

        } else {
          certifiedParents <- c()
        }




        if(length(certifiedParents) > 0){
          # Can be certified as unconfounded ancestor so cannot be a sibling
          siblings[v, certifiedParents] <- siblings[certifiedParents, v] <- 0


          # Update effect estimates

          update <- updateDebiased_symb(v, certifiedParents, Y.orig, Y, Sigma, directEffect)


          totalEffect <- update$totalEffect
          directEffect <- update$directEffect
          Y[[v]] <- update$res.v

          if(verbose){
            cat("\n")
            print(round(directEffect, 3))
          }

        }




      } # End if remaining siblings >= l for v
    } # end cycle through 1:p

    # If any updates were made to Y
    if(sum((old.directEffect - directEffect)^2) > 1e-8 ) {
      # reset counter to 1
      l <- 1
    } else {
      # else advance counter
      l <- l + 1
    }
  } # End big while for ANY sib(v) >= l



  #### Pruning Procedure #####
  ordering <- getOrdering(totalEffect)
  if(verbose){
    print("!!Pruning Procedure!!")
    print(directEffect)

    print(ordering)
  }

  for(z in ordering){
    parents <- which(directEffect[z, ] != 0)

    if( length(parents) > 0){
      if(verbose){
        cat(paste("Pruning:", z))
       cat("\nCurrent Parents: ")
       cat(parents)

      }


      pruningStat <- sapply(parents, pruneAncest_symb, parents, z, directEffect, Sigma, Y.orig, Y, K, errMomentList, cutoff)
      if(verbose){
        cat(pruningStat)

      }

      certifiedParents <- parents[which(!pruningStat)]



      if(length(certifiedParents) > 0){

        update <- updateDebiased_symb(z, certifiedParents, Y.orig, Y, Sigma, directEffect)

        totalEffect <- update$totalEffect
        directEffect <- update$directEffect
        Y[[z]] <- update$res.v

      } else {

        # If no parents, remove all from directEffect and update totalEffect
        directEffect[z, ] <- 0
        totalEffect <- solve(diag(rep(1, p)) -  directEffect)
        Y[[z]] <- Y.orig[[z]]
      }



    } else {
      if(verbose){
        cat(paste("Pruning:", z, "; no parents\n"))
      }
    }






  }


  return(list(totalEffect = totalEffect,
              directEffect = directEffect, dEdge = ifelse(directEffect != 0, 1, 0),
              bEdge = siblings + diag(rep(1, p))))
}





bang_symb2 <- function(Y, Sigma, errMomentList, K, cutoff = 1e-8, verbose = F){

  if(verbose){
    cat("======= Estimating BAP ========\n")
  }


  ### directEffect is B matrix
  ### totalEffect is (I-B)^{-1}
  totalEffect <- directEffect <- matrix(0, nrow = p, ncol = p)
  diag(totalEffect) <- 1

  siblings <- matrix(1, nrow = p, ncol = p)
  diag(siblings) <- 0

  Y.orig <- Y


  # counter for set size
  l <- 1

  while(any(rowSums(siblings) >= l)){
    old.directEffect <- directEffect




    # cycle through variables
    for (v in 1:p) {
      poss.pa <- which(siblings[v, ] != 0)

      if (verbose) {
        cat(paste("\n\n=== Checking ", v, " ===\n"))
        cat("Remaining siblings: ")
        cat(paste(poss.pa, collapse = " "))
        cat(paste("; sets of size: ", l, "\n"))
      }

      # prune anything that is independent
      if(length(poss.pa) > 0){

        ## Call to symbolic eval
        # takes in current Y, siblings
        # returns updated siblings matrix
        #
        siblings <- pruneSib_symb(v, Y, siblings, errMomentList, cutoff)
        if(verbose){
          cat("Unpruned siblings: ")
          cat(which(siblings[v,] != 0))
          cat("\n")
        }


      } # End independence testing


      # Set the possible parents to be any nodes which are not certified and are not ancestors of v
      poss.pa <- setdiff(which(siblings[v, ] != 0), which(totalEffect[, v] != 0))

      if(length(poss.pa) >= l){



        # See if we've already checked all the remaining sets
        if(length(poss.pa) == 1){
          # if only 1 possible parent
          # combn gives unwanted result, so manually make matrix
          condSets <- matrix(c(poss.pa,
                               which(directEffect[v, ] != 0)), ncol = 1)

        } else {

          ## Get all sets of possible parents of size l
          ## Add current parents
          condSets <- combn(poss.pa, l)
          condSets <- rbind(condSets, matrix(which(directEffect[v, ] != 0),
                                             nrow = length(which(directEffect[v, ] != 0)),
                                             ncol = dim(condSets)[2]))

        }

        ### Call to Symbolic Evaluation
        # certifyAncest takes in C, v, D, Sigma,
        ancestralCert <- apply(condSets, MAR = 2, certifyAncest_symb, v, directEffect, Sigma,
                               Y.orig, Y, K, errMomentList, cutoff)

        certifiedParents <- unique(c(condSets[, which(ancestralCert)]))

        if(verbose){
          cat("Testing Parents: ")
          cat(poss.pa)
          cat("\n")

          cat("Certified Parents:" )
          cat(certifiedParents)




        }

        if(length(certifiedParents) > 0){
          # Can be certified as unconfounded ancestor so cannot be a sibling
          siblings[v, certifiedParents] <- siblings[certifiedParents, v] <- 0


          # Update effect estimates

          update <- updateDebiased_symb(v, certifiedParents, Y.orig, Y, Sigma, directEffect)


          totalEffect <- update$totalEffect
          directEffect <- update$directEffect
          Y[[v]] <- update$res.v

          if(verbose){
            cat("\n")
            print(round(directEffect, 3))
          }

        }




      } # End if remaining siblings >= l for v
    } # end cycle through 1:p

    # If any updates were made to Y
    if(sum((old.directEffect - directEffect)^2) > 1e-8 ) {
      # reset counter to 1
      l <- 1
    } else {
      # else advance counter
      l <- l + 1
    }
  } # End big while for ANY sib(v) >= l



  #### Pruning Procedure #####
  ordering <- getOrdering(totalEffect)
  if(verbose){
    print("!!Pruning Procedure!!")
    print(directEffect)

    print(ordering)
  }

  for(z in ordering){
    parents <- which(directEffect[z, ] != 0)

    if( length(parents) > 0){
      if(verbose){
        cat(paste("Pruning:", z))
        cat("\nCurrent Parents: ")
        cat(parents)

      }


      pruningStat <- sapply(parents, pruneAncest_symb, parents, z, directEffect, Sigma, Y.orig, Y, K, errMomentList, cutoff)
      if(verbose){
        cat(pruningStat)

      }

      certifiedParents <- parents[which(!pruningStat)]



      if(length(certifiedParents) > 0){

        update <- updateDebiased_symb(z, certifiedParents, Y.orig, Y, Sigma, directEffect)

        totalEffect <- update$totalEffect
        directEffect <- update$directEffect
        Y[[z]] <- update$res.v

      } else {

        # If no parents, remove all from directEffect and update totalEffect
        directEffect[z, ] <- 0
        totalEffect <- solve(diag(rep(1, p)) -  directEffect)
        Y[[z]] <- Y.orig[[z]]
      }



    } else {
      if(verbose){
        cat(paste("Pruning:", z, "; no parents\n"))
      }
    }






  }

  return(list(totalEffect = totalEffect,
              directEffect = directEffect, dEdge = ifelse(directEffect != 0, 1, 0),
              bEdge = siblings + diag(rep(1, p))))
}

