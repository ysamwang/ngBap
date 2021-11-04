#'  Causal discovery with latent confounding and non-Gaussian data
#'
#' Given a matrix of estimated total effects, returns a causal ordering
#'
#'
#' @param totalEffect a p x p matrix of estimated total effects
#' @return a topological ordering which respects the total effects so that i < j if i is not a descendent of j
#'}
#' @examples
#' \dontrun{
#' dat <- ngBap::rBAP(n = 50000, p = 7, dist = "gamma", d = 3, b = 5, ancestral = F, shuffle = T, signs = T)
#' Y <- dat$Y
#' out <- ngBap::bang(Y, K = 3, level = .01, verbose = F, restrict = 1)
#' }
#'
#'
#' @export
#'
getOrdering <- function(totalEffect){
  p <- dim(totalEffect)[1]
  diag(totalEffect) <- 0

  ancestors <- lapply(split(totalEffect, rep(1:p, each = p)), function(x){which(x != 0)})


  ordering <- rep(0, p)
  ordering[1] <- 1
  for(k in 2:p){

    j <- k - 1
    while(j > 0 && !(ordering[j] %in% ancestors[[k]])){
      j <- j - 1
    }

    ## Insert k at index j ##
    if(j < (k - 1)){
      temp <- ordering[(j+1):(k-1)]
      ordering[(j+2):k] <- temp
    }

    ordering[j+1] <- k

  }
  return(ordering)
}

