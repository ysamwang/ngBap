### Helper to scale off diagonal elements until the matrix is positive definite
makePD <- function(mat, backTrack = .97, tol = 1e-2){
  alpha <- 1
  ret <- mat * alpha + diag(rep(1- alpha, dim(mat)[2]))

  while(min(eigen(ret, only.values = T)$values) < tol){
    alpha <- alpha * backTrack
    ret <- mat * alpha + diag(rep(1 - alpha , dim(mat)[2]))
  }
  # cat(eigen(ret, only.values = )$val)
  return(ret)
}


### Helper to invert 1- B and get total effects
getTotal <- function(B){
  p <- dim(B)[1]
  return(solve(diag(rep(1, p)) - B))
}


#' Causal discovery with latent confounding and non-Gaussian data
#'
#'
#' Samples from a linear SEM corresponding to a BAP
#'
#' @param n number of observations
#' @param p number of variables
#' @param dist the distribution of the errors
#' \itemize{
#' \item  "gauss" for normal distribution
#' \item "t" for T distribution with t.df degrees of freedom
#' \item "gamma" for gamma distribution
#' \item  "unif" for uniform distribution
#' \item "laplace" for Laplace distribution
#' \item "lognormal" for lognormal distribution
#' }
#' @param d number of directed edges
#' @param b number of bidirected edges
#' @param B the direct edges; leave as null to randomly sample
#' @param Omega the covariance of epsilon; leave as null to randomly sample
#' @param directLow the lower bound for the absolute value of a randomly sampled directed edge weight
#' @param directHigh the upper bound for the absolute value of a randomly sampled directed edge weight
#' @param bidirectLow the lower bound for the absolute value of a randomly sampled bidirected edge weight
#' @param bidirectHigh the upper bound for the absolute value of a randomly sampled bidirected edge weight
#' @param t.df degree of freedom if distribution is T
#' @param ancestral T restricts the graph to be ancestral; F allows for a BAP which may not be ancestral
#' @param shuffle Setting shuffle = T randomly relabels the nodes; if set to F then 1... p will always is a valid causal ordering
#' @param signs Setting signs = F, all parameters are positive; setting signs = T allows both positive and negative values
#' @param nonAncestral Setting signs = F, all parameters are positive; setting signs = F allows both positive and negative values
#' @return
#' \itemize{
#' \item B the directed edge coefficients
#' \item Omega covariance of varepsilon
#' \item directEdges the adjacency matrix of the directed edges (support of B)
#' \item bidirectEdgesthe adjacency matrix of the bidirected edges (support of Omega)
#' \item errs the sampled epsilon values
#' \item Y the sampled Y values
#' \item Sigma the population covariance
#' \item ord a topological ordering of the nodes
#' }
#' @export
rBAP <- function(n, p, dist, d = NULL, b = NULL, B = NULL, Omega = NULL, directLow = .6, directHigh = 1,
                 bidirectLow = .3, bidirectHigh = .5, t.df = 13, ancestral = F, shuffle = T, signs = T, nonAncestral = F){


  if(is.null(B)){

    B <- matrix(0, p, p)
    B[sample(which(lower.tri(B)), size = d)] <- 1
    B.edge <- sample(c(-1, 1), size = p^2, replace = T) * runif(p^2, directLow, directHigh) * B

    if(!signs){

      B.edge <- abs(B.edge)

    }

  } else {

    B.edge <- B
    B <- ifelse(B != 0, 1, 0)

  }

  # Adjacency matrix for bidirected edges
  if(is.null(Omega)){

    Omega <- matrix(0, p, p)
    if(ancestral){

      ## Enforce ancestral condition by checking total effects
      validEdges <- which(getTotal(B) == 0 & lower.tri(B))

    } else {

      ## Enforce bow-free condition
      validEdges <- which(B == 0 & lower.tri(B))
    }

    if(length(validEdges) > 1){

      Omega[sample(validEdges, size = min(b, length(validEdges) )) ] <- 1

    } else if (length(validEdges) == 1 & b > 0){

      Omega[validEdges] <- 1

    }


    Omega.edge <- sample(c(-1, 1), size = p^2, replace = T) * runif(p^2, bidirectLow, bidirectHigh) * Omega
    Omega.edge <- Omega.edge + t(Omega.edge)
    diag(Omega.edge) <- 1
    # ensure that Omega.edge is PD with smallest eigenvalue of at least 1e-2


    ## Make Omega symmetric with 1's on diag
    Omega <- Omega + t(Omega) + diag(rep(1,p))

    if(!signs){

      Omega.edge <- abs(Omega.edge)

    }



  } else {

    Omega.edge <- Omega
    Omega <- ifelse(Omega != 0, 1, 0)

  }

  Omega.edge <- ngBap::makePD(Omega.edge, tol = 1e-2)


  # if signs is not true, then let all parameters be positive



  if(dist == "gauss"){
    errs <- mvtnorm::rmvnorm(n = n, mean = rep(0, p), sigma = Omega.edge)

  } else if(dist == "t") {
    errs <- mvtnorm::rmvt(n, sigma = Omega.edge * (t.df - 2)/ (t.df), df = t.df)

  } else if (dist == "gamma") {
    errs <- lcmix::rmvgamma(n, corr = Omega.edge) - 1

  } else if (dist == "unif") {
    errs <- (MultiRNG::draw.d.variate.uniform(n, p, Omega.edge) - .5) * 2 * sqrt(3)

  } else if (dist == "laplace") {
    errs <- LaplacesDemon::rmvl(n, mu = rep(0, p), Sigma = Omega.edge)

  } else if (dist == "lognormal"){
    errs <- mvtnorm::rmvnorm(n = n, mean = rep(0, p), sigma = Omega.edge)
    errs <- (exp(errs) - exp(1/2)) / sqrt(exp(1) * (exp(1) - 1))
    # errs <- (exp(errs) - exp(1/2))

  }

  Y <- solve(diag(rep(1, p)) - B.edge, t(errs))
  Y <- t(Y)
  sigma <- solve(diag(rep(1, p)) - B.edge, Omega.edge) %*% t(solve(diag(rep(1,p)) - B.edge))

  ord <- 1:p
  if (shuffle){
    ord <- sample(p)
    B.edge <- B.edge[ord, ord]
    Omega.edge <- Omega.edge[ord, ord]
    B <- B[ord, ord]
    Omega <- Omega[ord, ord]
    errs <- errs[, ord]
    Y <- Y[, ord]
    sigma <- sigma[ord, ord]
  }

  return(list(B = B.edge, Omega = Omega.edge, directEdges = B, bidirectEdges = Omega,
              errs = errs, Y = Y, Sigma = sigma, ord = ord))
}


