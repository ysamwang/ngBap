#'  Causal discovery with latent confounding and non-Gaussian data
#'
#' Given a mixed graph, projects it into a PAG using fciPlus
#'
#'
#' @param B matrix of directed edges
#' @param O matrix of bidirected edges
#' @return the mag output of fciPlus given generic edgeweights
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
#'
mag2pag <- function(B, O){
  p <- dim(B)[1]
  B.weights <- runif(p^2, .05, .95) * B
  O.weights <- matrix(0, p, p)
  O.weights[lower.tri(O) & O == 1] <- runif(sum(lower.tri(O) & O == 1), 1 / (2 * p), 1 / (p))
  O.weights <- O.weights + t(O.weights) + diag(rep(1, p))
  sigma <- solve(diag(rep(1, p)) - B.weights) %*% O.weights %*% t(solve(diag(rep(1, p)) - B.weights))
  mag <- pcalg::fciPlus(suffStat = list(C = cov2cor(sigma), n = 10^9), alpha = .999, p = p, indepTest = pcalg::gaussCItest)
  return(mag)
}
