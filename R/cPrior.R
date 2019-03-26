#' cPrior
#'
#' Calculates the prior depending on the current edge state and edge type for
#' the mutate function. For example, gv-ge, gv-gv, ge-ge.
#'
#' @param edges A vector indicating the position in the prior vector of the two
#' edge states that the current edge can move to.
#'
#' @param edgeType A 0 or 1 indicating whether the edge is a gv-ge edge (1) or
#' a gv-gv or ge-ge edge (0).
#'
#' @param pmr Logical. If true the Metropolis-Hastings algorithm will use the
#' Principle of Mendelian Randomization, PMR. This prevents the direction of an
#' edge pointing from a gene expression node to a genetic variant node.
#'
#' @param prior A vector containing the prior probability of seeing each edge
#' direction.
#'
#' @return The probability of the specified edge state
#'
#' @export
#'
cPrior <- function (edges,
                    edgeType,
                    pmr,
                    prior) {

  if (pmr == TRUE) {

    # When using the pmr the probability of an edge, between a gv node and a ge
    # node, moving to edge state 1 is 0. The prior of edge state 2 is now the
    # prior of state 1 plus the prior of state 2.
    pmrPrior <- c(prior[[1]] + prior[[2]], 0, prior[[3]])

  }

  if (edgeType == 0) {

    # Calculate the probability of moving to the two possible edge states.
    priorProb <- c((prior[[edges[[1]]]] /
                      (prior[[edges[[1]]]] +
                         prior[[edges[[2]]]])),
                   (prior[[edges[[2]]]] /
                      (prior[[edges[[1]]]] +
                         prior[[edges[[2]]]])))

  } else {

    # Calculate the probability of moving to the two possible edge states when
    # using the pmr.
    priorProb <- c((pmrPrior[[edges[[1]]]] /
                      (pmrPrior[[edges[[1]]]] +
                         pmrPrior[[edges[[2]]]])),
                   (pmrPrior[[edges[[2]]]] /
                      (pmrPrior[[edges[[1]]]] +
                         pmrPrior[[edges[[2]]]])))

  }

  return (priorProb)

}
