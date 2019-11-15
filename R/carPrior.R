# carPrior
#
# Calculates the prior depending on the current edge state and edge type for
# the caRatio function. For example, gv-ge, gv-gv, ge-ge.
#
# @param edgeDir1 A scalar indicating the state of the edge for graph one.
#
# @param edgeDir2 A scalar indicating the state of the edge for graph two.
#
# @param edgeType A 0 or 1 indicating whether the edge is a gv-ge edge (1) or
# a gv-gv or ge-ge edge (1).
#
# @param pmr Logical. If true the Metropolis-Hastings algorithm will use the
# Principle of Mendelian Randomization, PMR. This prevents the direction of an
# edge pointing from a gene expression node to a genetic variant node.
#
# @param prior A vector containing the prior probability of seeing each edge
# direction.
#
# @return The probability of the specified edge state
#
carPrior <- function (edgeDir1,
                      edgeDir2,
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

    # Extract the prior for the edge state in graph one.
    priors <- log(prior[[edgeDir1 + 1]])

    # Calculate the probability of moving from the state in graph two to the
    # edge state in graph one.
    transition <- log(prior[[edgeDir2 + 1]] /
                        sum(prior[-(edgeDir1 + 1)]))

  } else {

    # Extract the prior for the edge state in graph one.
    priors <- log(pmrPrior[[edgeDir1 + 1]])

    # Calculate the probability of moving from the state in graph two to the
    # edge state in graph one.
    transition <- log(pmrPrior[[edgeDir2 + 1]] /
                        sum(pmrPrior[-(edgeDir1 + 1)]))

  }

  return (c(priors, transition))

}
