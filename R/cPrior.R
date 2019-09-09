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
