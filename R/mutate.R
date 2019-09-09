#' @importFrom stats rbinom
#'
mutate <- function (edgeType,
                    graph,
                    nEdges,
                    pmr,
                    prior,
                    ztbProb) {

  # Generate a zero truncated binomial random variable to determine the number
  # of edges to change states.
  nChange <- sample(x = 1:nEdges,
                    size = 1,
                    prob = ztbProb)

  # Select which edges will be changed from their current state to a different
  # state.
  wEdges <- sample(x = 1:nEdges,
                   size = nChange,
                   replace = FALSE)

  # Loop through the edges to be mutated
  for (e in 1:nChange) {

    if (graph[[wEdges[[e]]]] == 0) {

      # Calculate the probability of moving to edge state 1 or 2 (location 2 or
      # 3 in the prior vector)
      probability <- cPrior(edges = c(2, 3),
                            edgeType = edgeType[[wEdges[[e]]]],
                            pmr = pmr,
                            prior = prior)

      graph[[wEdges[[e]]]] <- sample(x = c(1, 2),
                                     size = 1,
                                     prob = probability)

    } else if (graph[[wEdges[[e]]]] == 1) {

      # Calculate the probability of moving to edge state 0 or 2 (location 1 or
      # 3 in the prior vector)
      probability <- cPrior(edges = c(1, 3),
                            edgeType = edgeType[[wEdges[[e]]]],
                            pmr = pmr,
                            prior = prior)

      graph[[wEdges[[e]]]] <- sample(x = c(0, 2),
                                     size = 1,
                                     prob = probability)

    } else {

      # Calculate the probability of moving to edge state 0 or 1 (location 1 or
      # 2 in the prior vector)
      probability <- cPrior(edges = c(1, 2),
                            edgeType = edgeType[[wEdges[[e]]]],
                            pmr = pmr,
                            prior = prior)

      graph[[wEdges[[e]]]] <- sample(x = c(0, 1),
                                     size = 1,
                                     prob = probability)

    }

  }

  return (graph)

}
