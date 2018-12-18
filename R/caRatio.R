#' caRatio
#'
#' Calculates the acceptance ratio for the Metropolis-Hastings and Genetic
#' Algorithm. It takes the ratio of the original vector of edge directions, old
#' individual, and the proposed vector of edge directions, new individual.
#'
#' @param current A vector of edge directions. This is the graph from
#' the start of the current interation.
#'
#' @param edgeType A 0 or 1 indicating whether the edge is a gv-ge edge (1) or
#' a gv-gv or ge-ge edge (1).
#'
#' @param m The length of the individual vector. This is the number of edges in
#' the graph and the fitness of the graph.
#'
#' @param pmr Logical. If true the Metropolis-Hastings algorithm will use the
#' Principle of Mendelian Randomization, PMR. This prevents the direction of an
#' edge pointing from a gene expression node to a genetic variant node.
#'
#' @param proposed A vector of edge directions. This the the proposed
#' vector of edge directions.
#'
#' @param prior A vector containing the prior probability of seeing each edge
#' direction.
#'
#' @param scoreFun A character string indicating the fitness to be used. The
#' fitness can either be BIC of the log likelihood.
#'
#' @return A vector of edge directions and the fitness of the individual
#'
#' @importFrom stats runif
#'
#' @export
#'
caRatio <- function (current,
                     edgeType,
                     m,
                     pmr,
                     proposed,
                     prior,
                     scoreFun) {

  # To calculate the acceptance ratio we need to first have the prio probability
  # of the edges in the old individual and the new individual. Because they will
  # be the same for most of the edges we only need to consider the probabilities
  # of the edges that are different.
  difference <- proposed[1:(m - 1)] != current[1:(m - 1)]

  # Ater getting the locations of the differences of the edge directions between
  # the current and proposed graphs we need to get the directions of the edges.
  edgeDirP <- proposed[1:(m - 1)][difference]
  edgeDirC <- current[1:(m - 1)][difference]

  # Extract the edge type for the edges that have changed.
  etDiff <- edgeType[difference]

  # After getting the edge directions for each edge that is different between
  # the current and proposed graphs we need to get the prior probability
  # associated with each edge direction.
  priorP <- c()
  priorC <- c()

  # The following vectors will contain the transition probabilities for the
  # current and proposed graphs. transProbP will hold the probabilities for
  # moving from the proposed graph to the current graph and transProbC will hold
  # the probabilities of moving from the current graph to the proposed graph.
  transProbP <- c()
  transProbC <- c()

  # Attach the correct prior probability to each edge direction for both the old
  # and new individuals.
  for(e in 1:length(edgeDirP)) {

    # I can use the edge direction (0, 1, or 2) to select the correct prior by
    # adding a 1 to it and using that number to subset the prior vector. For
    # example, if the edge direction is 0 the corresponding prior is in the
    # first position of the prior vector so prior[[0 + 1]] will give 0.05 which
    # is the default prior for an edge being 0.

    # Calculate the log(prior) for the current edge state and the probability of
    # moving from the proposed graph to the current graph.
    ptProposed <- carPrior(edgeDir1 = edgeDirP[[e]],
                           edgeDir2 = edgeDirC[[e]],
                           edgeType = etDiff[[e]],
                           pmr = pmr,
                           prior = prior)

    priorP[[e]] <- ptProposed[1]
    transProbP[[e]] <- ptProposed[2]

    # Calculate the log(prior) for the current edge state and the probability of
    # moving from the current graph to the proposed graph.
    ptCurrent <- carPrior(edgeDir1 = edgeDirC[[e]],
                          edgeDir2 = edgeDirP[[e]],
                          edgeType = etDiff[[e]],
                          pmr = pmr,
                          prior = prior)

    priorC[[e]] <- ptCurrent[1]
    transProbC[[e]] <- ptCurrent[2]

  }

  switch (scoreFun,

          'logLikelihood' = {

            ratio <- ((sum(priorP) + proposed[[m]] + sum(transProbP))
                      - (sum(priorC) + current[[m]] + sum(transProbC)))

            # Generate log uniform(0, 1) to compare to alpha which is
            # min(ratio, 0).
            logU <- log(runif(n = 1,
                              min = 0,
                              max = 1))

            alpha <- min(ratio, 0)

            # Determine if the proposed graph should be accepted.
            if (logU < alpha) {

              return (proposed)

            } else {

              return (current)

            }

          },

          'BIC' = {

            # Because we want to minimize the BIC the ratio will change from
            # new/old to old/new. This will change the caRatio function from
            # always accepting a higher value to always accepting a lower value.
            ratio <- ((sum(priorC) + current[[m]])
                      - (sum(priorP) + proposed[[m]]))

            # Generate log uniform(0, 1) to compare to alpha which is
            # min(ratio, 0).
            logU <- log(runif(n = 1,
                              min = 0,
                              max = 1))

            alpha <- min(ratio, 0)

            # Determine if the proposed graph should be accepted.
            if (logU < alpha) {

              return (proposed)

            } else {

              return (current)

            }

          })

}
