#' caRatio
#'
#' Calculates the acceptance ratio for the Metropolis-Hastings and Genetic
#' Algorithm. It takes the ratio of the original vector of edge directions, old
#' individual, and the proposed vector of edge directions, new individual.
#'
#' @param m The length of the individual vector. This is the number of edges in
#' the graph and the fitness of the graph.
#'
#' @param proposed A vector of edge directions. This the the proposed
#' vector of edge directions.
#'
#' @param current A vector of edge directions. This is the graph from
#' the start of the current interation.
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
caRatio <- function (m,
                     proposed,
                     current,
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
    priorP[[e]] <- log(prior[[edgeDirP[[e]] + 1]])
    transProbP[[e]] <- log(prior[[edgeDirC[[e]] + 1]] /
                             sum(prior[-(edgeDirP[[e]] + 1)]))

    priorC[[e]] <- log(prior[[edgeDirC[[e]] + 1]])
    transProbC[[e]] <- log(prior[[edgeDirP[[e]] + 1]] /
                             sum(prior[-(edgeDirC[[e]] + 1)]))

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
