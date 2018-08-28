#' caRatio
#'
#' Calculates the acceptance ratio for the Metropolis-Hastings and Genetic
#' Algorithm. It takes the ratio of the original vector of edge directions, old
#' individual, and the proposed vector of edge directions, new individual.
#'
#' @param m The length of the individual vector. This is the number of edges in
#' the graph and the fitness of the graph.
#'
#' @param newIndividual A vector of edge directions. This the the proposed
#' vector of edge directions.
#'
#' @param oldIndividual A vector of edge directions. This is the individual from
#' the start of the current interation/generation
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
                     newIndividual,
                     oldIndividual,
                     prior,
                     scoreFun) {

  # To calculate the acceptance ratio we need to first have the prio probability
  # of the edges in the old individual and the new individual. Because they will
  # be the same for most of the edges we only need to consider the probabilities
  # of the edges that are different.
  difference <- newIndividual[1:(m - 1)] != oldIndividual[1:(m - 1)]

  # Ater getting the locations of the differences of the edge directions between
  # the old and new individuals we need to get the directions of the edges.
  diffNew <- newIndividual[difference]
  diffOld <- oldIndividual[difference]

  # After getting the edge directions for each edge that if different between
  # the old and new individuals we need to get the prior probability associated
  # with each edge direction.
  priorNew <- c()
  priorOld <- c()


  # Attach the correct prior probability to each edge direction for both the old
  # and new individuals.
  for(e in 1:length(diffNew)) {

    if (diffNew[[e]] == 0) {

      priorNew[[e]] <- log(prior[[1]])

    } else if (diffNew[[e]] == 1) {

      priorNew[[e]] <- log(prior[[2]])

    } else {

      priorNew[[e]] <- log(prior[[3]])

    }

    if (diffOld[[e]] == 0) {

      priorOld[[e]] <- log(prior[[1]])

    } else if (diffOld[[e]] == 1) {

      priorOld[[e]] <- log(prior[[2]])

    } else {

      priorOld[[e]] <- log(prior[[3]])

    }

  }

  switch (scoreFun,

          'logLikelihood' = {

            ratio <- ((sum(priorNew) + newIndividual[[m]])
                      - (sum(priorOld) + oldIndividual[[m]]))

            # Generate log uniform(0, 1) to compare to alpha which is
            # min(ratio, 0).
            logU <- log(runif(n = 1,
                              min = 0,
                              max = 1))

            alpha <- min(ratio, 0)

            # Determine if the proposed graph should be accepted.
            if (logU < alpha) {

              return (newIndividual)

            } else {

              return (oldIndividual)

            }

          },

          'BIC' = {

            # Because we want to minimize the BIC the ratio will change from
            # new/old to old/new. This will change the caRatio function from
            # always accepting a higher value to always accepting a lower value.
            ratio <- ((sum(priorOld) + oldIndividual[[m]])
                      - (sum(priorNew) + newIndividual[[m]]))

            # Generate log uniform(0, 1) to compare to alpha which is
            # min(ratio, 0).
            logU <- log(runif(n = 1,
                              min = 0,
                              max = 1))

            alpha <- min(ratio, 0)

            # Determine if the proposed graph should be accepted.
            if (logU < alpha) {

              return (newIndividual)

            } else {

              return (oldIndividual)

            }

          })

}
