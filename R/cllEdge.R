#' cllEdge
#'
#' Calculates the log likelihood of the current node when running a genetic
#' algorithm on the direction of the edges in the graph.
#'
#' @param coordinates A matrix of the row and column numbers of the nonzero
#' elements in the adjacency matrix. The row numbers make up the first row of
#' the coordinates matrix and the column numbers make up the second row.
#'
#' @param data A matrix with the nodes across the columns and the observations
#' along the rows.
#'
#' @param individual The DNA of the current individual.
#'
#' @param scoreFun A character string indicating what method to use for
#' calculating the fitness.
#'
#' @return The log likelihood of the graph according to the orientation of the
#' edges denoted by the DNA of the current individual.
#'
#' @importFrom stats dnorm
#'
#' @export
#'
cllEdge <- function (individual,
                     coordinates,
                     data,
                     scoreFun) {

  # Convert the DNA of the current individual to an adjacency matrix
  adjMatrix <- toAdjMatrix(individual = individual,
                           coordinates = coordinates,
                           nNodes = ncol(data))

  # empty vector to hold the log likelihood for each node in the graph
  logLikelihood <- c()

  # Hold the number of parents per node. This will be used for calculating the
  # BIC if it is chosen as the fitness measure.
  k <- c()

  # Number of nodes in the graph. This will be used in for loops throughout the
  # function.
  nNodes <- ncol(adjMatrix)

  for (e in 1:nNodes) {

    # Get the number of parents, estimates for the mean and standard
    # deviation of the current node, and the data of the parent nodes
    ep <- estimates(adjMatrix = adjMatrix,
                    node = e,
                    data = data)

    k[[e]] <- ep$nParents

    if (ep$nParents == 0) {

      ll <- sum(dnorm(x = data[, e],
                      mean = ep$estimates[[1]],
                      sd = ep$estimates[[2]],
                      log = TRUE))

      logLikelihood[[e]] <- ll

    } else {

      ll <- sum(dnorm(x = data[, e],
                      mean = makeFormula(estimates = ep$estimates,
                                         nParents = ep$nParents,
                                         parentData = ep$parentData),
                      sd = ep$estimates[[length(ep$estimates)]],
                      log = TRUE))

      logLikelihood[[e]] <- ll

    }

  }

  switch(scoreFun,

         'logLikelihood' = {

           return (sum(logLikelihood))

         },

         'BIC' = {

           # Calculate the sample size. This will be used in calculating BIC.
           n <- dim(data)[2]

           # Vector to hold the BIC for each node in the graph.
           bic <- c()

           # Calculate BIC for each node in the graph.
           for (e in 1:nNodes) {

             bic[[e]] <- log(n) * (k[[e]] + 2) - 2 * logLikelihood[[e]]

           }

           return (sum(bic))

         })

}
