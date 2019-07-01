#' mhEdge
#'
#' The main function for the Metropolis-Hastings algorithm. It returns the
#' posterior distribution of the edge directions.
#'
#' @param adjMatrix The adjacency matrix returned by the MRPC algorithm. The
#' adjacency matrix is a matrix of zeros and ones. The ones represent an edge
#' and also indicates the direction of that edge.
#'
#' @param burnIn A number between 0 and 1 indicating the percentage of the
#' sample that will be discarded.
#'
#' @param data A matrix with the nodes across the columns and the observations
#' down the rows.
#'
#' @param iterations An integer for the number of iterations to run the MH
#' algorithm.
#'
#' @param nGV The number of genetic variants in the graph.
#'
#' @param pmr Logical. If true the Metropolis-Hastings algorithm will use the
#' Principle of Mendelian Randomization, PMR. This prevents the direction of an
#' edge pointing from a gene expression node to a genetic variant node.
#'
#' @param prior A vector containing the prior probability of seeing each edge
#' direction.
#'
#' @param thinTo An integer indicating the number of observations the sample
#' should be thinned to.
#'
#' @return A matrix where the first m columns are the edge directions and the
#' m + 1 column is the log likelihood. Each row represents one graph.
#'
#' @importFrom methods new
#'
#' @importFrom stats dbinom
#'
#' @examples
#' # Generate data under topology M1.
#' data_m1 <- m1gv(N = 200, p = 0.27, ss = 1)
#'
#' # Use the true edges as the inferred edges.
#' adjmatrix_m1 <- matrix(c(0, 1, 0,
#'                          0, 0, 1,
#'                          0, 0, 0),
#'                        byrow = TRUE,
#'                        nrow = 3)
#'
#' # Run the Metropolis-Hastings algorithm with the Principle of Mendelian
#' # Randomization (PMR).
#' mh_m1_pmr <- mhEdge(adjMatrix = adjmatrix_m1,
#'                     data = data_m1,
#'                     iterations = 500,
#'                     nGV = 1,
#'                     pmr = TRUE,
#'                     prior = c(0.05,
#'                               0.05,
#'                               0.9))
#'
#' # Run the Metropolis-Hastings algorithm on the same data as the previous
#' # example but with pmr set to FALSE.
#' mh_m1 <- mhEdge(adjMatrix = adjmatrix_m1,
#'                 data = data_m1,
#'                 iterations = 500,
#'                 nGV = 1,
#'                 pmr = FALSE,
#'                 prior = c(0.05,
#'                           0.05,
#'                           0.9))
#'
#' @export
#'
mhEdge <- function (adjMatrix,
                    burnIn = 0.2,
                    data,
                    iterations = 1000,
                    nGV,
                    pmr = FALSE,
                    prior = c(0.05,
                              0.05,
                              0.9),
                    thinTo = 200) {

  # Check for potential cycles in the graph and return the edge directions, edge
  # numbers, and the decimal numbers for each cycle if any exist.
  cp <- cyclePrep(adjMatrix)

  # Call the coordinates function to get the coordinates of all the nonzero
  # elements in the adjacency matrix and to create the DNA of the individuals.
  coord <- coordinates(adjMatrix = adjMatrix)

  # Get the number of edges in the graph.
  nEdges <- dim(coord)[2]

  # Get the edge types for each edge in the graph. The two edge types are gv-ge
  # or gv-gv, ge-ge.
  edgeType <- detType(coord = coord,
                      nEdges = nEdges,
                      nGV = nGV,
                      pmr = pmr)

  # Calculate the probability for 1:nEdges from a zero truncated binomial
  # distribution.
  ztbProb <- dbinom(x = 1:nEdges,
                    size = nEdges,
                    prob = 1 / nEdges) / (1 - (1 - 1/nEdges)^nEdges)

  # Get the number of nodes in the graph. This will be used in for loops
  # throughout the function.
  nNodes <- ncol(adjMatrix)

  # Create the length of the individual which is the number of edges plus the
  # fitness of the individual.
  m <- nEdges + 1

  # A vector of indices that indicate which iterations will be kept.
  if (burnIn == 0) {

    keep <- ceiling(seq(from = 1,
                        to = iterations,
                        length.out = thinTo))

  } else {

    keep <- ceiling(seq(from = burnIn * iterations,
                        to = iterations,
                        length.out = thinTo))

  }


  # A counter used to fill in the MarkovChain matrix.
  counter <- 0

  # Create a matrix to hold the accepted graph at each iteration of the MH
  # algorithm.
  MarkovChain <- matrix(nrow = thinTo,
                        ncol = m)

  # Create a vector to hold the decimal number for the accepted graph.
  graphDecimal <- vector(mode = 'integer',
                         length = thinTo)

  # Randomly generate a starting graph.
  graph <- sample(x = c(0, 1, 2),
                  size = nEdges,
                  replace = TRUE,
                  prob = prior)

  # Check if any gv nodes have ge node parents.
  if (pmr == TRUE) {

    # Change the direction of the ge -> gv edges to gv -> ge.
    graph[edgeType == 1 & graph == 1] <- 0

  }

  # Check the output of the cyclePrep function. If there are cycles in the graph
  # run the rmCycle function to remove any directed cycles.
  if (!is.null(cp$cycleDN)) {

    # Remove any directed cycles in the graph.
    graph <- rmCycle(cycleDN = cp$cycleDN,
                     edgeNum = cp$edgeNum,
                     edgeType = edgeType,
                     individual = graph,
                     pmr = pmr,
                     prior = prior)

  }

  # Calculate the log likelihood for the starting graph given the data.
  graph <- cllEdge(coordinates = coord,
                   data = data,
                   graph = graph,
                   m = m,
                   nEdges = nEdges,
                   nGV = nGV,
                   nNodes = nNodes)

  # Get the time that the for loop starts.
  startTime <- Sys.time()

  # Run through the iterations of the Metropolis Hastings algorithm.
  for (e in 1:iterations) {

    # Change the current graph to a different graph at each iteration.
    mGraph <- mutate(edgeType = edgeType,
                     graph = graph,
                     nEdges = nEdges,
                     pmr = pmr,
                     prior = prior,
                     ztbProb = ztbProb)

    # If there are potential directed cycles in the graph check if they are
    # present and remove them if they are.
    if (!is.null(cp$cycleDN)) {

      # Remove any directed cycles in the graph.
      mGraph <- rmCycle(cycleDN = cp$cycleDN,
                        edgeNum = cp$edgeNum,
                        edgeType = edgeType,
                        individual = mGraph,
                        pmr = pmr,
                        prior = prior)

    }

    # If the graph after mutation and removing directed cycles is the same as
    # the graph at the previous iteration then return that graph without
    # calculating the log likelihood.
    if (identical(mGraph[1:nEdges], graph[1:nEdges])) {

      # Store the accepted graph after burnin and thinning.
      if (e %in% keep) {

        # Update the counter if e is in the keep vector
        counter <- counter + 1

        # Store the accepted graph in the MarkovChain matrix
        MarkovChain[counter, ] <- graph

        # Calculate the decimal for the accepted graph.
        graphDecimal[[counter]] <- sum(graph[1:nEdges] * 3^(1:nEdges))

      }

      # Skip the rest of the for loop and go to the next iteration
      next

    }

    # Calculate the log likelihood of the changed graph
    mGraph <- cllEdge(coordinates = coord,
                      data = data,
                      graph = mGraph,
                      m = m,
                      nEdges = nEdges,
                      nGV = nGV,
                      nNodes = nNodes)

    # Determine whether to accept the proposed graph 'mGraph' or the original
    # graph 'graph'.
    graph <- caRatio(current = graph,
                     edgeType = edgeType,
                     m = m,
                     pmr = pmr,
                     proposed = mGraph,
                     prior = prior)

    # Store the accepted graph after burnin and thinning.
    if (e %in% keep) {

      # Update the counter if e is in the keep vector
      counter <- counter + 1

      # Store the accepted graph in the MarkovChain matrix
      MarkovChain[counter, ] <- graph

      # Calculate the decimal for the accepted graph.
      graphDecimal[[counter]] <- sum(graph[1:nEdges] * 3^(1:nEdges))

    }

  }

  # Get the time the for loop ends
  endTime <- Sys.time()

  # Return an object of class mcmc
  mcmcObj <- new('mcmc',
                 chain = MarkovChain[, 1:nEdges],
                 decimal = graphDecimal,
                 likelihood = MarkovChain[, m],
                 time = as.double((endTime - startTime),
                                  units = 'secs'))

  return (mcmcObj)

}
