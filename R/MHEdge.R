#' mhEdge
#'
#' The main function for the Metropolis-Hastings algorithm. It returns the
#' posterior distribution of the edge directions.
#'
#' @param data A matrix with the nodes across the columns and the observations
#' down the rows.
#'
#' @param adjMatrix An adjacency matrix indicating the edges that will be
#' considered by the Metropolis-Hastings algorithm. This can be the output from
#' another algorithm (e.g., PC). An adjacency matrix is a matrix of zeros and
#' ones. The ones represent an edge and its direction between two nodes.
#'
#' @param burnIn A number between 0 and 1 indicating the percentage of the
#' sample that will be discarded.
#'
#' @param iterations An integer for the number of iterations to run the MH
#' algorithm.
#'
#' @param nGV The number of genetic variants in the graph.
#'
#' @param pmr Logical. If true the Metropolis-Hastings algorithm will use the
#' Principle of Mendelian Randomization (PMR). This prevents the direction of an
#' edge pointing from a gene expression node to a genetic variant node.
#'
#' @param prior A vector containing the prior probability of each edge state.
#'
#' @param thinTo An integer indicating the number of observations the chain
#' should be thinned to.
#'
#' @return An object of class mcmc containing 9 elements:
#'
#' \itemize{
#'
#' \item burnIn -- The percentage of MCMC iterations that will be discarded from
#' the beginning of the chain.
#'
#' \item chain -- A matrix where each row contains the vector of edge states for
#' the accepted graph.
#'
#' \item decimal -- A vector of decimal numbers. Each element in the vector is
#' the decimal of the accepted graph.
#'
#' \item iterations -- The number of iterations for which the
#' Metropolis-Hastings algorithm is run.
#'
#' \item posteriorES -- A matrix of posterior probabilities for all three edge
#' states for each edge in the network.
#'
#' \item posteriorPM -- A posterior probability adjacency matrix.
#'
#' \item likelihood -- A vector of log likelihood values. Each element in the
#' vector is the log likelihood of the accepted graph.
#'
#' \item stepSize -- The number of MCMC iterations discarded between each
#' accepted graph.
#'
#' \item time -- The runtime of the Metropolis-Hastings algorithm in seconds.
#'
#' }
#'
#' @references Evan A Martin and Audrey Qiuyan Fu. A Bayesian approach to
#' directed acyclic graph with a candidate graph.
#'
#' @importFrom methods new
#'
#' @importFrom stats dbinom
#'
#' @importFrom utils setTxtProgressBar txtProgressBar
#'
#' @examples
#' # Generate data under topology m1_gv.
#' data_m1 <- simdata(b0 = 0,
#'                    N = 200,
#'                    s = 1,
#'                    graph = 'm1_gv',
#'                    ss = 1,
#'                    p = 0.27)
#'
#' # Use the true edges for the input.
#' am_m1 <- matrix(c(0, 1, 0,
#'                   0, 0, 1,
#'                   0, 0, 0),
#'                 byrow = TRUE,
#'                 nrow = 3)
#'
#' # Run the Metropolis-Hastings algorithm with the Principle of Mendelian
#' # Randomization (PMR).
#' mh_m1_pmr <- mhEdge(data = data_m1,
#'                     adjMatrix = am_m1,
#'                     burnIn = 0.2,
#'                     iterations = 1000,
#'                     nGV = 1,
#'                     pmr = TRUE,
#'                     prior = c(0.05,
#'                               0.05,
#'                               0.9),
#'                     thinTo = 200)
#'
#' summary(mh_m1_pmr)
#'
#' # Generate data under topology gn4.
#' data_gn4 <- simdata(b0 = 0,
#'                     N = 200,
#'                     s = 1,
#'                     graph = 'gn4',
#'                     ss = 1)
#'
#' # Use the true edges for the input.
#' am_gn4 <- matrix(c(0, 1, 1, 0,
#'                    0, 0, 0, 1,
#'                    0, 0, 0, 0,
#'                    0, 0, 1, 0),
#'                  byrow = TRUE,
#'                  nrow = 4)
#'
#' mh_gn4 <- mhEdge(data = data_gn4,
#'                  adjMatrix = am_gn4,
#'                  burnIn = 0.2,
#'                  iterations = 1000,
#'                  nGV = 0,
#'                  pmr = FALSE,
#'                  prior = c(0.05,
#'                            0.05,
#'                            0.9),
#'                  thinTo = 200)
#'
#' summary(mh_gn4)
#'
#' @export
#'
mhEdge <- function (data,
                    adjMatrix,
                    burnIn = 0.2,
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
  # elements in the adjacency matrix.
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

  #################################################
  # Beginning of the Metropolis-Hastings algorithm.
  #################################################

  # create progress bar
  pb <- txtProgressBar(min = 0, max = iterations, style = 3)

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

    # Update progress bar
    setTxtProgressBar(pb, e)

  }

  # Close the progress bar
  close(pb)

  #################################################
  # End of the Metropolis-Hastings algorithm.
  #################################################

  # Get the time the for loop ends
  endTime <- Sys.time()

  #################################################
  # Create the probability matrix for the three
  # edge states.
  #################################################

  # Create a matrix to hold the posterior probability of each edge.
  posteriorES <- as.data.frame(matrix(nrow = nEdges,
                                      ncol = 4))

  # zero - the edge points from the node with a smaller index to the
  # node with a larger index
  # one - the edge points from the node with a larger index to the
  # node with a smaller index
  # two - the edge is absent between the two nodes
  colnames(posteriorES) <- c('nodes', 'zero', 'one', 'two')

  rownames(posteriorES) <- paste('edge', 1:nEdges, sep = '')

  # Calculate the proportion of 0, 1, and 2 for each edge in the
  # network.
  for (e in 1:nEdges) {

    # calculate the posterior probability for state 0, 1, and 2
    edge_0 <- sum(MarkovChain[, e] == 0)
    edge_1 <- sum(MarkovChain[, e] == 1)
    edge_2 <- sum(MarkovChain[, e] == 2)

    posteriorES[e, 2:4] <- c(round(edge_0 / thinTo, 4),
                             round(edge_1 / thinTo, 4),
                             round(edge_2 / thinTo, 4))

    # Determine which nodes the edge is between
    posteriorES[e, 1] <- paste0(colnames(data)[coord[1, e]],
                                '-',
                                colnames(data)[coord[2, e]])

  }

  #################################################
  # Create the probabilistic adjacency matrix.
  #################################################

  # Probabilities for upper triangular matrix
  upper <- posteriorES[, 2]

  # Probabilities for lower triangular matrix
  lower <- posteriorES[, 3]

  # Initialize a pxp matrix to fill in later with the estimated probabilities
  # from bacyn.
  posteriorPM <- matrix(0, nrow = nNodes, ncol = nNodes)

  # loop through the row column coordinates of the true edge indicies to fill in
  # the posterior edge probabilities.
  for (v in 1:nEdges) {

    # Fill in the upper triangular matrix from the '0' column.
    posteriorPM[coord[1, v], coord[2, v]] <- upper[[v]]

    # Fill in hte lower triangular matrix from the '1' column.
    posteriorPM[coord[2, v], coord[1, v]] <- lower[[v]]

  }

  # name the rows and columns with the names from the data matrix.
  colnames(posteriorPM) <- colnames(data)
  rownames(posteriorPM) <- colnames(data)

  # Calculate the stepSize based on the burnIn, iterations, and thinTo
  stepSize <- ifelse(burnIn == 0,
                     ceiling(iterations / thinTo),
                     ceiling(((1 - burnIn) * iterations) / thinTo))

  # Return an object of class mcmc
  mcmcObj <- new('mcmc',
                 burnIn = burnIn * 100,
                 chain = MarkovChain[, 1:nEdges],
                 decimal = graphDecimal,
                 iterations = iterations,
                 posteriorES = posteriorES,
                 posteriorPM = posteriorPM,
                 likelihood = MarkovChain[, m],
                 stepSize = stepSize,
                 time = as.double((endTime - startTime),
                                  units = 'secs'))

  return (mcmcObj)

}
