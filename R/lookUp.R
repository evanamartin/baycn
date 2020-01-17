lookUp <- function (data,
                    adjMatrix,
                    nCPh,
                    nGV,
                    nNodes,
                    pmr) {

  # Add the adjacency matrix to its transpose to make all edges undirected.
  undirected <- adjMatrix + t(adjMatrix)

  # If a graph with undirected edges is passed to the removeCycles function then
  # the row sums could be 2 or greater even if that node is not in a potential
  # cycle because the cell in the undirected matrix for those nodes will be a 2.
  # To fix this problem we change any 2 in the undirected matrix to a 1.
  undirected[undirected == 2]  <- 1

  # Node Environments ----------------------------------------------------------

  # Create an environment that will hold all of the node environments. This
  # environment will be returned at the end of the function.
  allNodes <- new.env(parent = emptyenv())

  # Create a vector of names for each node
  nodeNames <- paste0('node', 1:nNodes)

  # Create an environment for each element in nodeNames.
  for (e in 1:nNodes) {

    assign(eval(nodeNames[[e]]),
           new.env(parent = emptyenv()),
           envir = allNodes)

  }

  # Parent Combinations/Log Likelihood -----------------------------------------

  # Create a list to hold the parent combinations matrix for each node in the
  # network.
  parents <- vector(mode = 'list',
                    length = nNodes)

  # Loop through each column in the adjacency matrix and create a matrix of all
  # possible parent combinations for each node.
  for (e in 1:nNodes) {

    parents[[e]] <- parentCombinations(pVector = undirected[, e],
                                       nNodes = nNodes,
                                       nGV = nGV,
                                       nodeIndex = e,
                                       pmr = pmr)

    # Calculate the log likelihood and store it in the environment for the
    # current node.
    for (v in 1:dim(parents[[e]])[2]) {

      assign(paste0('dn',
                    sum(parents[[e]][, v] * 2^(0:(nNodes - 1)))),
             cll(data = data,
                 nodeIndex = e,
                 v = v,
                 nCPh = nCPh,
                 nGV = nGV,
                 nNodes = nNodes,
                 pVector = parents[[e]][, v]),
             envir = eval(parse(text = paste0('allNodes$',
                                              nodeNames[[e]]))))

    }

  }

  return (allNodes)

}
