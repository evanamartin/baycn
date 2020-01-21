# lookup
#
# Creates an environment for each node in the graph and calculates the log
# likelihood for each node and for each combination of parent nodes, given the
# adjacency matrix.
#
# @param data A matrix with the variables across the columns and the
# observations down the rows. If there are genetic variants in the data these
# variables must come before the remaining variables. If there are clinical
# phenotypes in the data these variables must come after all other variables.
# For example, if there is a data set with one genetic variant variable, three
# gene expression variables, and one clinical phenotype variable the first
# column in the data matrix must contain the genetic variant data, the next
# three columns will contain the gene expression data, and the last column will
# contain the clinical phenotype data.
#
# @param adjMatrix An adjacency matrix indicating the edges that will be
# considered by the Metropolis-Hastings algorithm. This can be the output from
# another algorithm (e.g., PC). An adjacency matrix is a matrix of zeros and
# ones. The ones represent an edge and its direction between two nodes.
#
# @param nCPh The number of clinical phenotypes in the graph.
#
# @param nGV The number of genetic variants in the graph.
#
# @param nNodes The number of nodes in the graph.
#
# @param pmr Logical. If true the Metropolis-Hastings algorithm will use the
# Principle of Mendelian Randomization (PMR). This prevents the direction of an
# edge pointing from a gene expression or a clinical phenotype node to a
# genetic variant node.
#
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
