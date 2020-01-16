# lull
#
# A function to look up the likelihood for a given set of nodes based on the
# current parent configuration for each node.
#
# am - the adjacency matrix of the current graph
#
# likelihood - A vector containing the log likelihood for each node.
#
# llenv - the environment that has the likelihood for each parent combintation
# for each node.
#
# nNodes - The number of nodes in the graph.
#
# wNodes - the nodes for which the likelihood will be looked up.
#
lull <- function (am,
                  likelihood,
                  llenv,
                  nNodes,
                  wNodes) {

  # Check if wNodes has length zero. If it doesn't then the proposed and current
  # edge state vectors are different.
  if (length(wNodes) != 0) {

    # Loop through each node that has changed
    for (e in wNodes) {

      # Create the name for the current node. This will be used to find the
      # environment containing the likelihood values for the current node.
      node <- paste0('node', e)

      # Create the name for the current parent configuration. This will be used
      # to look up the likelihood based on the current parent configuration.
      parents <- paste0('dn',
                        sum(am[, e] * 2^(0:(nNodes - 1))))

      likelihood[[e]] <- llenv[[node]][[parents]]

    }

  }

  return (likelihood)

}
