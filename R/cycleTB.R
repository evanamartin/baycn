#' cylceTB
#'
#' Takes the branches created by the cycleBranches function and trims them. The
#' cycleBranches function creates a tree as deep as possible. Therefore, some of
#' the branches may have nodes that are not part of a cycle.
#'
#' @param branches A matrix that holds all the branches of a tree.
#'
#' @return A list. The first element, tBranches, is a list of trimmed
#' branches. The second element, nBranches, is an integer of the total number of
#' branches in the tree.
#'
#' @importFrom stats na.omit
#'
cycleTB <- function (branches) {

  # Get the number of branches in the tree.
  nBranches <- dim(branches)[2]

  # A list to hold each trimmed branch. A trimmed branch contains only the nodes
  # that appear in a cycle. The branches passed to this function may contain
  # nodes that are connected to the cycle but are not part of the cycle.
  tBranches <- vector('list', nBranches)

  # The following loop will take each branch, starting at the leaf, and move
  # toward the root along the nodes until it sees a node that matches the leaf.
  # These nodes will form a cycle and they will be stored in tBranches.
  for (e in 1:nBranches) {

    # A vector of the nodes that form a cycle. They will be in reverse order
    # from how they were listed in the branches matrix.
    branch <- c()

    # Get the number of nodes in the eth column of the branches matrix excluding
    # any NAs.
    nNodes <- length(na.omit(branches[, e]))

    # The leaf or last node in the eth column of the branches matrix will be
    # used to stop the while loop.
    sNode <- branches[nNodes, e]

    # Add the leaf to the branch vector.
    branch[[1]] <- branches[nNodes, e]

    # Add the parent of the leaf to the branch vector.
    branch[[2]] <- branches[nNodes - 1, e]

    # The previous two lines are run in order to get into the while loop

    # create a counter to fill in the branch vector.
    v <- 2

    # Continue adding nodes the the branch vector until we come to back to the
    # node we started with.
    while (branch[[v]] != branch[[1]]) {

      v <- v + 1

      # Continue adding nodes to the vector
      branch[[v]] <- branches[nNodes - (v - 1), e]

    }

    # Add the trimmed branch to the list
    tBranches[[e]] <- branch

  }

  return (list(nBranches = nBranches,
               tBranches = tBranches))

}
