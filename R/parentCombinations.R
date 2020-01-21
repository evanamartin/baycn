#' @importFrom utils combn
#'
parentCombinations <- function (pVector,
                                nNodes,
                                nGV,
                                nodeIndex,
                                pmr) {

  # total number of parents.
  totalParents <- sum(pVector)

  # Determine where the parents occur in the vector.
  whichParents <- which(pVector == 1)

  # Return a vector of zeros if the current node does not have any parents.
  if (totalParents == 0) {

    return (matrix(0,
                   nrow = nNodes,
                   ncol = 1))

  } else if (totalParents == 1) {

    # If the current node is a gv node remove any ge node parents.
    if (pmr & nodeIndex <= nGV & whichParents > nGV) {

      return (matrix(0,
                     nrow = nNodes,
                     ncol = 1))

    }

    allCombinations <- matrix(0,
                              nrow = nNodes,
                              ncol = 2)

    allCombinations[which(pVector == 1), 2] <- 1

    return (allCombinations)

  } else {

    # Get the total number of combinations for 2, 3, 4, ... parents.
    totalCombinations <- vector(mode = 'numeric',
                                length = totalParents)

    for (e in 1:totalParents) {

      totalCombinations[[e]] <- choose(totalParents, e)

    }

    # Create a matrix where each column is the parent vector for the current
    # node.
    allCombinations <- matrix(0,
                              nrow = nNodes,
                              ncol = 1 + sum(totalCombinations))

    # Next I need to determine all combinations of parent sets. For example,
    # 0 parents, 1 parent, 2 parents and so on
    # Create a counter that increases whenever v increases. Start the counter at
    # 1 because the first column represents no parents.
    counter <- 1
    for (e in 1:totalParents) {

      # Calculate all combinations of 1, 2, 3, ... parents
      combinations <- combn(whichParents, e)

      for (v in 1:totalCombinations[[e]]) {

        # increase the counter each time a new column vector is created.
        counter <- counter + 1

        # Change the parents from a 0 to a 1 in the column vector.
        allCombinations[combinations[, v], counter] <- 1

        # If the current node is a gv node remove any ge node parents.
        if (pmr & nodeIndex <= nGV) {

          allCombinations[(nGV + 1):nNodes, counter] <- 0

        }

      }

    }

    return (unique(allCombinations,
                   MARGIN = 2))

  }

}
