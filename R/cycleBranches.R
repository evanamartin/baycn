#' cycleBranches
#'
#' Turns the stem and leaf matrix returned by the cycleSaL function into a list
#' of matrices. Each matrix contains all the cycles for each child of the root
#' node.
#'
#' @param SaL A binary matrix containing the nodes that could potentially form a
#' directed cycle.
#'
#' @return A matrix that contains all the cycles for the current graph. The
#' cycles are displayed down the columns of the matrix.
#'
#' @importFrom stats na.omit
#'
#' @export
#'
cycleBranches <- function (SaL) {

  # Start by getting the name of the root node.
  root <- rownames(SaL)[1]

  # Create a vector of the children of the root node.
  rChildren <- names(SaL[1, SaL[1, ] == 1])

  # Get the number of children of the root node.
  nChildren <- length(rChildren)

  # Create a list to hold the matrices of cycles. There will be one matrix for
  # each child of the root node.
  tree <- vector(mode = 'list',
                 length = nChildren)

  # This for loop will go through each child of the root node and build a matrix
  # that contains all of the nodes that could form a directed cycle down the
  # columns of the matrix.
  for (e in 1:nChildren) {

    # A matrix that holds the cycles of the eth child of the root node.
    tree[[e]] <- matrix(nrow = nrow(SaL) + 1,
                         ncol = 1)

    # The name in the first row/column will always be the name of the root node.
    tree[[e]][1, 1] <- root

    # The name in the second row first, column will always be the name of the
    # eth child of the root node.
    tree[[e]][2, 1] <- rChildren[[e]]

    # Create a variable for the number of columns of the matrices in tree. This
    # variable will be updated at the end of the following repeat loop to
    # account for branches further down the tree.
    nCol <- ncol(tree[[e]])

    # An index to move across the columns of the tree[[e]] matrix
    v <- 1

    # The following loop creates branches on the tree as deep as possible for
    # all possible branches for each child of the root node. The loop breaks
    # when no new branches need to be added to the tree[[e]] matrix and the
    # last branch is completed.
    repeat {

      if (v == 1) {

        # Get the potential children of the eth child of the root node. They are
        # potential children because we will not keep all of the children of the
        # current node. One of the child nodes will be the root node.
        pChildren <- names(SaL[tree[[e]][2, 1],
                                   SaL[tree[[e]][2, 1], ] == 1])

        # Eliminate the root node from the children of the current node. This is
        # done becuase it would lead to creating branches that bounced back and
        # forth between two nodes.
        children <- pChildren[pChildren != tree[[e]][1, 1]]

        # If there is more than one child create additional columns to hold each
        # child and their descendants.
        if (length(children) > 1) {

          # This for loop creates new rows in the tree matrix for each child.
          for (i in 1:length(children)) {

            # Copy the current column and add it to the end of the tree[[e]]
            # matrix length(children) - 1 times.
            if (i != 1) {

              # Create a new column in the tree[[e]] matrix for an additional
              # branch in the tree.
              tree[[e]] <- cbind(tree[[e]], tree[[e]][, 1])

              # Add the children to their respective columns
              tree[[e]][3, i] <- children[[i]]

              # This is run the first time through the for loop, i = 1, because
              # the child node needs to be added to the first column of the
              # matrix and not one of the columns added for the additional
              # children of the current node.
            } else {

              # Add the children to their respective columns.
              tree[[e]][3, i] <- children[[i]]

            }

          }

          # This is run if the current node only has one child.
        } else {

          # Add the child node to the row beneath the parent node.
          tree[[e]][3, 1] <- children

        }

        # nRow is a counter for moving down the rows of tree[[e]]. It starts at
        # 3 because the first three rows have already been added to the matrix.
        # It is also used in the while loop to check if the current node has
        # appeared previously.
        nRow <- 3

        # This runs when v is not equal to 1, when it is not the first time
        # through the repeat loop.
      } else {

        # Get the number of elements in the current column to start filling in
        # the descendants of the current node.
        nRow <- length(na.omit(tree[[e]][, v]))

      }

      # Loops through the current column, v, of tree[[e]] and adds nodes until
      # the loop comes to a node that has already been added.
      while (! tree[[e]][nRow, v] %in% tree[[e]][1:(nRow - 1), v]) {

        # Increase the row index by one.
        nRow <- nRow + 1

        # Potential children of the current parent.
        pChildren <- names(SaL[tree[[e]][nRow - 1, v],
                                   SaL[tree[[e]][nRow - 1, v], ] == 1])

        # Eliminate the child that is the same as the grandparent.
        children <- pChildren[pChildren != tree[[e]][nRow - 2, v]]

        # If there is more than one child create additional columns to hold each
        # child and their descendants.
        if (length(children) > 1) {

          # This for loop adds one column for each of the children after the
          # first child.
          for (i in 1:length(children)) {

            # Copy the current column and add it to the end of the tree[[e]]
            # matrix length(children) - 1 times.
            if (i != 1) {

              # Create a new column in the tree[[e]] matrix for the additional
              # branch in the tree.
              tree[[e]] <- cbind(tree[[e]], tree[[e]][, v])

              # Add the children to their respective columns
              tree[[e]][nRow, ncol(tree[[e]])] <- children[[i]]

              # This is run when i in the for loop that loops through the
              # children of the current parent is not equal to one.
            } else {

              # Add the children to their respective columns
              tree[[e]][nRow, v] <- children[[i]]

            }

          }

          # This runs when the current parent only has one child.
        } else {

          # Add the child node to the row beneath the parent node.
          tree[[e]][nRow, v] <- children

        }

      }

      # Get the number of columns currently in the matrix tree[[e]]
      nCol <- ncol(tree[[e]])

      # Exit the loop when the loop has filled in the last column of the
      # tree[[e]] matrix and no more columns need to be added.
      if (nCol == v) break

      # Increase the counter used to loop through the columns of the tree[[e]]
      # matrix.
      v <- v + 1

    }

  }

  # Combine all of the matrices into one matrix that contains all the branches.
  branches <- tree[[1]]

  for (b in 2:length(rChildren)) {

    branches <- cbind(branches, tree[[b]])

  }

  return (branches)

}
