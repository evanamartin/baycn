# cycleSaL
#
# The first step in identifying potential directed cycles in the graph is to
# find the nodes that could form a cycle. Since a node must be connected to at
# least two other nodes to be a part of a cycle we start by summing all of the
# rows of the undirected adjacency matrix and creating a new matrix that only
# has the rows of the undirected adjacency matrix with a sum of two or higher.
# The function cycleSaL creates the matrix of nodes that could form a directed
# cycle.
#
# Creates a stem and leaf matrix from the adjacency matrix with the nodes that
# could potentially form a cycle. The row names of the matrix created by the
# cycleSaL function are the stems and the 1s in the matrix are the leaves.
#
# @param adjMatrix The adjacency matrix of the graph.
#
# @param nGV The number of genetic variants in the graph.
#
# @param pmr Logical. If true the Metropolis-Hastings algorithm will use the
# Principle of Mendelian Randomization (PMR). This prevents the direction of an
# edge pointing from a gene expression node to a genetic variant node.
#
# @return  A binary matrix where the 1s represent an edge between the row node
# and the column node. If there are no cycles in the graph the function returns
# NULL.
#
cycleSaL <- function (adjMatrix,
                      nGV,
                      pmr) {

  # Add the adjacency matrix to its transpose to make all edges undirected.
  undirected <- adjMatrix + t(adjMatrix)

  # If a graph with undirected edges is passed to the removeCycles function then
  # the row sums could be 2 or greater even if that node is not in a potential
  # cycle because the cell in the undirected matrix for those nodes will be a 2.
  # To fix this problem we change any 2 in the undirected matrix to a 1.
  undirected[undirected == 2]  <- 1

  # Name the rows and columns of the undirected matrix from 1 to the number of
  # nodes in the adjacency matrix. These names will be used to remove the rows
  # of nodes that are not part of a cycle.
  rownames(undirected) <- c(1:dim(undirected)[2])
  colnames(undirected) <- c(1:dim(undirected)[2])

  # If pmr is TRUE remove all GV rows and columns from the SaL matrix otherwise
  # leave all columns from the undirected matrix in the SaL matrix.
  if (pmr == TRUE) {

    # First remove the genetic variant rows and columns.
    SaL <- undirected[-c(1:nGV), -c(1:nGV)]

    # Sum each row across all columns. The nodes that can potentially form a
    # cycle will have a sum of at least 2.
    sums <- apply(X = SaL,
                  MARGIN = 1,
                  FUN = sum)

    # Next remove rows that have a some of less than 2.
    SaL <- SaL[which(sums >= 2), ]

  } else {

    # Sum each row across all columns. The nodes that can potentially form a
    # cycle will have a sum of at least 2.
    sums <- apply(X = undirected,
                  MARGIN = 1,
                  FUN = sum)

    # Create a matrix whose rows are the nodes that could potentially form a
    # cycle.
    SaL <- undirected[which(sums >= 2), ]

  }

  # If there is only one row in the SaL matrix then it is turned into a vector
  # and we cannot give row names to a vector. The following if statement checks
  # if SaL is a vector or if it is a matrix with a dimension of two. If the
  # number of rows is 2 or fewer there cannot be a cycle in the graph.
  if (is.null(dim(SaL))) {

    return (NULL)

    # I had to break up the if statement because if SaL is a vector then the
    # statement in else if would return logical(0)
  } else if (dim(SaL)[1] <= 2) {

    return(NULL)

  }

  # Turn the row and column names of the SaL matrix into a numeric vector. These
  # will be used later to determine which entries in the cycles matrix need to
  # remain a one and which need to be changed to zero. The cells that will be
  # changed to a zero are the column numbers that do not appear in the rows.
  rNames <- as.numeric(rownames(SaL))
  cNames <- as.numeric(colnames(SaL))

  # Go through each element in the SaL matrix and turn the ones to zeros if
  # the column number the one occurs in does not appear in the rownames vector.
  for (e in seq_along(rNames)) {

    for (v in seq_along(cNames)) {

      if (SaL[e, v] == 1) {

        # If the column number doesn't show up in the row names change the one
        # to a zero in the cycles matrix.
        if (!any(rNames == cNames[[v]])) {

          SaL[e, v] <- 0

        }

      }

    }

  }

  # get the new row sums after the nodes that can't make up a cycle are removed
  # from the cycles matrix.
  newSums <- apply(X = SaL,
                   MARGIN = 1,
                   FUN = sum)

  # The following while lookp will continue to trim down the SaL matrix when
  # there is a string of nodes that are connected to at least two nodes but none
  # of the nodes form a cycle.
  while (any(newSums <= 1)) {

    # Only keep the rows with a sum >= 2
    SaL <- SaL[which(newSums >= 2), ]

    # Get the new row names
    rNames <- as.numeric(rownames(SaL))

    # If there is only one row left then the dim function will return NULL
    if (is.null(rNames)) {

      return (NULL)

      # If there are 2 or fewer rows with a row sum of 2 or greater there are no
      # cycles in the graph and cycleSaL returns a zero.
    } else if (length(rNames) <= 2) {

      return (NULL)

    } else {

      # Go through each element in the SaL matrix and turn the ones to zeros if
      # the column number the one occurs in does not appear in the rownames
      # vector.
      for (e in seq_along(rNames)) {

        for (v in seq_along(cNames)) {

          if (SaL[e, v] == 1) {

            # If the column number doesn't show up in the row names change the
            # one to a zero in the cycles matrix.
            if (!any(rNames == cNames[[v]])) {

              SaL[e, v] <- 0

            }

          }

        }

      }

      # This will be used to check if there are any rows that have a sum that is
      # less than 2. If there are then the while loop will continue. If not then
      # the while loop finishes.
      newSums <- apply(X = SaL,
                       MARGIN = 1,
                       FUN = sum)

    }

  }

  # Remove any columns that do not contain any 1s
  SaL <- SaL[, apply(X = SaL,
                     MARGIN = 2,
                     FUN = sum) != 0]

  return (SaL)

}

# cycleBranches
#
# Turns the stem and leaf matrix returned by the cycleSaL function into a list
# of matrices. Each matrix contains all the cycles for each child of the root
# node.
#
# @param SaL A binary matrix containing the nodes that could potentially form a
# directed cycle.
#
# @return A matrix that contains all the cycles for the current graph. The
# cycles are displayed down the columns of the matrix.
#
#' @importFrom stats na.omit
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

    # The name in the second row first column will always be the name of the
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

# cylceTB
#
# Takes the branches created by the cycleBranches function and trims them. The
# cycleBranches function creates a tree as deep as possible. Therefore, some of
# the branches may have nodes that are not part of a cycle.
#
# @param branches A matrix that holds all the branches of a tree.
#
# @return A list. The first element, tBranches, is a list of trimmed
# branches. The second element, nBranches, is an integer of the total number of
# branches in the tree.
#
# @importFrom stats na.omit
#
cycleTB <- function (branches) {

  # Extract the size of all cycles.
  cycleSize <- unique(apply(branches,
                           2,
                           function (x) sum(!is.na(x)) - 1))

  # Order the cycle sizes from smallest to largest.
  cycleSize <- cycleSize[order(cycleSize)]


  # Determine the number of different cycle sizes that are possible. For
  # example, 3 node cycles, four node cycles, and so on.
  ncs <- length(cycleSize)

  # Get the number of branches in the tree.
  nBranches <- dim(branches)[2]

  # A list to hold each trimmed branch according to the number of nodes in the
  # cycle. A trimmed branch contains only the nodes that appear in a cycle. The
  # branches passed to this function may contain nodes that are connected to the
  # cycle but are not part of the cycle. The tBranches list has length 3:nNodes.
  tBranches <- vector(mode = 'list',
                      length = ncs)

  # Create a counter for each list in tBranches
  counter <- rep(0, ncs)

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

    # Add the trimmed branch to the appropriate list. For example, a three node
    # cycle is added to the three node list and so on. The for loop loops
    # through the possible cycle sizes.
    for (cs in 1:ncs) {

      if ((length(branch) - 1) == cycleSize[[cs]]) {

        # update the counter
        counter[[cs]] <- counter[[cs]] + 1

        # Add the trimmed branch to the list corresponding to the number of
        # nodes in the trimmed branch.
        tBranches[[cs]][[counter[[cs]]]] <- branch

      }

    }

  }

  return (list(nBranches = nBranches,
               tBranches = tBranches))

}

# cycleCED
#
# Extracts the coordinates of the adjacency matrix for each edge that is part
# of a directed cycle and takes the names of the nodes in each cycle and
# creates a vector of edge directions that form a directed cycle.
#
# @param nBranches An integer. The number of branches in the tree.
#
# @param tBranches A list of cycles. Each element in the list is a vector of
# nodes that could form a directed cycle.
#
# @return A list with the coordinates of the adjacency matrix for each edge in
# each cycle, a vector of columns visited by the current cycles, and the edge
# directions for each cycle.
#
cycleCED <- function (nBranches,
                      tBranches) {

  # Extract the number of cycle sizes.
  ncs <- length(tBranches)

  # Create a list that will have a matrix for each cycle of size n. The rows of
  # the matrix will contain the adjacency matrix coordinates for each edge.
  cCoord <- vector(mode = 'list',
                   length = ncs)

  # A list that will hold the column indices of each column visited in the SaL
  # matrix by the current cycle.
  Columns <- vector(mode = 'list',
                    length = ncs)

  # The edgeDir list will hold the vectors of the edge directions that will form
  # a cycle in the graph.
  edgeDir <- vector(mode = 'list',
                    length = ncs)

  # Loop through each trimmed branch or cycle present in the
  # graph.
  for (e in 1:ncs) {

    # Extract the number of cycles for the current cycle size.
    nCycles <- length(tBranches[[e]])

    # Create a list for each cycle of size n for the adjacency matrix
    # coordinates.
    cCoord[[e]] <- vector(mode = 'list',
                          length = nCycles)

    # Create a list for each cycle of size n for the columns visited.
    Columns[[e]] <- vector(mode = 'list',
                           length = nCycles)

    # Create a list for each cycle of size n for the edge directions.
    edgeDir[[e]] <- vector(mode = 'list',
                           length = nCycles)

    # Loop through each cycle for the current cycle size (number of nodes in the
    # cycle).
    for (v in 1:nCycles) {

      # Get the nubmer of edges in the current cycle. The index e represents the
      # number of cycle sizes and v represents the number of cycles for the
      # current cycle size.
      nEdges <- length(tBranches[[e]][[v]]) - 1

      # A sublist that holds the column indices for each column visited in the
      # current cycle.
      Columns[[e]][[v]] <- vector(mode = 'integer',
                                  length = nEdges)

      # A vector to hold the directions of each edge for the directed cycle.
      edgeDir[[e]][[v]] <- numeric(nEdges)

      # Create a matrix to hold the adjacency matrix coordinates in the rows.
      cCoord[[e]][[v]] <- matrix(nrow = nEdges,
                                 ncol = 2)

      # This for loop gets the row and column indices for each edge in the
      # cycle.
      for (a in 1:nEdges) {

        # The row index is for the parent node and the column index is for the
        # child node. These indices will be used to determine if a 0 or 1 should
        # be present in the adjacency matrix in order for the nodes to form a
        # directed cycle.
        cCoord[[e]][[v]][a, ] <- c(tBranches[[e]][[v]][[a]],     # Row index
                                   tBranches[[e]][[v]][[a + 1]]) # Column Index

        # Keep the columns visited by each edge in the current cycle. This will
        # be used to determine if there are disjoint cycles in the graph.
        Columns[[e]][[v]][[a]] <- as.numeric(tBranches[[e]][[v]][[a + 1]])

        # If the value of the first node is less than the value of the second
        # node the direction of the edge points from the node with a smaller
        # index to the node with a larger index.
        if (as.numeric(tBranches[[e]][[v]][[a]]) <
            as.numeric(tBranches[[e]][[v]][[a + 1]])) {

          edgeDir[[e]][[v]][[a]] <- 0

          # If the value of the first node is greater than the value of the
          # second node the direction of the edge points from the node with a
          # larger index to the node with a smaller index.
        } else {

          edgeDir[[e]][[v]][[a]] <- 1

        }

      }

    }

  }

  # Get the unique column numbers for each column visited
  uColumns <- unique(unlist(Columns))

  return (list(cCoord = cCoord,
               columns = uColumns,
               edgeDir = edgeDir))

}

# cycleDJ
#
# Finds all disjoint cycles in the network.
#
# @param adjMatrix The adjacency matrix of the graph.
#
# @param SaL A binary matrix containing the nodes that could potentially form a
# directed cycle.
#
# @param ced A list containing the coordinates in the adjacency matrix of the
# nodes that could form a cycle, the columns of the SaL matrix that have been
# visited by the current cycles, and the edge directions for each cycle.
#
# @return A list with the coordinates of the adjacency matrix for each edge in
# each cycle, a vector of columns visited by the current cycles, and the edge
# directions for each cycle.
#
cycleDJ <- function (adjMatrix,
                     SaL,
                     ced) {

  # Keep the row names of the original SaL matrix. This will be used to check if
  # every row in the matrix has been used.
  rowNamesSaL <- as.numeric(row.names(SaL))

  # Column names visited so far
  columnNamesSaL <- ced$columns

  # Check if all the cycles have been found. If TRUE all cycles have been found.
  # If FALSE there are disjoint cycles and the cycleSaL function needs to be run
  # again.
  disjoint <- setequal(rowNamesSaL, columnNamesSaL)

  # If there are not disjoint cycles return NULL.
  if (disjoint == TRUE) {

    return (NULL)

  }

  # Loop through the previous funcions, eliminating rows from the SaL matrix at
  # each iteration of the while loop, until all cycles are found.
  while (!disjoint) {

    # Remove the column names that have been visited previously from the SaL
    # matrix.
    SaL2 <- SaL[-c(columnNamesSaL), ]

    # Create a tree as deep as possible with the remaining columns of the SaL
    # matrix.
    branches <- cycleBranches(SaL2)

    # Trim the branches starting at the leaf node to only include the nodes that
    # belong to the current cycle.
    trimmedB <- cycleTB(branches)

    # Extract the coordinates in the adjacency matrix for each edge in each
    # directed cycle and the edge states that create each directed cycle.
    ced2 <- cycleCED(nBranches = trimmedB$nBranches,
                     tBranches = trimmedB$tBranches)

    # Update the ced list.
    ced$cCoord <- append(ced$cCoord, ced2$cCoord)
    ced$columns <- append(ced$columns, ced2$columns)
    ced$edgeDir <- append(ced$edgeDir, ced2$edgeDir)

    # Update the columns that have been visited.
    columnNamesSaL <- append(columnNamesSaL, ced$columns)

    # Update the disjoint test.
    disjoint <- setequal(rowNamesSaL, columnNamesSaL)

  }

  # Return the ced list for the disjoint cycles.
  return (ced)

}

# combineCED
#
# Combines lists of the same cycle size for the cCoord and edgeDir lists.
#
# @param ced A list containing the coordinates in the adjacency matrix of the
# nodes that could form a cycle, the columns of the SaL matrix that have been
# visited by the current cycles, and the edge directions for each cycle.
#
# @return A list (the updated ced list) with the coordinates of the adjacency
# matrix for each edge in each cycle, a vector of columns visited by the current
# cycles, and the edge directions for each cycle.
#
combineCED <- function (ced) {

  # Create a counter to keep track of the number of lists in ced
  counter <- 1

  # Loop through each list in cCoord and edgeDir and combine the output of the
  # subsequent lists to the output of the previous lists if they have the same
  # cycle size.  For example, if the first list has output for cycles with three
  # nodes and the second list also hase output for cycles with three nodes
  # combine the first and second lists.
  repeat {

    # Compare the lists after counter to the list whose index match the counter.
    for (e in (counter + 1):length(ced$edgeDir)) {

      # If the cycle size is the same between the two lists combine them.
      if (length(ced$edgeDir[[counter]][[1]])==length(ced$edgeDir[[e]][[1]])) {

        # Append the eth cCoord list to the counter cCoord list.
        ced$cCoord[[counter]] <- append(ced$cCoord[[counter]],
                                        ced$cCoord[[e]])

        # Append the eth edgeDir list to the counter edgeDir list.
        ced$edgeDir[[counter]] <- append(ced$edgeDir[[counter]],
                                         ced$edgeDir[[e]])

        # Add e to a vector that will be used to remove lists from the cCoord
        # and edgeDir list.
        if (counter == 1) {

          rmList <- e

        } else {

          rmList <- append(rmList, e)

        }

      }

    }

    # Remove the eth list from cCoord becuase this was added to the list
    # with the counter index.
    ced$cCoord <- ced$cCoord[-rmList]

    # Remove the eth list from edgeDir becuase it was added to the list with
    # the counter index.
    ced$edgeDir <- ced$edgeDir[-rmList]

    # Update the counter
    counter <- counter + 1

    # Exit the loop if the remaining number of lists is the same as the counter.
    if (counter >= length(ced$edgeDir)) {

      break

    }

  }

  return (ced)

}

# cycleEID
#
# Extracts the edge index of each edge in the graph.
#
# @param adjMatrix The adjacency matrix of the graph.
#
# @param cCoord A list containing the coordinates in the adjacency matrix of
# the nodes that could form a cycle.
#
# @return A list containing the edge indices for the edges that could form a
# directed cycle for each cycle in the graph.
#
cycleEID <- function (adjMatrix,
                      cCoord) {

  # This line returns the coordinates of the adjacency matrix in a matrix where
  # the column numbers correspond to the position in the DNA of each edge. This
  # will be used to subset the individual to remove the possibility of a cycle
  # forming.
  position <- coordinates(adjMatrix)

  # Store the number of edges in the graph.
  m <- dim(position)[2]

  # Number of cycle sizes.
  ncs <- length(cCoord)

  # Stores the index of each edge for each cycle in the tBranches list
  edgeID <- vector(mode = 'list',
                   length = ncs)

  # Loop through each of the cycle sizes.
  for (e in 1:ncs) {

    # Number of directed cycles in the graph.
    nCycles <- length(cCoord[[e]])

    # Create a sublist of edgeID for each cycle in the current cycle size.
    edgeID[[e]] <- vector(mode = 'list',
                          length = nCycles)

    # Loop through the cycles in the current cycle size.
    for (v in 1:nCycles) {

      # Get the number of edges in the current cycle.
      nEdges <- dim(cCoord[[e]][[v]])[1]

      # Create a vector to store the edge number for each edge in the cycle.
      edgeID[[e]][[v]] <- vector(mode = 'numeric',
                                 length = nEdges)

      # Loop through each edge in the current cycle.
      for (a in 1:nEdges) {

        # This loops through all of the coordinates for each edge in the graph
        # and compares the coordinates of the current edge from the cycle to the
        # coordinates of each edge in the graph until it finds a match. When it
        # finds a match it stores the column number in the coordinates matrix
        # which is the edge number of the current edge in the graph.
        for (n in 1:m) {

          if (setequal(position[, n], cCoord[[e]][[v]][a, ])) {

            edgeID[[e]][[v]][[a]] <- n

            # Stop checking for the edge number once it is found.
            break

          }

        }

      }

    }

  }

  return (edgeID)

}

# cycleUNQ
#
# Only keeps the unique cycle decimals and IDs.
#
# @param cycleDN A list containing the decimal number of each cycle in each
# cycle size.
#
# @param edgeID A list containing the edge indices for the edges that could form
# a directed cycle for each cycle in the graph.
#
# @return A list containing the unique cycle decimals and IDs.
#
cycleUNQ <- function (cycleDN,
                      edgeID) {

  # Create a list of unique decimal numbers.
  uniqueDN <- vector(mode = 'list',
                     length = length(cycleDN))

  # Loop through each of the cycle sizes and only keep the unique decimals.
  for (e in 1:length(cycleDN)) {

    # Only keep the unique decimal numbers.
    uniqueDN[[e]] <- unique(cycleDN[[e]])

    # Keep the corresponding edgeID.
    edgeID[[e]] <- edgeID[[e]][match(uniqueDN[[e]], cycleDN[[e]])]

  }

  # Return the unique decimals and IDs.
  return (list(cycleDN = uniqueDN,
               edgeID = edgeID))

}

# cyclePrep
#
# Carries out all of the steps needed in order to remove a cycle.
#
# @param adjMatrix The adjacency matrix of the graph.
#
# @return A list containing the unique decimal numbers for the cycles in the
# graph and the indices of each edge in every cycle. If there are no cycles
# in the graph the function returns NULL.
#
cyclePrep <- function (adjMatrix,
                       nGV = 0,
                       pmr = FALSE) {

  # Check for potential cycles in the graph.
  SaL <- cycleSaL(adjMatrix = adjMatrix,
                  nGV = nGV,
                  pmr = pmr)

  # If there aren't any cycles return NULL.
  if (is.null(SaL)) {

    return (list(cycleDN = NULL,
                 edgeID = NULL))

  }

  # Creates a branch as deep as possible for each child of the parent node.
  branches <- cycleBranches(SaL)

  # Trim the branches starting at the leaf node to only include the nodes that
  # belong to the current cycle.
  trimmedB <- cycleTB(branches)

  # Extract the coordinates in the adjacency matrix for each edge in each
  # directed cycle and the edge states that create each directed cycle.
  ced <- cycleCED(nBranches = trimmedB$nBranches,
                  tBranches = trimmedB$tBranches)

  # If there are disjoint cycles find them.
  ced_temp <- cycleDJ(adjMatrix = adjMatrix,
                      SaL = SaL,
                      ced = ced)

  # check if ced_temp has any output. If it does then there were directed cycles
  # and ced will be replaced by ced_tmp.
  if (!is.null(ced_temp)) {

    # Update the ced list with ced_temp which contains the output for disjoint
    # cycles.
    ced <- combineCED(ced = ced_temp)

  }

  # Extract the indicies of the edges that form each cycle.
  edgeID <- cycleEID(adjMatrix,
                     ced$cCoord)

  # Convert each directed cycle into a decimal number.
  cycleDN <- decimal(ced$edgeDir,
                     edgeID)

  # Remove duplicate cycles
  dn_id <- cycleUNQ(cycleDN = cycleDN,
                    edgeID = edgeID)

  # Return the unique decimals and directed cycles
  return (dn_id)

}

# rmCycle
#
# Changes the edge direction of the current graph until there are no directed
# cylces.
#
# @param cycleDN A vector of decimal numbers for each directed cycle in the
# graph.
#
# @param edgeID A list containing the index for each of the edges in each of
# the potential cycles in the graph.
#
# @param edgeType A 0 or 1 indicating whether the edge is a gv-ge edge (1) or
# a gv-gv or ge-ge edge (1).
#
# @param individual A vector of the edge directions, and the log likelihood
# when the function is called within the populate function, of an individual.
#
# @param pmr Logical. If true the Metropolis-Hastings algorithm will use the
# Principle of Mendelian Randomization, PMR. This prevents the direction of an
# edge pointing from a gene expression node to a genetic variant node.
#
# @param prior A vector containing the prior probability of seeing each edge
# direction.
#
# @return A vector of the DNA of the individual with the cycles removed.
#
rmCycle <- function (cycleDN,
                     edgeID,
                     edgeType,
                     individual,
                     pmr,
                     prior) {

  # Get the decimal numbers for each potential cycle in the current individual.
  dnCI <- decimalCI(edgeID,
                    individual)

  # Get the number of cycle sizes.
  ncs <- length(cycleDN)

  # Create a list for the cycles in the current edge states vector (individual).
  whichCycle <- vector(mode = 'list',
                       length = ncs)

  # If there is a cycle change the direction of an edge invovled in the cycle.
  while (any(unlist(dnCI) == unlist(cycleDN))) {

    # Loop through each cycle size and check for cycles for that cycle size.
    for (e in 1:ncs) {

      # Get the edges that create a cycle from the current edge states and the
      # eth cycle size.
      whichCycle[[e]] <- which(dnCI[[e]] == cycleDN[[e]])

      if (length(whichCycle[[e]]) == 0) {

        # Exit the for loop
        break

      }

      # Loop through each cycle for the current cycle size.
      for (v in 1:length(whichCycle[[e]])) {

        # Get the edges involved in the eth cycle
        curEdges <- edgeID[[e]][[whichCycle[[e]][[v]]]]

        # Choose one of these edges and change it.
        whichEdge <- sample(x = curEdges,
                            size = 1)

        # If there is a directed cycle the edge state is either a 0 or 1.
        if(individual[[whichEdge]] == 0) {

          # Calculate the probability of moving to edge state 1 or 2 (location 2
          # or 3 in the prior vector)
          probability <- cPrior(edges = c(2, 3),
                                edgeType = edgeType[[whichEdge]],
                                pmr = pmr,
                                prior = prior)

          individual[[whichEdge]] <- sample(x = c(1, 2),
                                            size = 1,
                                            prob = probability)

        } else {

          # Calculate the probability of moving to edge state 0 or 2 (location 1
          # or 3 in the prior vector)
          probability <- cPrior(edges = c(1, 3),
                                edgeType = edgeType[[whichEdge]],
                                pmr = pmr,
                                prior = prior)

          individual[[whichEdge]] <- sample(x = c(0, 2),
                                            size = 1,
                                            prob = probability)

        }

      }

    }

    # Get the decimal numbers for each potential cycle in the current individual
    # after an edge in each cycle has been changed.
    dnCI <- decimalCI(edgeID,
                      individual)

  }

  return (individual)

}
