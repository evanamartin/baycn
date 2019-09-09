cycleSaL <- function (adjMatrix) {

  # Add the adjacency matrix to its transpose to make all edges undirected.
  undirected <- adjMatrix + t(adjMatrix)

  # If a graph with undirected edges is passed to the removeCycles function then
  # the row sums could be 2 or greater even if that node is not in a potential
  # cycle because the cell in the undirected matrix for those nodes will be a 2.
  # To fix this problem we change any 2 in the undirected matrix to a 1.
  undirected[undirected == 2]  <- 1

  # Sum each row across all columns. The nodes that can potentially form a cycle
  # will have a sum greater than 2.
  sums <- apply(X = undirected,
                MARGIN = 1,
                FUN = sum)

  # Create a matrix whose rows are the nodes that could potentially form a
  # cycle.
  SaL <- undirected[which(sums >= 2), ]

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

  # Name the rows according to their placement in the undirected matrix. The
  # name of each row is the row number of the undirected adjacency matrix where
  # there was a sum of 2 or greater.
  rownames(SaL) <- which(sums >= 2)

  # The column names are also the column numbers of the undirected adjacency
  # matrix. At this point in the function we keep all of the columns in the
  # undirected adjacency matrix.
  colnames(SaL) <- c(1:dim(undirected)[2])

  # Get the number of rows and columns in the matrix of potential cycles. These
  # will be used later to determine which entries in the cycles matrix need to
  # remain a one and which need to be changed to zero. The cells that will be
  # changed to a zero are the column numbers that don't also appear in the rows.
  nRows <- dim(SaL)[1]
  nCols <- dim(SaL)[2]

  # Go through each element in the cycles matrix and turn the ones to zeros if
  # the column number the one occurs in does not appear in the rownames vector.
  for (e in 1:nRows) {

    for (v in 1:nCols) {

      if (SaL[e, v] == 1) {

        # If the column number doesn't show up in the row names change the one
        # to a zero in the cycles matrix.
        if (!any(as.numeric(rownames(SaL)) == v)) {

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

    # Get the new row dimensions
    nRows <- dim(SaL)[1]

    # If there is only one row left then the dim function will return NULL
    if (is.null(nRows)) {

      return (NULL)

      # If there are 2 or fewer rows with a row sum of 2 or greater there are no
      # cycles in the graph and cycleSaL returns a zero.
    } else if (nRows <= 2) {

      return (NULL)

    } else {

      # Go through each element in the cycles matrix and turn the ones to zeros
      # if the column number the one occurs in does not appear in the rownames
      # vector.
      for (e in 1:nRows) {

        for (v in 1:nCols) {

          if (SaL[e, v] == 1) {

            # If the column number doesn't show up in the row names change the
            # one to a zero in the cycles matrix.
            if (!any(as.numeric(rownames(SaL)) == v)) {

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

cycleCED <- function (nBranches,
                      tBranches) {

  # Create a list that will have a list for each edge in each directed cycle.
  # The first level of lists are for each cycle and the sublists are for each
  # edge in the cycle.
  cCoord <- vector(mode = 'list',
                   length = nBranches)

  # A list that will hold the column indices of each column visited in the SaL
  # matrix by the current cycle.
  Columns <- vector(mode = 'list',
                    length = nBranches)

  # The edgeDir list will hold the vectors of the edge directions that will form
  # a cycle in the graph.
  edgeDir <- vector(mode = 'list',
                    length = nBranches)

  # Loop through each trimmed branch or cycle present in the
  # graph.
  for (e in 1:nBranches) {

    # Get the nubmer of edges in the current cycle
    nEdges <- length(tBranches[[e]]) - 1

    # A sublist that holds the columns indices for each column visited in the
    # current cycle.
    Columns[[e]] <- vector(mode = 'integer',
                           length = nEdges)

    # A vector to hold the directions of each edge for the directed cycle.
    edgeDir[[e]] <- numeric(nEdges)

    # This for loop gets the row and column indices for each edge in the cycle.
    for (v in 1:nEdges) {

      # The row index is for the parent node and the column index is for the
      # child node. These indices will be used to determine if a 0 or 1 should
      # be present in the adjacency matrix in order for the nodes to form a
      # directed cycle.
      #                        Row index            Column index
      cCoord[[e]][[v]] <- c(tBranches[[e]][[v]], tBranches[[e]][[v + 1]])

      # Keep the columns visited by each edge in the current cycle. This will be
      # used to determine if there are disjoint cycles in the graph.
      Columns[[e]][[v]] <- as.numeric(tBranches[[e]][[v + 1]])

      # If the value of the first node is less than the value of the second node
      # the direction of the edge points from the node with a smaller index to
      # the node with a larger index.
      if (as.numeric(tBranches[[e]][[v]]) <
          as.numeric(tBranches[[e]][[v + 1]])) {

        edgeDir[[e]][[v]] <- 0

        # If the value of the first node is greater than the value of the second
        # node the direction of the edge points from the node with a larger
        # index to the node with a smaller index.
      } else {

        edgeDir[[e]][[v]] <- 1

      }

    }

  }

  # Get the unique column numbers for each column visited
  uColumns <- unique(unlist(Columns))

  return (list(cCoord = cCoord,
               columns = uColumns,
               edgeDir = edgeDir))

}

cycleEN <- function (adjMatrix,
                     cCoord) {

  # This line returns the coordinates of the adjacency matrix in a matrix where
  # the column numbers correspond to the position in the DNA of each edge. This
  # will be used to subset the individual to remove the possibility of a cycle
  # forming.
  position <- coordinates(adjMatrix)

  # Number of directed cycles in the graph.
  nCycles <- length(cCoord)

  # Stores the location of each edge for each cycle in the tBranches list
  edgeNum <- vector(mode = 'list',
                    length = nCycles)

  # The first for loop loops through each of the cycles present in the graph.
  for (e in 1:nCycles) {

    nEdges <- length(cCoord[[e]])

    # Create a list to store the edge number for each edge in the cycle.
    edgeNum[[e]] <- vector(mode = 'numeric',
                           length = nEdges)

    # This for loop is the length of the current cycle, the number of edges in
    # the cycle.
    for (v in 1:nEdges) {

      # This loops through all of the coordinates for each edge in the graph and
      # compares the coordinates of the current edge from the cycle to the
      # coordinates of each edge in the graph until it finds a match. When it
      # finds a match it stores the column number in the coordinates matrix
      # which is the edge number of the current edge in the graph.
      for (a in 1:dim(position)[2]) {

        if (setequal(position[, a], cCoord[[e]][[v]])) {

          edgeNum[[e]][[v]] <- a

          # Stop checking for the edge number once it is found.
          break

        }

      }

    }

  }

  return (edgeNum)

}

cyclePrep <- function (adjMatrix) {

  # Check for potential cycles in the graph.
  SaL <- cycleSaL(adjMatrix)

  # If there aren't any cycles return NULL.
  if (is.null(SaL)) {

    return (list(cycleDN = NULL,
                 edgeNum = NULL))

  }

  # Creates a branch as deep as possible for each child of the parent node.
  branches <- cycleBranches(SaL)

  # Trim the branches starting at the leaf node to only the nodes that belong to
  # the cycle.
  tBranches <- cycleTB(branches)

  # Extract the coordinates in the adjacency matrix for each edge in each
  # directed cycle and the edge states that create a directed cycle.
  ced <- cycleCED(nBranches = tBranches$nBranches,
                  tBranches = tBranches$tBranches)

  # Extract the indicies of the edges that form a cycle.
  edgeNum <- cycleEN(adjMatrix,
                     ced$cCoord)

  # Convert each directed cycle into a decimal number.
  cycleDN <- decimal(ced$edgeDir,
                     edgeNum)

  # Keep the row names of the original SaL matrix. This will be used to check if
  # every row in the matrix has been used.
  rowNamesSaL <- as.numeric(row.names(SaL))

  # Column names visited so far
  columnNamesSaL <- ced$columns

  # Check if all the disjoint cycles have been found (if any).
  disjoint <- setequal(rowNamesSaL, columnNamesSaL)

  # Loop through the previous funcions, eliminating rows from the SaL matrix at
  # each iteration of the while loop, until all cycles are found.
  while (!disjoint) {

    SaL2 <- SaL[-c(columnNamesSaL), ]
    branches <- cycleBranches(SaL2)
    tBranches <- cycleTB(branches)
    ced <- cycleCED(nBranches = tBranches$nBranches,
                    tBranches = tBranches$tBranches)
    edgeNum2 <- cycleEN(adjMatrix,
                        ced$cCoord)
    cycleDN2 <- decimal(ced$edgeDir,
                        edgeNum2)

    # Add the new cycles, columns, and decimals to edgeNum and cycleDN
    edgeNum <- append(edgeNum, edgeNum2)
    columnNamesSaL <- append(columnNamesSaL, ced$columns)
    cycleDN <- append(cycleDN, cycleDN2)

    # Update the disjoint test
    disjoint <- setequal(rowNamesSaL, columnNamesSaL)

  }

  # Keep only the unique values of the decimal numbers.
  uniqueDN <- unique(cycleDN)

  # Return the unique decimals and directed cycles
  return (list(cycleDN = uniqueDN,
               edgeNum = edgeNum[match(uniqueDN, cycleDN)]))

}

rmCycle <- function (cycleDN,
                     edgeNum,
                     edgeType,
                     individual,
                     pmr,
                     prior) {

  # Get the decimal numbers for each potential cycle in the current individual.
  dnCI <- decimalCI(edgeNum,
                    individual)

  # If there is a cycle change the direction of an edge invovled in the cycle.
  while (any(dnCI == cycleDN)) {

    # Get the edges that create a cycle in the current individual.
    whichCycle <- which(dnCI == cycleDN)

    for (e in 1:length(whichCycle)) {

      # Get the edges involved in the eth cycle
      curEdges <- edgeNum[[whichCycle[[e]]]]

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

      # Check for cycles after each iteration of the for loop. If there are no
      # cycles break out of the for loop.
      dnCI <- decimalCI(edgeNum,
                        individual)

      if (!any(dnCI == cycleDN)) {

        break

      }

    }

    # Get the decimal numbers for each potential cycle in the current individual
    # after an edge in each cycle has been changed.
    dnCI <- decimalCI(edgeNum,
                      individual)

  }

  return (individual)

}
