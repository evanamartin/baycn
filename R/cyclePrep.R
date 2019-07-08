#' cyclePrep
#'
#' Carries out all of the steps needed in order to remove a cycle.
#'
#' @param adjMatrix The adjacency matrix of the graph.
#'
#' @return A list containing the unique decimal numbers for the cycles in the
#' graph and the indices of each edge in every cycle. If there are no cycles
#' in the graph the function returns NULL.
#'
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
