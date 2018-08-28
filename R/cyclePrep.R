#' cyclePrep
#'
#' Carries out all of the steps needed in order to remove a cycle.
#'
#' @param adjMatrix The adjacency matrix of the graph.
#'
#' @return A list containing the unique decimal numbers for the cycles in the
#' graph, the number of each edge in each cycle, and the location in the
#' original cycleDN list where each unique cycle occurs. If there are no cycles
#' in the graph the function returns NULL.
#'
#' @export
#'
cyclePrep <- function (adjMatrix) {

  # Check for potential cycles in the graph.
  SaL <- cycleSaL(adjMatrix)

  # If there aren't any cycles return NULL.
  if (is.null(SaL)) {

    return (NULL)

  }

  # If there are potential cycles in the graph then carry out the remaining
  # functions to remove them.
  branches <- cycleBranches(SaL)
  tBranches <- cycleTB(branches)
  cCoord <- cycleCoord(tBranches)
  edgeDir <- cycleED(tBranches)
  edgeNum <- cycleEN(adjMatrix,
                     cCoord)
  cycleDN <- decimal(edgeDir,
                     edgeNum)

  # Keep only the unique values of the decimal numbers.
  dnUnique <- unique(cycleDN)

  # Get the location in the cycleDN vector where each unique number occurs.
  # If there are multiple matches take only the first match.
  wCycle <- c()
  for (e in 1:length(dnUnique)) {

    wCycle[[e]] <- which(cycleDN == dnUnique[[e]])[[1]]

  }

  return (list(dnUnique = dnUnique,
               edgeNum = edgeNum,
               wCycle = wCycle))

}
