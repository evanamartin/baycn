divergent <- function (coordinates,
                       currentES,
                       proposedES){

  dEdges <- which(currentES != proposedES)

  dNodes <- unique(as.vector(coordinates[, dEdges]))

  return (list(dEdges = dEdges,
               dNodes = dNodes))

}
