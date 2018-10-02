context('Principle of Medelian Randomization')

test_that('the toED function returns the correct vector of edge directions',{

  # The adjacency matrix below has the following structure
  # T1 <- V1 <- T2 -> V2 -> T3 -> T4
  # The two V nodes are the genetic variant nodes and the four T nodes are the
  # gene expression nodes. The two edges where the T nodes point to the V nodes
  # should be removed by the toED function.
  adjMatrix <- matrix(c(0, 0, 1, 0, 0, 0,
                        0, 0, 0, 0, 1, 0,
                        0, 0, 0, 0, 0, 0,
                        1, 1, 0, 0, 0, 0,
                        0, 0, 0, 0, 0, 1,
                        0, 0, 0, 0, 0, 0),
                      byrow = TRUE,
                      nrow = 6)

  # The coordinates of the adjacency matrix for each edge in the graph.
  coord <- coordinates(adjMatrix)

  # The edge directions of the graph according to the adjacency matrix.
  currentGraph <- c(0, 1, 1, 0, 0)

  newGraph <- toED(adjMatrix = adjMatrix,
                   coordinates = coord,
                   graph = currentGraph,
                   nGV = 2,
                   nNodes = 6)

  # The vector on the right is what we expect to see after running the toED
  # function. It should remove the edges between nodes V1 and T2 and between
  # V2 and T2.
  expect_identical(newGraph, c(0, 2, 2, 0, 0))

  # The adjacency matrix for topology M1 with edge one pointing from T1 to V.
  adjMatrix <- matrix(c(0, 0, 0,
                        1, 0, 1,
                        0, 0, 0),
                      byrow = TRUE,
                      nrow = 3)

  # Extract the coordinates of the edges from the adjacency matrix.
  coord <- coordinates(adjMatrix)

  # The edge directions of the graph according to the adjacency matrix.
  currentGraph <- c(1, 0)

  newGraph <- toED(adjMatrix = adjMatrix,
                   coordinates = coord,
                   graph = currentGraph,
                   nGV = 1,
                   nNodes = 3)

  # The vector on the right is what we expect to see after running the toED
  # function. It should remove the edge between nodes V and T1.
  expect_identical(newGraph, c(2, 0))

})
