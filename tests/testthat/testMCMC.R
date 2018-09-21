context('Test MH algorithm')

test_that('the MHEdge function returns the correct matrix',{

  # Load the standards to test against.
  system.file('testdata', 'standards.RData', package = 'edgeFrequency')

  set.seed(338)

  # Generate data for topology M1 with no genetic variants.
  dataM1ge_200_1 <- m1ge(N = 200, ss = 1)

  # Adjacency matrix for toplogy M1
  adjMatrixM1 <- matrix(c(0, 1, 0,
                          0, 0, 1,
                          0, 0, 0),
                        byrow = TRUE,
                        nrow = 3)

  # Run the MH algorithm with the edges from the true graph.
  mhm1ge_200_1 <- MHEdge(adjMatrix = adjMatrixM1,
                         data = dataM1ge_200_1,
                         iterations = 500,
                         mutationRate = 1/2,
                         prior = c(0.05,
                                   0.05,
                                   0.9),
                         scoreFun = 'logLikelihood')

  # Generate data under topology H2
  dataH2_200_1 <- h2(N = 200, ss = 1)

  # Adjacency matrix for topology H2
  adjMatrixH2 <- matrix(c(0, 1, 1, 0, 0,
                          0, 0, 0, 1, 0,
                          0, 0, 0, 0, 1,
                          0, 0, 0, 0, 1,
                          0, 0, 0, 0, 0),
                        byrow = TRUE,
                        nrow = 5)

  # Run the MH algorithm with the edges from the true graph.
  mhh2_200_1 <- MHEdge(adjMatrix = adjMatrixH2,
                       data = dataH2_200_1,
                       iterations = 500,
                       mutationRate = 1/5,
                       prior = c(0.05,
                                 0.05,
                                 0.9),
                       scoreFun = 'logLikelihood')

  expect_identical(standardH2, mhh2_200_1)
  expect_identical(standardM1ge, mhm1ge_200_1)

})
