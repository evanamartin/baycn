context('MCMC no genetic variants')

test_that('the MHEdge function returns the correct matrix',{

  # Load the standards to test against.
  load(system.file('testdata',
                   'standards.RData',
                   package = 'edgeFrequency'))

  set.seed(338)

  # Generate data for topology M1 with no genetic variants.
  m1ge_200_1_data <- m1ge(N = 200, ss = 1)

  # Adjacency matrix for toplogy M1
  adjMatrixM1 <- matrix(c(0, 1, 0,
                          0, 0, 1,
                          0, 0, 0),
                        byrow = TRUE,
                        nrow = 3)

  # Run the MH algorithm with the edges from the true graph.
  m1ge_200_1_mh <- MHEdge(adjMatrix = adjMatrixM1,
                          data = m1ge_200_1_data,
                          iterations = 500,
                          mutationRate = 1/2,
                          nGV = 0,
                          pmr = FALSE,
                          prior = c(0.05,
                                    0.05,
                                    0.9),
                          scoreFun = 'logLikelihood')

  # Generate data under topology H2
  h2_200_1_data <- h2(N = 200, ss = 1)

  # Adjacency matrix for topology H2
  adjMatrixH2 <- matrix(c(0, 1, 1, 0, 0,
                          0, 0, 0, 1, 0,
                          0, 0, 0, 0, 1,
                          0, 0, 0, 0, 1,
                          0, 0, 0, 0, 0),
                        byrow = TRUE,
                        nrow = 5)

  # Run the MH algorithm with the edges from the true graph.
  h2_200_1_mh <- MHEdge(adjMatrix = adjMatrixH2,
                        data = h2_200_1_data,
                        iterations = 500,
                        mutationRate = 1/5,
                        nGV = 0,
                        pmr = FALSE,
                        prior = c(0.05,
                                  0.05,
                                  0.9),
                        scoreFun = 'logLikelihood')

  expect_identical(m1ge_standard, m1ge_200_1_mh)
  expect_identical(h2_standard, h2_200_1_mh)

})
