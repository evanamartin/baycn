context('MCMC no genetic variants')

test_that('the MHEdge function returns the correct matrix',{

  # Load the standards to test against.
  load(system.file('testdata',
                   'standards.RData',
                   package = 'baycn'))

  set.seed(338)

  # Generate data for topology M1 with no genetic variants.
  data_m1ge_200_1 <- m1ge(N = 200, ss = 1)

  # Adjacency matrix for toplogy M1
  adjMatrix_m1 <- matrix(c(0, 1, 0,
                           0, 0, 1,
                           0, 0, 0),
                         byrow = TRUE,
                         nrow = 3)

  # Run the MH algorithm with the edges from the true graph.
  mh_m1ge_200_1 <- mhEdge(adjMatrix = adjMatrix_m1,
                          data = data_m1ge_200_1,
                          iterations = 100,
                          mutationRate = 1/2,
                          nGV = 0,
                          pmr = FALSE,
                          prior = c(0.05,
                                    0.05,
                                    0.9),
                          scoreFun = 'logLikelihood')

  # Generate data under topology H2
  data_h2_200_1 <- h2(N = 200, ss = 1)

  # Adjacency matrix for topology H2
  adjMatrix_h2 <- matrix(c(0, 1, 1, 0, 0,
                           0, 0, 0, 1, 0,
                           0, 0, 0, 0, 1,
                           0, 0, 0, 0, 1,
                           0, 0, 0, 0, 0),
                         byrow = TRUE,
                         nrow = 5)

  # Run the MH algorithm with the edges from the true graph.
  mh_h2_200_1 <- mhEdge(adjMatrix = adjMatrix_h2,
                        data = data_h2_200_1,
                        iterations = 100,
                        mutationRate = 1/5,
                        nGV = 0,
                        pmr = FALSE,
                        prior = c(0.05,
                                  0.05,
                                  0.9),
                        scoreFun = 'logLikelihood')

  expect_identical(standard_m1ge, mh_m1ge_200_1)
  expect_identical(standard_h2, mh_h2_200_1)

})
