context('MCMC with genetic variants')

test_that('MHEdge returns the correct matrix with and without pmr',{

  # Load the standards to test against.
  load(system.file('testdata',
                   'standardsGV.RData',
                   package = 'baycn'))

  set.seed(22)

  # Generate data under topology M1 with one genetic variant.
  data_m1gv_200_1 <- m1gv(N = 200,
                          p = 0.45,
                          ss = 1)

  # Adjacency matrix for toplogy M1
  adjMatrix_m1gv <- matrix(c(0, 1, 0,
                             0, 0, 1,
                             0, 0, 0),
                           byrow = TRUE,
                           nrow = 3)

  # Run the MH algorithm with the edges from the true graph with out pmr.
  mh_m1gv_200_1 <- mhEdge(adjMatrix = adjMatrix_m1gv,
                          data = data_m1gv_200_1,
                          iterations = 100,
                          nGV = 1,
                          pmr = FALSE,
                          prior = c(0.05,
                                    0.05,
                                    0.9))

  # Run the MH algorithm with the edges from the true graph with pmr.
  mh_m1gv_200_1_pmr <- mhEdge(adjMatrix = adjMatrix_m1gv,
                              data = data_m1gv_200_1,
                              iterations = 100,
                              nGV = 1,
                              pmr = TRUE,
                              prior = c(0.05,
                                        0.05,
                                        0.9))

  # Generate data under topology M3 with one genetic variant.
  data_m3gv_200_1 <- m3gv(N = 200,
                          p = 0.45,
                          ss = 1)

  # Adjacency matrix for topology M3
  adjMatrix_m3gv <- matrix(c(0, 1, 1,
                             0, 0, 0,
                             0, 0, 0),
                           nrow = 3,
                           byrow = TRUE)

  # Run the MH algorithm with the edges from the true graph without pmr.
  mh_m3gv_200_1 <- mhEdge(adjMatrix = adjMatrix_m3gv,
                          data = data_m3gv_200_1,
                          iterations = 100,
                          nGV = 1,
                          pmr = FALSE,
                          prior = c(0.05,
                                    0.05,
                                    0.9))

  # Run the MH algorithm with the edges from the true graph with pmr.
  mh_m3gv_200_1_pmr <- mhEdge(adjMatrix = adjMatrix_m3gv,
                              data = data_m3gv_200_1,
                              iterations = 100,
                              nGV = 1,
                              pmr = TRUE,
                              prior = c(0.05,
                                        0.05,
                                        0.9))

  expect_identical(standard_m1gv, mh_m1gv_200_1)
  expect_identical(standard_m1gv_pmr, mh_m1gv_200_1_pmr)
  expect_identical(standard_m3gv, mh_m3gv_200_1)
  expect_identical(standard_m3gv_pmr, mh_m3gv_200_1_pmr)

})
