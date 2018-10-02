context('MCMC with genetic variants')

test_that('MHEdge returns the correct matrix with and without pmr',{

  # Load the standards to test against.
  load(system.file('testdata',
                   'standardsGV.RData',
                   package = 'edgeFrequency'))

  set.seed(22)

  # Generate data under topology M1 with one genetic variant.
  m1gv_200_1_data <- m1gv(N = 200,
                          p = 0.45,
                          ss = 1)

  # Adjacency matrix for toplogy M1
  m1gv_adjMatrix <- matrix(c(0, 1, 0,
                             0, 0, 1,
                             0, 0, 0),
                           byrow = TRUE,
                           nrow = 3)

  # Run the MH algorithm with the edges from the true graph with out pmr.
  m1gv_200_1_mh <- MHEdge(adjMatrix = m1gv_adjMatrix,
                          data = m1gv_200_1_data,
                          iterations = 100,
                          mutationRate = 1/2,
                          nGV = 1,
                          pmr = FALSE,
                          prior = c(0.05,
                                    0.05,
                                    0.9),
                          scoreFun = 'logLikelihood')

  # Run the MH algorithm with the edges from the true graph with pmr.
  m1gv_200_1_mh_pmr <- MHEdge(adjMatrix = m1gv_adjMatrix,
                              data = m1gv_200_1_data,
                              iterations = 100,
                              mutationRate = 1/2,
                              nGV = 1,
                              pmr = TRUE,
                              prior = c(0.05,
                                        0.05,
                                        0.9),
                              scoreFun = 'logLikelihood')

  # Generate data under topology M3 with one genetic variant.
  m3gv_200_1_data <- m3gv(N = 200,
                          p = 0.45,
                          ss = 1)

  # Adjacency matrix for topology M3
  m3gv_adjMatrix <- matrix(c(0, 1, 1,
                             0, 0, 0,
                             0, 0, 0),
                           nrow = 3,
                           byrow = TRUE)

  # Run the MH algorithm with the edges from the true graph without pmr.
  m3gv_200_1_mh <- MHEdge(adjMatrix = m3gv_adjMatrix,
                          data = m3gv_200_1_data,
                          iterations = 100,
                          mutationRate = 1/2,
                          nGV = 1,
                          pmr = FALSE,
                          prior = c(0.05,
                                    0.05,
                                    0.9),
                          scoreFun = 'logLikelihood')

  # Run the MH algorithm with the edges from the true graph with pmr.
  m3gv_200_1_mh_pmr <- MHEdge(adjMatrix = m3gv_adjMatrix,
                              data = m3gv_200_1_data,
                              iterations = 100,
                              mutationRate = 1/2,
                              nGV = 1,
                              pmr = TRUE,
                              prior = c(0.05,
                                        0.05,
                                        0.9),
                              scoreFun = 'logLikelihood')

  expect_identical(m1gv_standard, m1gv_200_1_mh)
  expect_identical(m1gv_standard_pmr, m1gv_200_1_mh_pmr)
  expect_identical(m3gv_standard, m3gv_200_1_mh)
  expect_identical(m3gv_standard_pmr, m3gv_200_1_mh_pmr)

})
