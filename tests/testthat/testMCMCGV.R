context('MCMC with genetic variants')

test_that('MHEdge returns the correct matrix with and without pmr',{

  # Load the standards to test against.
  load(system.file('testdata',
                   'standardsGV.RData',
                   package = 'baycn'))

  set.seed(22)

  # Generate data under topology M1 with one genetic variant.
  data_m1gv_200_1 <- simdata(b0 = 0,
                             N = 200,
                             s = 1,
                             graph = 'm1_gv',
                             ss = 1,
                             p = 0.45)

  # Adjacency matrix for toplogy M1
  am_m1gv <- matrix(c(0, 1, 0,
                      0, 0, 1,
                      0, 0, 0),
                    byrow = TRUE,
                    nrow = 3)

  # Run the MH algorithm with the edges from the true graph with out pmr.
  mh_m1gv_200_1 <- mhEdge(adjMatrix = am_m1gv,
                          burnIn = 0,
                          data = data_m1gv_200_1,
                          iterations = 100,
                          nGV = 1,
                          pmr = FALSE,
                          prior = c(0.05,
                                    0.05,
                                    0.9),
                          progress = TRUE,
                          thinTo = 100)

  # Run the MH algorithm with the edges from the true graph with pmr.
  mh_m1gv_200_1_pmr <- mhEdge(adjMatrix = am_m1gv,
                              burnIn = 0,
                              data = data_m1gv_200_1,
                              iterations = 100,
                              nGV = 1,
                              pmr = TRUE,
                              prior = c(0.05,
                                        0.05,
                                        0.9),
                              progress = FALSE,
                              thinTo = 100)

  expect_identical(standard_m1gv@chain, mh_m1gv_200_1@chain)
  expect_identical(standard_m1gv_pmr@chain, mh_m1gv_200_1_pmr@chain)

  expect_equal(standard_m1gv@likelihood, mh_m1gv_200_1@likelihood)
  expect_equal(standard_m1gv_pmr@likelihood, mh_m1gv_200_1_pmr@likelihood)

})
