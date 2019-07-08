context('MCMC no genetic variants')

test_that('the MHEdge function returns the correct matrix',{

  # Load the standards to test against.
  load(system.file('testdata',
                   'standards.RData',
                   package = 'baycn'))

  set.seed(338)

  # Generate data for topology M1 with no genetic variants.
  data_m1ge_200_1 <- simdata(b0 = 1,
                             N = 200,
                             s = 1,
                             graph = 'm1_ge',
                             ss = 1)

  # Adjacency matrix for toplogy M1
  am_m1 <- matrix(c(0, 1, 0,
                    0, 0, 1,
                    0, 0, 0),
                  byrow = TRUE,
                  nrow = 3)

  # Run the MH algorithm with the edges from the true graph.
  mh_m1ge_200_1 <- mhEdge(adjMatrix = am_m1,
                          burnIn = 0,
                          data = data_m1ge_200_1,
                          iterations = 100,
                          nGV = 0,
                          pmr = FALSE,
                          prior = c(0.05,
                                    0.05,
                                    0.9),
                          thinTo = 100)

  # Generate data under topology H2
  data_h2_200_1 <- simdata(b0 = 1,
                           N = 200,
                           s = 1,
                           graph = 'gn5',
                           ss = 1)

  # Adjacency matrix for topology H2
  am_h2 <- matrix(c(0, 1, 1, 0, 0,
                    0, 0, 0, 1, 0,
                    0, 0, 0, 0, 1,
                    0, 0, 0, 0, 1,
                    0, 0, 0, 0, 0),
                  byrow = TRUE,
                  nrow = 5)

  # Run the MH algorithm with the edges from the true graph.
  mh_h2_200_1 <- mhEdge(adjMatrix = am_h2,
                        burnIn = 0,
                        data = data_h2_200_1,
                        iterations = 100,
                        nGV = 0,
                        pmr = FALSE,
                        prior = c(0.05,
                                  0.05,
                                  0.9),
                        thinTo = 100)

  expect_identical(standard_m1ge@chain, mh_m1ge_200_1@chain)
  expect_identical(standard_h2@chain, mh_h2_200_1@chain)

  expect_equal(standard_m1ge@likelihood, mh_m1ge_200_1@likelihood)
  expect_equal(standard_h2@likelihood, mh_h2_200_1@likelihood)

})
