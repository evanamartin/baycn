context('baycn without genetic variants')

test_that('the mhEdge function infers the correct graph',{

  # Adjacency matrix for M1 - M4 -----------------------------------------------

  # Fully connected adjacency matrix for M1 - M4
  am_m <- matrix(c(0, 1, 1,
                   0, 0, 1,
                   0, 0, 0),
                 byrow = TRUE,
                 nrow = 3)

  # Run baycn on data simulated for M2 -----------------------------------------

  set.seed(5)

  data_m2 <- simdata(b0 = 0,
                     N = 500,
                     s = 1,
                     ss = 1,
                     graph = 'm2_ge')

  baycn_m2 <- mhEdge(adjMatrix = am_m,
                     burnIn = 0.2,
                     data = data_m2,
                     iterations = 1000,
                     nGV = 0,
                     pmr = FALSE,
                     prior = c(0.05,
                               0.05,
                               0.9),
                     progress = FALSE,
                     thinTo = 200)

  # Run baycn on data simulated for M1 -------------------------------------------

  data_m1 <- simdata(b0 = 0,
                     N = 500,
                     s = 1,
                     ss = 1,
                     graph = 'm1_ge')

  baycn_m1 <- mhEdge(adjMatrix = am_m,
                     burnIn = 0.2,
                     data = data_m1,
                     iterations = 1000,
                     nGV = 0,
                     pmr = FALSE,
                     prior = c(0.05,
                               0.05,
                               0.9),
                     progress = TRUE,
                     thinTo = 200)

  # Calculate the MSE for M1 and M2 --------------------------------------------

  # Expected probabilities for M2
  ep_m2 <- matrix(c(1, 0, 0,
                    0, 0, 1,
                    0, 1, 0),
                  byrow = TRUE,
                  nrow = 3)

  mse_m2 <- sum((baycn_m2@posteriorES[, 2:4] - ep_m2)^2)

  # Expected probabilities for M1
  ep_m1 <- matrix(c(0.33, 0.66, 0,
                    0,    0,    1,
                    0.66, 0.33, 0),
                  byrow = TRUE,
                  nrow = 3)

  mse_m1 <- sum((baycn_m1@posteriorES[, 2:4] - ep_m1)^2)

  expect_true(mse_m2 < 0.1)
  expect_true(mse_m1 < 0.1)

})
