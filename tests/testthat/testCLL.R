context('Log likelihood')

test_that('the likelihood functions calculate the log likelihood correctly',{

  # Simulate data for gn4, m2_gv, and m1_cph -----------------------------------

  set.seed(33)

  data_gn4 <- simdata(b0 = 0,
                      N = 200,
                      s = 1,
                      ss = 1,
                      graph = 'gn4')

  data_m2_gv <- simdata(b0 = 0,
                        N = 200,
                        s = 1,
                        ss = 1,
                        q = 0.45,
                        graph = 'm2_gv')

  data_m1_cph <- simdata(b0 = 0,
                         N = 200,
                         s = 1,
                         ss = 1,
                         p = 0.6,
                         q = 0.45,
                         graph = 'm1_cph')

  # Adjacency matrices ---------------------------------------------------------

  am_gn4 <- matrix(c(0, 1, 1, 0,
                     0, 0, 0, 1,
                     0, 0, 0, 0,
                     0, 0, 1, 0),
                   byrow = TRUE,
                   nrow = 4)

  diag(am_gn4) <- 0

  am_m2 <- matrix(1,
                  nrow = 3,
                  ncol = 3)

  diag(am_m2) <- 0

  am_m1_cph <- matrix(c(0, 1, 0, 0,
                        0, 0, 1, 1,
                        0, 0, 0, 0,
                        0, 0, 0, 0),
                      byrow = TRUE,
                      nrow = 4)

  # Calculate the log likelihood for gn4 ---------------------------------------

  ll_gn4 <- lookUp(data = data_gn4,
                   adjMatrix = am_gn4,
                   nCPh = 0,
                   nGV = 0,
                   nNodes = 4,
                   pmr = FALSE)

  # Log likelihood for the true graph.
  ll_tg_gn4 <- sum(ll_gn4$node1$dn0,
                   ll_gn4$node2$dn1,
                   ll_gn4$node3$dn9,
                   ll_gn4$node4$dn2)

  expect_equal(round(ll_tg_gn4, 3), -1121.318)

  # Calculate the log likelihood for m2_gv -------------------------------------

  ll_m2_gv <- lookUp(data = data_m2_gv,
                     adjMatrix = am_m2,
                     nCPh = 0,
                     nGV = 1,
                     nNodes = 3,
                     pmr = TRUE)

  # Log likelihood for the true graph.
  ll_tg_m2_gv <- sum(ll_m2_gv$node1$dn0,
                     ll_m2_gv$node2$dn5,
                     ll_m2_gv$node3$dn0)

  expect_equal(round(ll_tg_m2_gv, 3), -774.512)

  # Calculate the log likelihood for m1_cph ------------------------------------

  ll_m1_cph <- lookUp(data = data_m1_cph,
                      adjMatrix = am_m1_cph,
                      nCPh = 1,
                      nGV = 1,
                      nNodes = 4,
                      pmr = TRUE)

  # Log likelihood for the true graph.
  ll_tg_m1_cph <- sum(ll_m1_cph$node1$dn0,
                      ll_m1_cph$node2$dn1,
                      ll_m1_cph$node3$dn2,
                      ll_m1_cph$node4$dn2)

  expect_equal(round(ll_tg_m1_cph, 3), -841.734)

})
