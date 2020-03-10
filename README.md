<img src='man/figures/logo.png' align="right" height="139" /></a>

[![CRAN status](https://www.r-pkg.org/badges/version/baycn)](https://cran.r-project.org/package=baycn)
[![Build Status](https://travis-ci.org/evanamartin/baycn.svg?branch=master)](https://travis-ci.org/evanamartin/baycn)

# baycn

baycn is a Bayesian hybrid approach for inferring Directed Acyclic Graphs (DAGs) for continuous, discrete, and mixed data. The algorithm can use the graph inferred by another more efficient graph inference method as input; the input graph may contain false edges or undirected edges, but can help reduce the search space to a more manageable size. A Bayesian Metropolis-Hastings Markov chain Monte Carlo algorithm is then used to infer the probability of direction and absence for the edges in the network.

## Installation

``` r
# The released version of baycn can be installed from CRAN:
install.packages('baycn')

# The development (and most up to date) version can be installed from GitHub
# using devtools:
# If devtools is not installed first run install.packages('devtools')
devtools::install_github('evanamartin/baycn')
```

## Examples

Here we show a few toy examples to demonstrate how to use baycn to infer a network and interpret the output.

### Continuous data

The following R code simulates and analyzes data from topology GN4 which has four continuous random variables. The following figure shows the true graph for topology GN4.

<img src='man/figures/gn4.jpg' align="center" width="100" />

This graph has two other Markov equivalent graphs, which form a Markov equivalence class. The graphs in the Markov equivalence class for GN4 are shown in the following figure. 

<img src='man/figures/gn4_mec.jpg' align="center" width="300" />

The edges in orange denote the edges whose direction cannot be uniquely determined due to Markov equivalence. We can calculate the expected probability for the orange edges by counting the number of times they appear for each direction and dividing by the number of Markov equivalent graphs. For example, we expect the T1-T2 edge to point from T1 to T2 1/3 of the time and the T2-T4 edge to point from T2 to T4 2/3 of the time.

```r
set.seed(5)

# Generate data for topology GN4.
data_gn4 <- simdata(b0 = 0,
                         N = 600,
                         s = 1,
                         ss = 1,
                         graph = 'gn4')

# Fully connected adjacency matrix for topology GN4.
am_gn4 <- matrix(c(0, 1, 1, 1,
                   1, 0, 1, 1,
                   1, 1, 0, 1,
                   1, 1, 1, 0),
                byrow = TRUE,
                nrow = 4)

# Run baycn on the simulated data.
baycn_gn4 <- mhEdge(adjMatrix = am_gn4,
                  burnIn = 0.2,
                  data = data_gn4,
                  iterations = 30000,
                  nCPh = 0,
                  nGV = 0,
                  pmr = FALSE,
                  prior = c(0.05,
                            0.05,
                            0.9),
                  progress = TRUE,
                  thinTo = 200)
```

Below is the output from baycn on graph GN4 with a fully connected graph used as the input. Edges 1, 2, 5, and 6 are the true edges used to simulate the data. From the output we can see that edge 1 points from T1 to T2 (column 'zero') about 1/3 of the time and edge 5 points from T2 to T4 (column 'zero') about 2/3 of the time which is what we expected to see based on the Markov equivalent graphs for GN4. In addition baycn infers edges 2 and 6 to be oriented in the correct direction almost 100 percent of the time and edges 3 and 4 to be absent (column 'two') with a high probability which is consistent with the true graph. 

```r
# Display a summary of the output from baycn.
summary(baycn_gn4)

#< Posterior probability: 
#<       nodes  zero   one   two
#< edge1 T1-T2 0.380 0.620 0.000
#< edge2 T1-T3 0.980 0.020 0.000
#< edge3 T1-T4 0.115 0.070 0.815
#< edge4 T2-T3 0.050 0.005 0.945
#< edge5 T2-T4 0.675 0.325 0.000
#< edge6 T3-T4 0.030 0.970 0.000

#< Log likelihood: 
#<    Min    1Q Median    3Q      Max
#<  -3410 -3410  -3410 -3410 -3409.39

#< Number of unique graphs: 3
#< Run time in seconds: 5.974201
#< Iterations: 30000
#< Burn in: 20%
#< Step size: 120
```

### Mixed data

The following R code simulates and analyzes data from topology M1 which has one discrete random variable U and two continuous random variables. In this example the variable U represents a genetic variant and the two T nodes represent gene expression variables. For this analysis we will make the assumption that a gene expression node cannot be the parent of a genetic variant node by setting pmr = TRUE. The pmr argument can be set to TRUE anytime it is reasonable to assume the continuous variables cannot be the parents of the discrete variables (see the documentation for the mhEdge function). The following figure shows the true graph for topology M1.

<img src='man/figures/m1_gv.jpg' align="center" width="100" />

Because we assume that the T nodes cannot be the parents of the U node this graph does not have any other Markov equivalent graphs. We expect the probability of the U-T1 edge to always be oriented U -> T1 and the T1-T2 edge to always be oriented T1 -> T2.

```r
set.seed(72)

# Simulate data for topology M1 with one genetic variant.
data_m1_200_1 <- simdata(b0 = 0,
                         N = 200,
                         s = 1,
                         ss = 1,
                         graph = 'm1_gv')

# Fully connected adjacency matrix for topology M1.
am_m1 <- matrix(c(0, 1, 1,
                  0, 0, 1,
                  0, 0, 0),
                byrow = TRUE,
                nrow = 3)

# Run baycn on the simulated data with pmr = TRUE.
baycn_m1_gv <- mhEdge(adjMatrix = am_m1,
                  burnIn = 0.2,
                  data = data_m1_200_1,
                  iterations = 20000,
                  nCPh = 0,
                  nGV = 1,
                  pmr = TRUE,
                  prior = c(0.05,
                            0.05,
                            0.9),
                  progress = TRUE,
                  thinTo = 200)
```

Below is the output from baycn for graph M1 with a fully connected graph used as the input. The posterior probability estimated by baycn for each edge and edge state is close to what we expect to see. For this analysis we expect edges 1 and 3 to have a high posterior for state zero because M1 does not have any Markov equivalent graphs. In addition, we do not expect to see any posterior probability for state one for either of these edges because of the restriction placed on them (pmr = TRUE). Finally, for edge 2 we expect the posterior to be close to 1 for state two because this edge was not used in simulating the data. 

```r
# Display a summary of the output from baycn.
summary(baycn_m1_gv)

#< Posterior probability: 
#<       nodes  zero  one   two
#< edge1  U-T1 1.000 0.00 0.000
#< edge2  U-T2 0.145 0.00 0.855
#< edge3 T1-T2 0.910 0.09 0.000

#< Log likelihood: 
#<      Min      1Q  Median      3Q     Max
#<  -665.22 -665.22 -665.22 -665.22 -665.22

#< Number of unique graphs: 3
#< Run time in seconds: 3.456078
#< Iterations: 20000
#< Burn in: 20%
#< Step size: 80
```
