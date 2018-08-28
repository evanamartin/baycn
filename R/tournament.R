#' selectBest
#'
#' Selects the individual with the highest log likelihood among a small subset
#' of the entire population.
#'
#' @param population A matrix where the first m columns are the DNA of the
#' individuals. Where m is the number of edges in the graph. The final column
#' is the log likelihood of the individual. The number of rows is the population
#' size.
#'
#' @param tournamentSize An integer. This is the number of individuals that will
#' be selected from the population and a winner, the individual with the highest
#' log likelihood, will be selected to be cloned and mutated. Then it will
#' replace the individual selected by the selectWorst function. In general the
#' tournament size is around one tenth of the population size.
#'
#' @param scoreFun A character string indicating what method to use for
#' calculating the fitness.
#'
#' @return A vector of length m. Where m is the number of edges in the graph.
#'
#' @export
#'
selectBest <- function (population,
                        tournamentSize,
                        scoreFun) {

  n <- nrow(population)
  m <- ncol(population)

  switch(scoreFun,

         'logLikelihood' = {

           # Select the location of an individual to start
           winner <- sample(x = 1:n,
                            size = 1)

           for (e in 1:tournamentSize) {

             # pick the location of a contender for the tournament
             contender <- sample(x = 1:n,
                                 size = 1)

             # replace the location of the winner with the location of the
             # contender if the contender has a higher log likelihood.
             if (population[contender, m] > population[winner, m]) {

               winner <- contender

             }

           }

           return (winner)

         },

         'BIC' = {

           # The individual who has a lower BIC is a better individual. Below I
           # change the winner from the individual who has the highest fitness
           # to the individual who has the lowest fitness.

           # Select the location of an individual to start
           winner <- sample(x = 1:n,
                            size = 1)

           for (e in 1:tournamentSize) {

             # pick the location of a contender for the tournament
             contender <- sample(x = 1:n,
                                 size = 1)

             # replace the location of the winner with the location of the
             # contender if the contender has a lower BIC
             if (population[contender, m] < population[winner, m]) {

               winner <- contender

             }

           }

           return (winner)

         })

}

#' selectWorst
#'
#' Selects the individual with the lowest log likelihood among a small subset
#' of the entire population.
#'
#' @param population A matrix where the first m columns are the DNA of the
#' individuals. Where m is the number of edges in the graph. The final column
#' is the log likelihood of the individual. The number of rows is the population
#' size.
#'
#' @param tournamentSize An integer. This is the number of individuals that will
#' be selected from the population and a loser, the individual with the lowest
#' log likelihood, will be replaced by the individual that was selected by the
#' selectBest function after the winner was cloned and mutated. In general the
#' tournament size is around one tenth of the population size.
#'
#' @param scoreFun A character string indicating what method to use for
#' calculating the fitness.
#'
#' @return A vector of length m. Where m is the number of edges in the graph.
#'
#' @export
#'
selectWorst <- function (population,
                         tournamentSize,
                         scoreFun) {

  n <- nrow(population)
  m <- ncol(population)

  switch (scoreFun,

    'logLikelihood' = {

      # Select the location of an individual to start
      loser <- sample(x = 1:n,
                      size = 1)

      for (e in 1:tournamentSize) {

        # pick the location of a contender for the tournament
        contender <- sample(x = 1:n,
                            size = 1)

        # replace the location of the loser with the location of the contender if
        # the contender has a lower fitness
        if (population[contender, m] < population[loser, m]) {

          loser <- contender

        }

      }

      return (loser)

    },

    'BIC' = {

      # The individual who has a higher BIC is a worse individual. Below I
      # change the loser from the individual who has the lowest fitness to the
      # individual who has the highest fitness.

      # Select the location of an individual to start
      loser <- sample(x = 1:n,
                      size = 1)

      for (e in 1:tournamentSize) {

        # pick the location of a contender for the tournament
        contender <- sample(x = 1:n,
                            size = 1)

        # replace the location of the loser with the location of the contender
        # if the contender has a higher BIC
        if (population[contender, m] > population[loser, m]) {

          loser <- contender

        }

      }

      return (loser)

    })

}
