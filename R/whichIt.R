whichIt <- function (burnIn,
                     iterations,
                     thinTo) {

  # Determine which iterations will be kept.
  if (burnIn == 0) {

    return (ceiling(seq(from = 1,
                        to = iterations,
                        length.out = thinTo)))

  } else {

    if ((1 - burnIn) * iterations < thinTo) {

      return ('Error: The number of iterations left after the burn-in is less than the thinTo value.')

    }

    return (ceiling(seq(from = burnIn * iterations,
                        to = iterations,
                        length.out = thinTo)))

  }

}
