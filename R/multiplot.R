#' @import grid
#'

# A function from the Cookbook for R website that prints multiple plots on one
# page.

multiplot <- function(...,
                      cols = 1) {

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...))

  numPlots = length(plots)

  if (numPlots == 1) {
    print(plots[[1]])

  } else {

    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))

    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout),
                                               ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain
      # this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]],
            vp = viewport(layout.pos.row = matchidx$row,
                          layout.pos.col = matchidx$col))
    }
  }
}
