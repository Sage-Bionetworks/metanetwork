#' Multiple plot function
#'
#' This function plots multiple ggplot objects passed in as either ..., or to plotlist
#' (as a list of ggplot objects). 
#' 
#' @param ... Required. Multiple ggplot objects, or a plotlist as a list of 
#' ggplot objects
#' @param plotlist Optional. A list of ggplot objects. (Default = NULL)
#' @param file Optional. A file name to save the plot as, currently non-functional.
#' (Default=NULL)
#' @param cols Optional. (Default = 1)
#' @param layout Optional. The plot layout as a matrix. eg. if 
#' `layout = matrix(c(1,2,3,3), nrow=2, byrow=TRUE)`, then plot 1 will go in the 
#' upper left, 2 will go in the upper right, and 3 will go all the way across 
#' the bottom. 
#'
#' @return A ggplot image
#' 
#' @export
multiplot <- function(..., plotlist=NULL, file = NULL, cols=1, layout=NULL) {
  #library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = grid::viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}