my_colors <- brewer.pal(5, "Spectral")
my_colors <- colorRampPalette(my_colors)(100)

#' Map of signature scores over the annotation areas of the tissue
#'
#' @param data summarized experiment object with signatures scores by row, 
#' coordinates and column.label in colData
#' @param whichSign name of the signature of interest
#'
#' @return ggplot object
#' @exports
#'
#' @examples
plotSign <- function(data, whichSign, size = 1.5){
  scores <- assay(data)[rownames(data) == whichSign,]
  Coord <- as.data.frame(spatialCoords(data))
  colnames(Coord) <- c("y","x")
  Coord$group <- colData(data)$column.label
  tmp <- cbind(scores,Coord)
  
  ggplot(tmp, aes(x,y)) +
    geom_point(data = select(tmp, -group), shape = "circle", size = size/0.8, col = "#E1CCCA", alpha =0.2) +
    geom_point(data = tmp, mapping = aes(x, y,col = scores), shape = "circle", size = size) +
    scale_color_gradient2(low = "#2B3CBA",
                            mid = "#FEFDBD",
                            high = "#D7191C",
                            midpoint = 0, 
                            n.breaks = 7) +
    scale_x_continuous(trans = "reverse") +
    coord_flip() +
    xlab("") + ylab("") +
    ggtitle(whichSign) +
    facet_wrap(vars(group)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "#F3F0F1"),
          axis.text.x = element_blank(),#element_text(size=18),
          axis.text.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title = element_blank(),#element_text(size=18),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 8, face = "plain"),
          plot.title = element_text(hjust = 0.5, size = 12),
          title = element_text(size=8, face = "bold"),
          strip.text = element_text(size=10, face = "bold"),
          plot.margin= grid::unit(c(3,2,3,2), "mm"))
}
