% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PlotSpartaco.R
\name{plot.spartaco}
\alias{plot.spartaco}
\title{Plot SpaRTaCo
This function returns the ggplots of the mean and the spatial signal-to-noise ratios from a SpaRTaCo model.}
\usage{
\method{plot}{spartaco}(
  x,
  type = 1,
  gene.name = readline("gene name: "),
  k = NULL,
  r = 1:ncol(x$mu),
  manual.palette = NULL,
  display.all.spots = T,
  remove.outliers = F,
  coef = 1,
  range = NULL,
  ...
)
}
\arguments{
\item{x}{a `spartaco` object;}

\item{type}{the type of plot you want to show.
- `1` displays the co-cluster mean levels;
- `2` displays the co-cluster spatial signal-to-noise ratios;
- `3` displays the map of the spots colored with respect to the estimated column clusters;
- `4` if `!is.null(k)`, it displays the average gene expression in the clusters `k` and `r`, otherwise it displays the expression of the gene given in `gene.name`.}

\item{gene.name}{(used only when `type == 4`); it receives the name of the gene to display. If `!is.null(k)`, it displays the average expression of the gene clusters given by `k`.}

\item{k}{(used only when `type == 4`) the gene clusters to plot.}

\item{r}{(used only when `type == 4`) the spot clusters to plot.}

\item{manual.palette}{a vector of colors used when `type == 3`.}

\item{display.all.spots}{if `TRUE` (default) and `type == 4`, it displays the entire grid of spots.}

\item{remove.outliers}{(used only when `type == 4)` if TRUE remove spots with extreme gene expression value (lower-whisker or upper-whisker)}

\item{coef}{(used only when `type == 4` and if `remove.outliers == TRUE`) in the selection of outliers, spots with gene expression value grater that upper-whisker value times coef are removed.}

\item{range}{vector to set the minimum and maximum of the scale of value to display (the default is `NULL` and the minimum and maximum of the value is used).}
}
\value{
The requested plot is displayed. In addition, if assigned to an object, it will return the `ggplot` object.
}
\description{
Plot SpaRTaCo
This function returns the ggplots of the mean and the spatial signal-to-noise ratios from a SpaRTaCo model.
}
