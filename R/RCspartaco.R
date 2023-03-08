#' SpaRTaCo with Column Cluster assigned
#'
#' This function returns the estimated model parameters and the row-clustering labels, given column-cluster labels.
#'
#' @import SpatialExperiment
#' @export
#'
#' @param data either a `SpatialExperiment` object or a matrix containing the experiment;
#' @param coordinates if `is.matrix(data)`, it takes the matrix of spatial coordinates of dimension `ncol(data)` x 2.
#' @param assay if `class(data) == "SpatialExperiment"`, it takes either the name or the index of the assay to be used;
#' @param column.labels is a vector of length `ncol(data)` containing the column clustering labels.
#' @param K the number of row clusters (only when `input.values == NULL`).
#' @param Delta.constr the constraint on the Delta matrix (default is 10; see **Details**).
#' @param max.iter the maximum number of iterations the estimation algorithm is run.
#' @param estimate.iterations the maximum number of iterations within each M Step.
#' @param input.values the starting points of the estimation process (see **Details**). If passed, the values of `K` is taken from it. The output of a previous model estimation can be passed here.
#' @param conv.criterion a list containing the parameters that define a converge criterion (see **Details**).
#' @param verbose different ways to display the on-going estimation process (see **Details**).
#' @param verbose.parallel.label an additional label that is used to display the on-going estimation process when `verbose == "parallel"`.
#' @param save.options a list for specifying the saving parameters (see **Details**).
#' @param seed set the interval seed of the function.
#'
#'
#' @return An object of class `spartaco` with the parameter estimates, the clustering labels, the log-likelihood value at each iteration and the data, the ICL, the data matrix and the coordinates matrix.

#'
#' @details `Delta.constr` gives the quantity \deqn{c = \tau_{kr} + \xi_{kr},} where \eqn{\tau_{kr}} and \eqn{\xi_{kr}} are the spatial variance and the nugget effect of block \eqn{(k,r)}.
#'
#' The algorithm can be initiated from a given set of starting values. To do so, `input.values` receives a list of the form
#' `list(mu, tau, xi, alpha, beta, phi, Cs)`, where:
#' - `mu`, `tau`, `xi`, `alpha` and `beta` are `K` x `R` matrices;
#' - `phi` is a vector of length `R`;
#' - `Cs` is a vector of length `nrow(data)` containing the row clustering labels;
#'
#' If the algorithm is initiated from some starting values,  `K` is set automatically according to the input values.
#' If an object of class `spartaco` is passed to `input.values`, the estimation starts from the final estimate of the previous run (see **Examples**).
#'
#' If `conv.criterion == NULL`, the algorithm is stopped after `max.iter` itereations, otherwise it is stopped when the increment of the log-likelihood is smaller than a certain threshold `conv.criterion$epsilon` for `conv.criterion$iterations` times in a row.
#'
#' The function allows also to save the results even if the estimation is not completed. To do so, `save.options` receives a list of two parameters:
#' `after` gives the number of iterations after which the results are saved, `file.name` contains the path where the results are saved.
#'
#' If `verbose == "full"`, the on-going estimation procedure is displayed. If `verbose == "progress"`, a dynamic progress bar will display the percentage of iterations completed.
#' If `verbose == F`, then nothing is displayed in console.

# \deqn{p(x) = \frac{\lambda^x e^{-\lambda}}{x!}}{%p(x) = \lambda^x exp(-\lambda)/x!} for \eqn{x = 0, 1, 2, \ldots}

RCspartaco <- function(data,
                     coordinates = NULL,
                     assay = NULL,
                     column.labels,
                     K = NULL,
                     Delta.constr = 10,
                     max.iter = 1000,
                     estimate.iterations = 100,
                     input.values = NULL,
                     conv.criterion = list(iterations = 10, epsilon = 1e-4),
                     verbose = T,
                     verbose.parallel.label = NULL,
                     save.options = NULL,
                     seed = NULL,
                     lambda.ridge = 0,
                     lambda.lasso = 0
) {
  if(class(data)[1] == "SpatialExperiment"){
    if(is.numeric(assay)) which.assay <- assay
    else which.assay <- which(names(data@assays@data) == assay)
    x <- as.matrix(data@assays@data[[which.assay]])
    row.names(x) <- rowData(data)$gene_name
    coordinates <- as.matrix(spatialCoords(data))
  } else {
    x <- data
  }

  set.seed(seed = seed)

  if(!is.null(input.values)){
    K <- nrow(input.values$mu)
  }

  main(x = x,
       coordinates = coordinates,
       K = K,
       column.labels = column.labels,
       Delta.constr = Delta.constr,
       max.iter = max.iter,
       estimate.iterations = estimate.iterations,
       conv.criterion = conv.criterion,
       input.values = input.values,
       save.options = save.options,
       verbose = verbose,
       verbose.parallel.label = verbose.parallel.label,
       lambda.ridge = lambda.ridge,
       lambda.lasso = lambda.lasso)
}
