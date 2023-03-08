#' Multiple runs of SpaRTaCo
#'
#' This function returns the estimated model parameters and the co-clustering labels obtained after running SpaRTaCo multiple times (parallel options available).
#'
#' @import SpatialExperiment
#' @import future.apply
#' @import future

#' @export
#'
#' @param data either a `SpatialExperiment` object or a matrix containing the experiment;
#' @param coordinates if `is.matrix(data)`, it takes the matrix of spatial coordinates of dimension `ncol(data)` x 2.
#' @param assay if `class(data) == "SpatialExperiment"`, it takes either the name or the index of the assay to be used;
#' @param K the number of row clusters (only when `input.values == NULL`);
#' @param R the number of column clusters (only when `input.values == NULL`);
#' @param nstart the number of parallel runs of the estimation algorithm;
#' @param Delta.constr the constraint on the Delta matrix (default is 10; see **Details**).
#' @param max.iter the maximum number of iterations the estimation algorithm is run.
#' @param estimate.iterations the maximum number of iterations within each M Step.
#' @param conv.criterion a list containing the parameters that define a converge criterion (see **Details**).
#' @param verbose logical; it `TRUE`, it displays the estimation process through a dynamic progress bar.
#' @param verbose.display.intervals integer; if `verbose == T`, the dynamic progress bar is updated after `max.iter/verbose.display.intervals` iterations.
#' @param compute.uncertainty if `TRUE` (default), it computes the clustering uncertainty of the rows and of the columns.
#'
#' @return An object of class `spartaco` with the parameter estimates, the clustering labels, the log-likelihood value at each iteration and the data, the ICL, the data matrix and the coordinates matrix, and the clustering uncertainty.
#'
#' @details This function allows to run the `spartaco` model starting from multiple starting points simultaneously.
#' It is possible to run this function using multiple cores; to do so, use the `multicore` function package `future` (see **Examples**).
#'
#' If `verbose == T`, it displays the estimation process through a progress bar. Note that in this case the final output will be based just on the last `max.iter/verbose.display.intervals` iterations.
#' For details about the rest of the parameters, check [spartaco::spartaco()]

#' @examples
#' library(spartaco)
#'
#' # First, create the data matrix:
#' n <- p <- 300
#' K <- R <- 3
#' x <- matrix(runif(n*p), n, p)
#' coordinates <- matrix(runif(2*p), p, 2)
#'
#' # Set the number of cores to be used for the computations. In this example, we use 3 cores.
#' future::plan(future::multisession(workers = 3))
#' output <- spartaco(data = x, coordinates = coordinates, K = K, R = R, max.iter = 1000, verbose.display.intervals = 100)
#'
#' # according to this setup, the progress bar will be updated once every 10 iterations are performed.

spartacoRC_multirun <- function(data,
                              coordinates = NULL,
                              assay = NULL,
                              column.labels,
                              K = NULL,
                              nstart = 5,
                              Delta.constr = 10,
                              max.iter = 1000,
                              estimate.iterations = 100,
                              conv.criterion = list(iterations = 10, epsilon = 1e-4),
                              verbose = T,
                              verbose.display.intervals = 10,
                              compute.uncertainty = TRUE,
                              lambda.ridge = 0,
                              lambda.lasso = 0)
{
  if(class(data)[1] == "SpatialExperiment"){
    if(is.numeric(assay)) which.assay <- assay
    else which.assay <- which(names(data@assays@data) == assay)
    x <- as.matrix(data@assays@data[[which.assay]])
    row.names(x) <- rowData(data)$gene_name
    coordinates <- as.matrix(spatialCoords(data))
  } else {
    x <- data
  }


  if(verbose){
    progressr::handlers(global = T)
    P <- progressr::progressor(along = 1:max.iter)
    results <- future_lapply(1:nstart, FUN = function(l){
      for(j in 1:verbose.display.intervals){
        P()
        if(j == 1) input.values <- NULL else input.values <- s
        s <- RCspartaco(data = x, coordinates = coordinates, column.labels, K = K, assay = NULL,
                      input.values = input.values,
                      Delta.constr = Delta.constr, max.iter = max.iter/verbose.display.intervals,
                      estimate.iterations = estimate.iterations, conv.criterion = conv.criterion,
                      verbose = F, save.options = NULL,
                      lambda.ridge = lambda.ridge,
                      lambda.lasso = lambda.lasso)
      }
      s
    }
    )
  }
  if(!verbose){
    results <- future_lapply(1:nstart, FUN = function(l)
      RCspartaco(data = x, coordinates = coordinates, column.labels, K = K, assay = NULL,
               Delta.constr = Delta.constr, max.iter = max.iter,
               estimate.iterations = estimate.iterations, conv.criterion = conv.criterion,
               verbose = F, save.options = NULL,
               lambda.ridge = lambda.ridge,
               lambda.lasso = lambda.lasso)
    )
  }
  output <- CombineSpartacoRC(results, compute.uncertainty)
  return(output)
}
