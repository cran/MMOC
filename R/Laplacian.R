#' @title Calculate the graph Laplacian from a given adjacency matrix
#' @name Laplacian
#' @description Calculate the graph laplacian from a given kernel matrix that represents the full graph weighted adjacency matrix
#'
#'
#' @param A An n by n kernel matrix, where n is the sample size, that represents your initial adjacency matrix. Kernel matrices are symmetric, positive semi-definite distance matrices
#' @param laplacian One of `"shift"`, `"Ng"`, `"rw"` or `"sym"`. See details for description
#' @param grf.type Type of graph to calculate: `"full"` for adjacency matrix equal to the kernel, `"knn"` for a k-nearest neighbors graph, `"e-graph"` for an "epsilon graph"
#' @param k An integer value for `k` in the k-nearest neighbors graph. Only the `k` largest edges (most similar neighbors) will be kept
#' @param rho A value for the dispersion parameter in the Gaussian kernel. It is in the denominator of the exponent, so higher values correspond to lower similarity. By default it is the median pairwise Gaussian distance
#' @param epsilon The cutoff value for the `e-graph`. Edges lower than this value will be removed
#' @param mutual Make a "mutual" knn graph. Only keeps edges when two nodes are both in each others k-nearest set
#' @param binary.grf Set all edges >0 to 1
#' @param plots Whether or not to plot the final graph, a heatmap of calculated kernel, and the eigen values of the Laplacian
#' @details The four Lapalacians are defined as \eqn{L_{shift}=I+D^{-1/2}AD^{-1/2}}, \eqn{L_{Ng}=D^{-1/2}AD^{-1/2}}, \eqn{L_{sym}=I-D^{-1/2}AD^{-1/2}}, and \eqn{L_{rw}=I-D^{-1}A}. The shifted Laplacian, \eqn{L_{shift}=I+D^{-1/2}AD^{-1/2}}, is recommended for multi-view spectral clustering.
#' @return An n\eqn{\times}n matrix where `n` is the number of rows in `dat`.
#' @references \url{https://academic.oup.com/bioinformatics/article/36/4/1159/5566508#199177546}
#' @examples
#'
#' ## Generating data with 3 distinct clusters
#' ## Note that 'clustStruct' returns a list
#' dd <- clustStruct(n=120, p=30, k=3, noiseDat='random')[[1]]
#'
#' ## Gaussian kernel
#' rho <- median(dist(dd))
#' A <- exp(-(1/rho)*as.matrix(dist(dd, method = "euclidean", upper = TRUE)^2))
#'
#' Laplacian(A, laplacian='shift', grf.type = 'knn')
#' @export
Laplacian <- function(A, laplacian = c('shift', 'Ng', 'sym', 'rw'),
                      grf.type = c('full', 'knn', 'e-graph'),
                      k=5,
                      rho=NULL,
                      epsilon=0,
                      mutual=FALSE,
                      binary.grf=FALSE,
                      plots=TRUE){

  if(missing(laplacian)) laplacian <- 'shift'
  if(missing(grf.type)) grf.type <- 'full'

  diag(A) <- 0 ## Adjacency matrices need diag of 0

  if(grf.type == 'e-graph'){
    A[A < epsilon] <- 0
  }

  if(grf.type == 'knn'){
    A <- knnGrf(A, k, mutual=mutual)
  }

  if(binary.grf) A[A >0] <- 1

  deg <- rowSums(A)

  ## 'Shifted' Laplacian
  if(laplacian == "shift"){
    ds <- ifelse(deg>0, 1/sqrt(deg), 0)
    # L = I + D^(-1/2) A D^(-1/2)
    L <- diag(nrow = nrow(A)) + diag(ds) %*% A %*% diag(ds)
  }

  ## Laplacian used by Ng (2002) 'On Spectral Clustering'
  if(laplacian == "Ng"){
    ds <- ifelse(deg>0, 1/sqrt(deg), 0)
    # L = D^(-1/2) A D^(-1/2)
    L <- diag(ds) %*% A %*% diag(ds)
  }

  if(laplacian == 'sym'){
    ds <- ifelse(deg>0, 1/sqrt(deg), 0)
    # L = I - D^(-1/2) A D^(-1/2)
    L <- diag(nrow = nrow(A)) - diag(ds) %*% A %*% diag(ds)
  }

  ## Random walk Laplacian (best one)
  if(laplacian == 'rw'){
    ds <- ifelse(deg>0, 1/deg, 0)
    # L = I - D^-1 A
    L <- diag(nrow = nrow(A)) - diag(ds) %*% A
  }

  if(plots){
    stats::heatmap(A, symm=TRUE)

    plot(eigen(L)$values, col='dodgerblue3', pch=20, type="b",
         ylab="Eigen values", main = "Eigen values of the Laplacian")

    gg <- igraph::graph_from_adjacency_matrix(A, mode='undirected', weighted = T, diag=F)
    plot(gg, vertex.label=NA)
  }

  L
}


