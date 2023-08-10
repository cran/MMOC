#' @title Calculate the graph Laplacian of a given data set
#' @name kernelLaplacian
#' @description Calculate the graph laplacian from a given data set with subjects as rows and features as columns.
#'
#' @importFrom Spectrum CNN_kernel
#'
#' @param dat A matrix like object with subjects as rows and features as columns.
#' @param kernel The type of kernel used to calculate the graph's adjacency matrix: `"Gaussian"` for the standard Gaussian kernel, `"ZM"` for the Zelnik-Manor kernel, `"Spectrum"` for the spectrum kernel, `"Linear"` for the linear kernel (dot product), and `"Cor"` for a kernel of pairwise correlations. See references for more details.
#' @param laplacian One of `"shift"`, `"Ng"`, `"rw"` or `"sym"`. See details for description
#' @param grf.type Type of graph to calculate: `"full"` for adjacency matrix equal to the kernel, `"knn"` for a k-nearest neighbors graph, `"e-graph"` for an "epsilon graph"
#' @param k An integer value for `k` in the k-nearest neighbors graph. Only the `k` largest edges (most similar neighbors) will be kept
#' @param p An integer value for the p-nearest neighbor in the `ZM` kernel
#' @param rho A value for the dispersion parameter in the Gaussian kernel. It is in the denominator of the exponent, so higher values correspond to lower similarity. By default it is the median pairwise Gaussian distance
#' @param epsilon The cutoff value for the `e-graph`. Edges lower than this value will be removed
#' @param mutual Make a "mutual" knn graph. Only keeps edges when two nodes are both in each others k-nearest set
#' @param binary.grf Set all edges >0 to 1
#' @param plots Whether or not to plot the final graph, a heatmap of calculated kernel, and the eigen values of the Laplacian
#' @param verbose Whether or not to give some summary statistics of the pairwise distances
#'
#' @details The four Lapalacians are defined as \eqn{L_{shift}=I+D^{-1/2}AD^{-1/2}}, \eqn{L_{Ng}=D^{-1/2}AD^{-1/2}}, \eqn{L_{sym}=I-D^{-1/2}AD^{-1/2}}, and \eqn{L_{rw}=I-D^{-1}A}. The shifted Laplacian, \eqn{L_{shift}=I+D^{-1/2}AD^{-1/2}}, is recommended for multi-view spectral clustering.
#'
#' @return An n\eqn{\times}n matrix where `n` is the number of rows in `dat`.
#' @references \url{https://academic.oup.com/bioinformatics/article/36/4/1159/5566508#199177546}
#' @examples
#'
#' ## Generating data with 3 distinct clusters
#' ## Note that 'clustStruct' returns a list
#' dd <- clustStruct(n=120, p=30, k=3, noiseDat='random')[[1]]
#'
#' kernelLaplacian(dd, kernel="Spectrum")
#' @export
kernelLaplacian <- function(dat,
                      kernel=c('Gaussian', 'ZM', 'Spectrum', 'Linear'),
                      laplacian = c('shift', 'Ng', 'sym', 'rw'),
                      grf.type = c('full', 'knn', 'e-graph'),
                      k=5,
                      p=5,
                      rho=NULL,
                      epsilon=0,
                      mutual=FALSE,
                      binary.grf=FALSE,
                      plots=TRUE,
                      verbose=TRUE){

  if(missing(kernel)) kernel <- 'Spectrum'
  if(missing(laplacian)) laplacian <- 'shift'
  if(missing(grf.type)) grf.type <- 'full'

  if(is.null(rho)) rho <- stats::median(stats::dist(dat))
  if(kernel == 'Gaussian') A <- Gaussian_kernel(dat, rho = rho)

  if(kernel == 'ZM') A <- Zhang_kernel(dat, p)
  if(kernel == 'Spectrum'){
    dtd <- as.data.frame( t(dat) )
    A <- Spectrum::CNN_kernel(dtd)
  }

  if(kernel == 'Linear') A <- dat%*%t(dat)

  if(verbose){
    message(paste('Distances calculated with', kernel, 'kernel.'))
    message(paste("Returning", laplacian, "Laplacian from a", grf.type, "graph."))
    message("Summary of parwise distances:")
    print(summary(A[lower.tri(A)]))
  }

  Laplacian(A, laplacian = laplacian, grf.type = grf.type,
            k=k, epsilon=epsilon, mutual=mutual,
            binary.grf=binary.grf, plots=plots)
}






