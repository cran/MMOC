#' @title Calculate the Flag mean of multiple subspaces
#' @name flagMean
#' @description Calculate the flag-mean of multiple subspaces. This method allows you to find the extrinsic mean of a finite set of subspaces. You can think of this as a median subspace. This method is also able to handle subspaces with different dimensions. See the references for more details
#'
#' @param LapList A list of Laplacian matrices
#' @param k A vector indicating how many eigenvectors to take from each Laplacian, i.e., the number of clusters in each view
#' @param laplacian One of `"shift"`, `"Ng"`, `"rw"` or `"sym"`. Should be the same type used to calculate your Laplacians
#' @param plots Whether or not to plot the singular values from SVD
#'
#' @return The output from a singular value decomposition. See \link[base]{svd}
#' @details Despite the complex linear algebra to achieve this result, the opperation is very simple. This function concatonates (cbind) the given subspaces and then performs singular value decomposition on the resulting matrix. This gives the 'median' subspace of the given set of subspaces. We would then cluster on the columns of the `U` matrix just as we do in standard spectral clustering
#'
#' @references
#' https://www.semanticscholar.org/paper/Flag-Manifolds-for-the-Characterization-of-in-Large-Marrinan-Beveridge/7d306512b545df98243f87cb8173df83b4672b18
#' https://www.sciencedirect.com/science/article/pii/S0024379514001669?via%3Dihub
#'
#' @examples
#'
#' ## Generating data with 2 and 3 distinct clusters
#' ## Note that 'clustStruct' returns a list
#' n=120; k <- c(2,3)
#' dd <- clustStruct(n=n, p=30, k=k, noiseDat='random')
#'
#' ## Laplacians
#' L_list <- lapply(dd, kernelLaplacian, kernel="Spectrum", plots=FALSE, verbose=FALSE)
#'
#' ## Calculating the flag mean
#' fm <- flagMean(L_list, k=k, laplacian='shift')
#'
#' ## Knowing the true structure makes it much easier to know how
#' ## many right singular vectors to grab. There are 4 distinct
#' ## groups in these data from 'clustStruct'
#'
#' trueGroups(n=n, k=k)
#'
#' kmeans(fm$u[, 1:4], centers=4)
#' @export
flagMean <- function(LapList, k,
                     laplacian = c('shift', 'Ng', 'sym', 'rw'),
                     plots=TRUE){

  if(missing(laplacian)) laplacian <- 'shift'

  ## Eigen decomposition of Laplacians
  EigList <- lapply(LapList, eigen)
  Ul <- mapply(smlEig, LapList, k,
               lap.type=laplacian, SIMPLIFY = FALSE)
  svdX <- FM(Ul)

  if(plots) plot(svdX$d, col='dodgerblue3', pch=20, type="b",
                 ylab="Singular Values", main="Singular values of the flag mean")

  message(paste("Returning flag mean from", laplacian, "Laplacians"))

  svdX
}

FM <- function(Ul){
  ## Check orthogonal
  # if(!all(sapply(Ul, orthoCheck))) stop("All matrices must have orthogonal columns")

  r <- max(vapply(Ul, ncol, 0))

  ## Concatenating matrices
  X <- do.call(cbind, Ul)
  svd(X)
}

smlEig <- function(L, k, lap.type){
  n <- nrow(L)

  if(lap.type %in% c('sym', 'rw')){
    ind <- (n-k+1):n
  } else ind <- 1:k

  eigen(L)$vectors[,ind, drop=FALSE]
}
