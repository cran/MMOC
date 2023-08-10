#' @title Generate multi-view data sets with simple cluster structures
#' @name clustStruct
#' @description Generates multiple data sets from a multivariate normal distribution using the \link[MASS]{mvrnorm} function from the `MASS` package.
#'
#' @param n An integer, the sample size for all generated data sets
#' @param p An integer, the number of columns (features) in each generated data set
#' @param k An integer or vector, the number of distinct clusters in each generated data set. `n/k` must be an integer for all values of `k`
#' @param noiseDat Either the character string `'random'`, indicating the covariance matrix is a diagonal matrix with `randNoise` along the diagonal, or a valid covariance matrix
#' @param randNoise The value along the diagonal when `noiseDat='random'`
#'
#' @details The function accepts `k` as a vector. It splits data into `k` groups with means `c(0, 2^( 1:(kk-1) ) )`, e.g., when `k=3` the data will be split into 3 groups with means 0, 2, and 4, respectively. The covariance matrix is either a diagonal matrix with `randNoise` (an integer) along the diagonal, or a given matrix.
#' @returns A list of n\eqn{\times}p data frames with the specified number of groups
#' @examples
#'
#' ## A single view with 30 variables and 3 groups
#' s1 <- clustStruct(n=120, p=30, k=3, noiseDat='random')[[1]]
#'
#' ## Multiple views with 30 variables
#' ## View 1 has 2 groups and View 2 has 3 groups
#' s2 <- clustStruct(n=120, p=30, k=c(2,3), noiseDat='random')
#'
#' ## Multiple views with 30 variables
#' ## View 1 has 2 groups, View 2 has 3, and View 3 has 3 groups
#' s3 <- clustStruct(n=120, p=30, k=c(2,3,3), noiseDat='random')
#'
#' ## Three view study.
#' # View 1: 2 groups, 30 variables, random noise = 5
#' # View 2: 3 groups, 60 variables, random noise = 2
#' # View 3: 4 groups, 45 variables, random noise = 4
#'
#' s4 <- clustStruct(n=120, k=c(2,3,4), p=c(30,60,45), randNoise=c(5,2,4))
#'
#' @export
clustStruct <- function(n, p, k, noiseDat='random', randNoise=2){
  if(any(n%%k!=0)) stop("n must be divisible by k.")
  if(any(k <= 0)) stop("k must be all positive integers")
  if(any(as.integer(k)!=k)) stop("k must be all positive integers")

  mapply(function(kk, pp, rn){

    means <- c(0, 2^( 1:(kk-1) ) )
    datL <- lapply(means, function(mm){
      MASS::mvrnorm(n/kk, rep(mm,pp), diag(pp))
    })

    dat <- do.call(rbind, datL)

    if(!is.null(noiseDat)){
      if(is.character(noiseDat)){
        S <- rn*diag(pp)
        noiseDat <- MASS::mvrnorm(n=n, mu=rep(0,pp), Sigma=S)
      }

      if(!all(dim(dat) == dim(noiseDat))) stop("'noiseDat' must have same 'n' and 'p' as simulated data")

      dat <- dat+noiseDat
    }

    dat
  }, kk=k, pp=p, rn=randNoise, SIMPLIFY=FALSE)

}







