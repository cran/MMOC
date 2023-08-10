###############################
##
## Project: MMOC
##
## Purpose: Helper functions for the package
##
## Author: Charlie Carpenter
## Email: charles.carpenter@cuanschutz.edu
##
## Date Created: 2022-08-16
##
## ---------------------------
## Notes:
##
##
## ---------------------------


# Kernels -----------------------------------------------------------------

Gaussian_kernel <- function(Z, rho){
  exp(-(1/rho)*as.matrix(stats::dist(Z, method = "euclidean", upper = TRUE)^2))
}

Zhang_kernel <- function(Z, p){
  d <- as.matrix( stats::dist(Z) )

  si <- apply(d, 2, function(s) sort(s)[p+1] )
  db <- si %o% si
  dt <- -d^2/db

  exp(dt)
}


# knn Graph ---------------------------------------------------------------

knnGrf <- function(A, k, mutual=FALSE){
  knn <- function(aa, k){
    aa[order(aa)[(k+2):length(aa)]] <- 0
    aa[aa > 0] <- 1
    aa
  }

  akn <- apply(A,2,knn,k=k)

  if(mutual){
    for(i in 1:nrow(akn)){
      for(j in 1:ncol(akn)){
        akn[i,j] <- ifelse(akn[i,j]==akn[j,i],
                           akn[i,j],
                           0)
      }
    }
  } else{
    for(i in 1:nrow(akn)){
      for(j in 1:ncol(akn)){
        akn[i,j] <- ifelse(akn[i,j]!=akn[j,i],
                           max(akn[i,j],akn[j,i]),
                           akn[i,j])
      }
    }
  }
  akn
}


# Basic Matrix Ops --------------------------------------------------------

## Trace of a matrix
tr <- function(x){
  if(!is.matrix(x)) stop("x must be a matrix")
  sum(diag(x))
}

## Check if matrix is orthogonal
orthoCheck <- function(x){
  xtx <- crossprod(x)
  I <- diag(nrow = nrow(xtx))

  # Same tolerance as 'isSymmetric'
  sum((I-xtx)^2) < 100*.Machine$double.eps
}


# Data Generation ---------------------------------------------------------

# clustStruct <- function(n, p, k, noiseDat=NULL, randNoise=2){
#   if(any(n%%k!=0)) stop("n must be divisible by k.")
#
#   lapply(k, function(kk){
#
#     means <- c(0, 2^( 1:(kk-1) ) )
#     datL <- lapply(means, function(mm){
#       MASS::mvrnorm(n/kk, rep(mm,p), diag(p))
#     })
#
#     dat <- do.call(rbind, datL)
#
#     if(!is.null(noiseDat)){
#       if(is.character(noiseDat)){
#         S <- randNoise*diag(p)
#         noiseDat <- MASS::mvrnorm(n=n, mu=rep(0,p), Sigma=S)
#       }
#       dat <- dat+noiseDat
#     }
#
#     dat
#   })
#
# }

getGroups <- function(n,k){
  ttt <- numeric(0)
  for(i in 1:length(k)){
    tk <- n/k[i]
    ttt <- cbind(ttt, rep(1:k[i], each=tk))
  }
  as.data.frame(ttt)
}


