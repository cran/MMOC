## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(MMOC)

## For plotting
library(plotly)

## -----------------------------------------------------------------------------
n <- 120 # Sample size
k <- c(2,3,4) # Number of clusters
p <- c(60, 90, 120) # Number of features
rn <- 4 # Amount of extra noise added to each cluster

dd <- clustStruct(n=n, p=p, k=k)

## -----------------------------------------------------------------------------
trueGroups(n,k)

## -----------------------------------------------------------------------------
## shifted Laplacians from a Spectrum Kernel on first data set
kl1 <- kernelLaplacian(dd[[1]], laplacian='shift')

## -----------------------------------------------------------------------------
## shifted Laplacians from a Spectrum Kernel on every data set
LapList <- lapply(dd, kernelLaplacian, laplacian='shift', plots=FALSE, verbose=FALSE)

## -----------------------------------------------------------------------------
## Calculating the approximate Laplacain
LR <- Lapprox(LapList, k=k, laplacian='shift')

## Clustering on the first K eigenvectors (K=6)
set.seed(234)
km <- kmeans(LR$vectors[,1:6], centers = 6, nstart = 25)
clusters <- factor(km$cluster)

## Plotting our estimated clusters in the first 3 eigenvectors
LR$vectors %>% as.data.frame %>% 
  mutate(clusters=clusters) %>%
  plot_ly() %>% 
  add_markers(x=~V1, y=~V2, z=~V3,
              color=~clusters, colors='Set1')

## -----------------------------------------------------------------------------
fm <- flagMean(LapList, k=k, laplacian = 'shift')

set.seed(234)
km <- kmeans(fm$u[,1:6], centers = 6, nstart = 25)
clusters <- factor(km$cluster)

## Plotting our estimated clusters in the first 3 left singular vectors
fm$u %>% as.data.frame %>% 
  mutate(clusters=clusters) %>%
  plot_ly() %>% 
  add_markers(x=~V1, y=~V2, z=~V3,
              color=~clusters, colors='Set1')

## ---- eval=FALSE--------------------------------------------------------------
#  LRk <- NbClust(LR$vectors[,1:sum(k)], method = 'keans',
#                 index='alllong')
#  
#  fmk <- NbClust(fm$u, method = 'keans', index='alllong')

## -----------------------------------------------------------------------------
library(Spectrum)

dt <- lapply(dd, t)
dt <- lapply(dt, function(x){
  rownames(x) <- 1:nrow(x)
  colnames(x) <- 1:ncol(x)
  x
})

spec <- Spectrum(dt, method=2)

fm$u %>% as.data.frame %>% 
  mutate(clusters=factor(spec$assignments)) %>%
  plot_ly() %>% 
  add_markers(x=~V1, y=~V2, z=~V3,
              color=~clusters, colors='Set1')

## -----------------------------------------------------------------------------
library(SNFtool)
K <- 20;		# number of neighbors, usually (10~30)
alpha <- 0.5;  	# hyperparameter, usually (0.3~0.8)
iters <- 20; 	# Number of Iterations, usually (10~20)

Wl <- lapply(dd, function(x){
  d <- (dist2(as.matrix(x),as.matrix(x)))^(1/2)
  affinityMatrix(d, K, alpha)
})

W <- SNF(Wl, K, iters)

## Telling SNF the true number of groups
snfGroup <- spectralClustering(W, K=6)

fm$u %>% as.data.frame %>% 
  mutate(clusters=factor(snfGroup)) %>%
  plot_ly() %>% 
  add_markers(x=~V1, y=~V2, z=~V3,
              color=~clusters, colors='Set1')


