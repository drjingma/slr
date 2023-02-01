alrinv <- function(y) {
  x <- cbind(exp(y),1)
  x / rowSums(x)
}

sigmoid = function(x){
  1 / (1 + exp(-x))
}

graph.laplacian <- function(W, normalized = TRUE, zeta=0.01){
  stopifnot(nrow(W) == ncol(W))

  n = nrow(W)    # number of vertices
  # We perturb the network by adding some links with low edge weights
  W <- W + zeta * mean(colSums(W))/n * tcrossprod(rep(1,n))
  g <- colSums(W) # degrees of vertices

  if(normalized){
    D_half = diag(1 / sqrt(g) )
    return(D_half %*% W %*% D_half )
  } else {
    return(W)
  }
}

#'
#' @description Performs spectral clustering for a given similarity matrix and pre-specified number of clusters.
#'
#' @param W The similarity matrix.
#' @param k Number of clusters.
#' @param zeta A small positive number used to perturb the similarity matrix by adding some links with low edge weights in the calculation of the normalized Laplacian.
#'
#' @return cluster assignments
#'
#' @details For spectral clustering, it can be beneficial to perturb the similarity matrix with a small positive value to improve the clustering performance. This leads to regularized spectral clustering.
#' For more details on regularized spectral clustering of networks, see Amini et al. (13').
#'
#' @author Jing Ma and Kristyn Pantoja.
#'
#' Maintainer: Jing Ma (\url{jingma@fredhutch.org})
#'
#' @references
#' Amini, A. A., Chen, A., Bickel, P. J., & Levina, E. (2013). Pseudo-likelihood methods for community detection in large sparse networks. Annals of Statistics. 41(4): 2097-2122
#'
spectral.clust <- function(W, k, zeta = 0) {
  L = graph.laplacian(W,zeta = zeta) # Compute graph Laplacian
  ei = eigen(L, symmetric = TRUE)    # Compute the eigenvectors and eigenvalues of L
  # we will use k-means to cluster the eigenvectors corresponding to the leading smallest eigenvalues
  ei$vectors <- ei$vectors[,base::order(abs(ei$values),decreasing=TRUE)]
  obj <- stats::kmeans(
    ei$vectors[, 1:k], centers = k, nstart = 100, algorithm = "Lloyd")
  if (k==2){
    cl <- 2*(obj$cluster - 1) - 1
  } else {
    cl <- obj$cluster
  }
  names(cl) <- rownames(W)
  # return the cluster membership
  return(cl)
}

AitchVar = function(x, y){
  stats::var(log(x) - log(y))
}

AitchVarVec = Vectorize(AitchVar)

getAitchisonVar = function(x){
  outer(X = x, Y = x, FUN = AitchVarVec)
}

sbp.fromHclust <- function(hclust){

  if(class(hclust) != "hclust"){

    stop("This function expects an 'hclust' object.")
  }

  labels <- hclust$labels
  merge <- hclust$merge

  out <- matrix(0, nrow(merge), nrow(merge) + 1)
  if(is.null(labels)){
    colnames(out) <- 1:ncol(out)
  }else{
    colnames(out) <- labels
  }

  branches <- vector("list", nrow(merge) - 1)
  for(i in 1:nrow(merge)){

    # Assign +1 to branch 1
    branch1 <- merge[i,1]
    if(branch1 < 0){
      include1 <- -1 * branch1
    }else{
      include1 <- branches[[branch1]]
    }
    out[i, include1] <- 1

    # Assign -1 to branch 2
    branch2 <- merge[i,2]
    if(branch2 < 0){
      include2 <- -1 * branch2
    }else{
      include2 <- branches[[branch2]]
    }
    out[i, include2] <- -1

    # Track base of branch
    branches[[i]] <- c(include1, include2)
  }

  # Sort balances by tree height
  sbp <- t(out[nrow(out):1,])
  colnames(sbp) <- paste0("z", 1:ncol(sbp))
  sbp
}


# Get a contrast vector from a given binary partition vector
#' @description Fill zeros to a given binary partition vector to obtain a contrast.
#'
#' @param bp Binary partition returned by the \code{slr} function.
#' @param X A sample by variable matrix of strictly positive relative abundances (\code{n} by \code{p}).
#'
#' @return A contrast vector of length \code{ncol(X)}.
#'
getContrast <- function(bp, X){
  contrast = matrix(0, nrow = ncol(X), ncol = 1)
  rownames(contrast) = colnames(X)
  contrast[match(names(bp), rownames(contrast))] = bp
  contrast
}

# Balance from a contrast vector
#' @param X A sample by variable matrix of strictly positive relative abundances (\code{n} by \code{p}). This matrix should not contain any zeros. See details.
#' @param contrast Contrast vector that takes values from -1, 0, 1. This vector should be have length \code{ncol(X)}.
#'
#' @return The balance predictor corresponding to \code{contrast}.
balance.fromContrast <- function(X, contrast){

  if(length(contrast) != ncol(X)) stop("Contrast must have length ncol(x).")
  if(any(!contrast %in% c(-1, 0, 1))) stop("Contrast must contain [-1, 0, 1] only.")

  logX <- log(X)
  ipos <- rowMeans(logX[, contrast == 1, drop = FALSE])
  ineg <- rowMeans(logX[, contrast == -1, drop = FALSE])

  ipos - ineg
}
