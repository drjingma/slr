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

spectral.clust <- function(W, k, zeta = 0) {
  L = graph.laplacian(W,zeta = zeta) # Compute graph Laplacian
  ei = eigen(L, symmetric = TRUE)    # Compute the eigenvectors and eigenvalues of L
  # we will use k-means to cluster the eigenvectors corresponding to the largest eigenvalues in absolute value
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

  if(!inherits(hclust, "hclust")){
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

balance.fromBP <- function(X, bp){

  if(length(bp) > ncol(X)) stop("bp must have length no greater than ncol(x).")
  if(any(!bp %in% c(-1, 1))) stop("bp must contain [-1, 1] only.")
  if(length(bp) < 2) stop("bp must have length no smaller than 2.")

  logX <- log(X[, match(names(bp),colnames(X))])
  ipos <- rowMeans(logX[, bp == 1, drop = FALSE])
  ineg <- rowMeans(logX[, bp == -1, drop = FALSE])

  ipos - ineg
}

#' Load the HIV data set
#' @description Load the HIV data set from the selbal R package.
#'
#' @details This function is a wrapper to import the \code{HIV} data set from the selbal R package.
#' Please make sure the package is installed by visiting \url{https://github.com/malucalle/selbal}.
#'
#' \code{HIV} is a data frame with 155 rows (samples) and 62 columns (variables).
#' The first 60 variables measure the counts of bacterial species at the genus taxonomy rank.
#' The last column \code{HIV_Status} is a factor indicating the HIV infection status:
#' Pos/Neg if an individual is HIV1 positive/negative. There is also a column \code{MSM}
#' which is an HIV risk factor, Men who has sex with men (MSM) or not (nonMSM).
#'
#' @references
#'
#' \url{https://pubmed.ncbi.nlm.nih.gov/27077120/}
#'
#' @export
#'
load_data <- function() {
  # check if package is installed
  if (requireNamespace("selbal", quietly = TRUE)) {
    # get name of random dataset
    x <- utils::data(list = "HIV", package = "selbal", envir = environment())
    return(get(x))
  } else {
    stop("Install package from https://github.com/malucalle/selbal first.")
  }
}
