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

getAitchisonVar = function(x){
  compositions::variation(x = compositions::acomp(x))
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

# Transform Samples with the ilr of a Balance
#
#  x A matrix with rows as samples (N) and columns as components (D).
#  contrast A vector. One column of a serial binary partition matrix
#  with values [-1, 0, 1] describing D components.
#
#  A transformation of samples for the balance provided.
balance.fromContrast <- function(x, contrast){

  if(length(contrast) != ncol(x)) stop("Contrast must have length ncol(x) = D.")
  if(any(!contrast %in% c(-1, 0, 1))) stop("Contrast must contain [-1, 0, 1] only.")

  lpos <- sum(contrast == 1)
  lneg <- sum(contrast == -1)
  const <- sqrt((lpos*lneg)/(lpos+lneg))

  logX <- log(x)
  ipos <- rowMeans(logX[, contrast == 1, drop = FALSE])
  ineg <- rowMeans(logX[, contrast == -1, drop = FALSE])

  const * log(exp(ipos) / exp(ineg))
}

# Compute Balances from an SBP Matrix
#
#  x A matrix with rows as samples (N) and columns as components (D).
#  y A serial binary partition matrix with rows as components (D) and
#  columns as balances (D-1).
#
# A transformation of samples for each balance in the SBP matrix.
balance.fromSBP <- function(x, y){

  if(!identical(colnames(x), rownames(y))){

    stop("Component names for data matrix and balance matrix do not match.")
  }

  x <- as.matrix(x)

  if(any(x == 0)){

    message("Alert: Replacing 0s with next smallest value to calculate balances.")
    zeros <- x == 0
    x[zeros] <- min(x[!zeros])
  }

  res <- apply(y, 2, function(z) balance.fromContrast(x, z))
  rownames(res) <- as.character(1:nrow(res))
  return(res)
}

# Make response for cox regression
#
# Internal function to make the response y passed to slr suitable
# for cox regression (i.e. with family = "cox"). Sanity checks are performed
# here too.
#
# If y is a class "Surv" object, this function returns y with no changes. If
# y is a two-column matrix with columns named 'time' and 'status', it is
# converted into a "Surv" object.
#
#  y Response variable. Either a class "Surv" object or a two-column
# matrix with columns named 'time' and 'status'.
#
#  A class "Surv" object.
#
#  survival Surv
response.coxnet <- function(y) {
  if (any(is.na(y))) stop(paste0("NAs encountered in response, not allowed"))

  # if Surv object, check that it is of correct type and perform sanity checks
  # One sanity check is that it have column names. If Surv() is called with a one-column matrix for
  # time, the name is lost
  # if all good, return with no changes
  if (survival::is.Surv(y)) {
    if (attr(y, "type") == "right") {
      if (any(y[, 1] <= 0))
        stop("Non-positive event times encountered; not permitted for Cox family")
      colnames(y) <- c("time","status")
      return(y)
    } else if (attr(y, "type") == "counting") {
      if (any(y[, 1] < 0) || any(y[, 2] <= 0))
        stop(paste("Negative start/non-positive stop times encountered;",
                   "not permitted for Cox family"))
      if (any(y[, 1] >= y[, 2]))
        stop("Some rows have start time >= stop time; not permitted")
      colnames(y) <- c("start","stop","status")
      return(y)
    } else {
      stop("cox.path() only supports 'Surv' objects of type 'right' or 'counting'")
    }
  }

  # if two-column matrix passed, make it into a Surv object
  if (!is.matrix(y) || !all(match(c("time","status"),dimnames(y)[[2]],0)))
    stop(paste0("Cox model requires a matrix with columns 'time' (>0) and ",
                "'status' (binary) as a response; a 'Surv' object suffices"),
         call. = FALSE)
  ty <- as.double(y[,"time"])
  tevent <- as.double(y[,"status"])
  if (any(ty <= 0))
    stop("negative event times encountered; not permitted for Cox family")
  yob <- survival::Surv(ty, tevent)
  colnames(yob) <- c("time","status")
  return(yob)
}

# get log contrast coefficients from a regression coefficient in the balance regression
balReg2llc = function(theta, sbp, normalized = TRUE){
  # theta defined in manuscript, i.e., coefficient of normalized balance
  if(!is.matrix(sbp)) sbp = matrix(sbp)
  if(ncol(sbp) != length(theta)){
    stop("SBP and coefficients do not match!")
  }
  kplus = apply(sbp, 2, function(col) sum(col == 1))
  kminus = apply(sbp, 2, function(col) sum(col == -1))
  const <- sqrt( (kplus * kminus)/(kplus + kminus) )

  reciprocals = matrix(
    0, nrow = nrow(sbp), ncol = ncol(sbp))
  for(i in 1:ncol(sbp)){
    if (normalized){
      reciprocals[sbp[, i] == 1, i] = const[i] / kplus[i]
      reciprocals[sbp[, i] == -1, i] = - const[i] / kminus[i]
    } else {
      reciprocals[sbp[, i] == 1, i] = 1 / kplus[i]
      reciprocals[sbp[, i] == -1, i] = -1 / kminus[i]
    }
  }

  ReciprocalstimesCoeffs = matrix(
    NA, nrow = nrow(sbp), ncol = ncol(sbp))
  for(i in 1:ncol(ReciprocalstimesCoeffs)){
    ReciprocalstimesCoeffs[, i] = reciprocals[, i] * theta[i]
  }
  beta = rowSums(ReciprocalstimesCoeffs)
  names(beta) = rownames(sbp)
  return(beta)
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


#' Microbiome composition related to ulcerative colitis
#'
#' @name UC
#' @docType data
#' @description
#'
#' A processed microbiome data set consisting of the counts of 447 different genera in a group of 132 subjects.
#' There are two independent cohorts, which is given by the column \code{cohort}.
#' There are 34 controls and 53 UC subjects in the discovery cohort (\code{cohort=='PRISM'}),
#' and 22 controls and 23 UC subjects in the validation cohort (\code{cohort=='Validation'}).
#' For details on data preprocessing, please see documentation on GitHub \url{https://github.com/drjingma/LogRatioReg/tree/master}.
#'
#' @usage data(UC)
#' @format This data set is organized as a \code{data.frame} with the first 447 columns being genera counts, the next
#' column indicating the sample's status (UC vs Control), and the last column indicating the cohort.
#'
#' @author Jing Ma (\email{jingma@fredhutch.org})
#'
#' @references
#'
#' \url{https://doi.org/10.1038/s41564-018-0306-4}
#'
#' @keywords data
NULL
