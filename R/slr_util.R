getGammaFromTheta = function(theta, sbp){
  # theta defined in manuscript, i.e., coefficient of non-normalized balance
  if(!is.matrix(sbp)) sbp = matrix(sbp)
  if(ncol(sbp) != length(theta)){
    stop("getGammaFromTheta: SBP and coefficients don't match")
  }
  kplus = apply(sbp, 2, function(col) sum(col == 1))
  kminus = apply(sbp, 2, function(col) sum(col == -1))
  reciprocals = matrix(
    0, nrow = nrow(sbp), ncol = ncol(sbp))
  for(i in 1:ncol(sbp)){
    reciprocals[sbp[, i] == 1, i] = 1 / kplus[i]
    reciprocals[sbp[, i] == -1, i] = -1 / kminus[i]
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

#' Get the coefficients of a balance regression model in terms of the linear log contrast model
#'
#' @param coefs coefficients of balance regression model (may include intercept)
#' @param sbp fully-specified SBP vector/matrix
#'
#' @return the coefficients of a balance regression model in terms of the linear log contrast model
#' @export
#'
#' @examples
getCoefsBM = function(coefs, sbp){
  # balance model defined in manuscript, i.e., without normalizing balance
  if(names(coefs)[1] == "(Intercept)"){ # model with intercept
    a0 = coefs[1]
    if(length(coefs) == 1){
      bm.coefs = 0
    } else{
      bm.coefs = coefs[-1]
    }
  } else{ # model without intercept
    a0 = NA
    bm.coefs = coefs
  }
  llc.coefs = getGammaFromTheta(bm.coefs, sbp)
  return(list(
    a0 = a0,
    bm.coefs = bm.coefs,
    llc.coefs = llc.coefs
  ))
}

#' Get a fully-specified SBP vector for an slr model
#'
#' @param slr.sbp sbp returned by slr object
#' @param x data matrix
#'
#' @return a fully-specified SBP vector
#' @export
#'
#' @examples
getSlrSBP = function(slr.sbp, x){
  p = ncol(x)
  slrspec.fullSBP = matrix(0, nrow = p, ncol = 1)
  rownames(slrspec.fullSBP) = colnames(x)
  slrspec.fullSBP[match(
    names(slr.sbp), rownames(slrspec.fullSBP))] = slr.sbp
}
