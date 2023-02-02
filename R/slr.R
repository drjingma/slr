getFeatureScores = function(X, y, screen.method, response.type, s0.perc){
  n = length(y)

  ## Compute univariate coefficient
  xclr <- apply(X, 1, function(a) log(a) - mean(log(a)))
  if (screen.method=='wald'){
    xclr.centered <- base::scale(t(xclr),center=TRUE,scale=TRUE)
    if (response.type=='continuous'){
      y.centered <- y - mean(y)
      sxx <- diag(crossprod(xclr.centered))
      sxy <- crossprod(xclr.centered,y.centered)
      syy <- sum(y.centered^2)
      numer <- sxy/sxx
      sd <- sqrt((syy/sxx - numer^2)/(n - 2)) # variance estimate

      if (is.null(s0.perc)) {
        fudge <- stats::median(sd)
      }
      if (!is.null(s0.perc)) {
        if (s0.perc >= 0) {
          fudge <- stats::quantile(sd, s0.perc)
        }
        if (s0.perc < 0) {
          fudge <- 0
        }
      }
      # this is the wald-statistic
      fs <- numer/(sd + fudge)

      # t-distribution with n-2 df
      fs <- stats::pt(abs(fs), df = n-2)
    } else if (response.type=='binary'){
      fs <- rep(0,ncol(xclr.centered))
      for (j in 1:ncol(xclr.centered)){
        fit <- stats::glm(y~x,data=data.frame(
          x=xclr.centered[,j],y=as.factor(y)),family=stats::binomial(link='logit'))
        fs[j] <- stats::coef(summary(fit))[2,3]
      }
      fs <- stats::pnorm(abs(fs))
    }
  }
  if (screen.method=='correlation'){
    fs <- stats::cor(t(xclr),y)
  }
  return(fs)
}


#' Supervised Log Ratio
#'
#' @param X A sample by variable matrix of strictly positive relative abundances (\code{n} by \code{p}). This matrix should not contain any zeros. See details.
#' @param y The response vector of length \code{n}.
#' @param screen.method Method of variable screening: could be correlation based ("correlation") or wald score based ("wald"). Default is \code{wald}.
#' @param cluster.method Method of clustering: "spectral" or "hierarchical". Default is \code{spectral}.
#' @param response.type Type of the response variable: could be "continuous" or "binary".
#' @param threshold A nonnegative constant between 0 and 1. If \code{NULL}, then no variable screening is performed.
#' @param s0.perc The percentile of standard deviation values added to the denominator of the wald score statistic. Default is 0. See details.
#' @param zeta A small positive value used to perturb the Aitchison similarity matrix when spectral clustering is used. Default is 0. See details.
#' @param positive.slope Logical flag indicating whether to define the balance such that its corresponding slope estimate is positive. Default is \code{TRUE}.
#'
#' @description \code{slr} fits a balance regression model in which the balance is defined by a sparse set of input variables.
#'
#' @details \code{slr} first uses a screening procedure to identify the active variables correlated with the response \code{y}.
#' Essentially, it computes the univariate regression coefficients for centered log ratio transformed \code{X} and
#' forms a reduced data matrix with variables whose univariate coefficients exceed a \code{threshold} in absolute value
#' (the threshold is chosen via cross-validation). Then, it performs clustering of the active variables on a suitable dissimilarity
#' derived from the reduced data matrix to get 2 clusters. Finally, \code{slr} uses the resulting balance to predict \code{y}.
#'
#' The design matrix \code{X} should not contain any zeros. In general, zeros in raw data should be imputed using a pseudocount or Bayesian method prior to applying \code{slr}.
#'
#' When the wald score statistic is used,
#'
#' When spectral clustering is used, it can be beneficial to perturb the Aitchison similarity matrix with a small positive value to improve the clustering performance. This leads to regularized spectral clustering.
#' For more details on regularized spectral clustering of networks, see Amini et al. (13').
#'
#' @return An object of class \code{"slr"} is returned, which is a list with
#'
#' \item{bp}{The binary partition of selected variables. A positive value of 1 indicates the variable is in the numerator, while -1 indicates a denominator variable.}
#' \item{feature.scores}{The feature scores for all variables.}
#' \item{glm.fit}{The generalized linear model fit from a univariate balance regression. }
#' @export
#'
#' @author Jing Ma and Kristyn Pantoja.
#'
#' Maintainer: Jing Ma (\url{jingma@fredhutch.org})
#'
#' @seealso \code{cv.slr}
#'
#' @references
#' Amini, A. A., Chen, A., Bickel, P. J., & Levina, E. (2013). Pseudo-likelihood methods for community detection in large sparse networks. Annals of Statistics. 41(4): 2097-2122
#'
#' @examples
#'
#' HIV <- load_data() # Load HIV data
#' X <- HIV[,1:60]
#' y <- ifelse(HIV[,62] == "Pos", 1, 0)
#' X.adjusted <- sweep(X+1,rowSums(X+1),MARGIN = 1, FUN='/') # zero handling
#'
#' # Run slr ----
#' fit <- slr(X.adjusted, y, screen.method='wald', cluster.method ='spectral',
#'            response.type = 'binary', threshold = 0.9, positive.slope = TRUE)
#' fit$bp
slr = function(
    X,
    y,
    screen.method=c('correlation','wald'),
    cluster.method = c('spectral', 'hierarchical'),
    response.type=c('survival','continuous','binary'),
    threshold,
    s0.perc=0,
    zeta=0,
    positive.slope = FALSE
){
  this.call <- match.call()
  screen.method <- match.arg(screen.method)
  response.type <- match.arg(response.type)
  if(!("data.frame" %in% class(X))) X = data.frame(X)

  n <- length(y)

  feature.scores = getFeatureScores(X, y, screen.method, response.type, s0.perc)
  which.features <- (abs(feature.scores) >= threshold)
  if (sum(which.features)<2){
    # Fit an intercept only regression model
    if (response.type=='binary'){
      model.train <- stats::glm(y~.,data=data.frame(
        y=as.factor(y)),family=stats::binomial(link='logit'))
    } else if (response.type=='continuous'){
      model.train <- stats::lm(y~.,data=data.frame(y=y))
    }
    object <- list(bp=NULL, Aitchison.var = NULL, cluster.mat = NULL)
  } else {
    x.reduced <- X[,which.features] # reduced data matrix
    Aitchison.var = getAitchisonVar(x.reduced)
    rownames(Aitchison.var) <- colnames(Aitchison.var) <- colnames(x.reduced)
    if(cluster.method == "spectral" | nrow(Aitchison.var) == 2){
      Aitchison.sim <- max(Aitchison.var) - Aitchison.var
      ## Perform spectral clustering
      bp.est <- spectral.clust(Aitchison.sim, k=2, zeta = zeta)
      cluster.mat = Aitchison.sim
    } else if(cluster.method == "hierarchical"){
      ## Perform hierarchical clustering
      htree.est <- stats::hclust(stats::dist(Aitchison.var))
      bp.est <- sbp.fromHclust(htree.est)[, 1] # grab 1st partition
      cluster.mat = Aitchison.var
    } else{
      stop("invalid cluster.method arg was provided!!")
    }
    balance <- balance.fromBP(x.reduced, bp.est) # predict from labeled data
    # model fitting
    if (response.type=='binary'){
      model.train <- stats::glm(y~balance,data=data.frame(balance=balance,y=as.factor(y)),
                                family=stats::binomial(link='logit'))
      if(positive.slope){
        if(stats::coef(model.train)[2] < 0){
          bp.est = - bp.est
          balance <- balance.fromBP(x.reduced, bp.est)
          model.train <- stats::glm(y~balance,data=data.frame(balance=balance,y=as.factor(y)),
                                    family=stats::binomial(link='logit'))
        }
      }
    } else if (response.type=='continuous'){
      model.train <- stats::lm(y~balance,data=data.frame(balance=balance,y=y))
      if(positive.slope){
        if(stats::coef(model.train)[2] < 0){
          bp.est = - bp.est
          balance <- balance.fromBP(x.reduced, bp.est)
          model.train <- stats::lm(y~balance,data=data.frame(balance=balance,y=y))
        }
      }
    }
    object <- list(bp = bp.est)
  }
  object$feature.scores <- feature.scores
  object$glm.fit <- model.train

  class(object) <- 'slr'
  return(object)
}

slr.predict <- function(
    object, newdata = NULL, response.type=c('continuous','binary')
){
  # prediction will be based on the canonical space
  if (missing(newdata) || is.null(newdata)) {
    stop('No new data provided!')
  } else {
    if (is.null(object$bp)){
      if (response.type=='binary'){
        predictor <- sigmoid(rep(1,nrow(newdata)) * as.numeric(stats::coef(object$glm.fit)))
      } else {
        predictor <- rep(1,nrow(newdata)) * as.numeric(stats::coef(object$glm.fit))
      }
    } else {
      newdata.reduced <- newdata[,colnames(newdata) %in% names(object$bp)]
      new.balance <- balance.fromBP(newdata.reduced,object$bp)
      if (response.type=='binary'){
        fitted.results <- stats::predict(
          object$glm.fit,newdata=data.frame(balance=new.balance),type='response')
        predictor = fitted.results
      } else if (response.type=='continuous'){
        predictor <- cbind(1,new.balance) %*% as.numeric(stats::coef(object$glm.fit))
      }
    }
    as.numeric(predictor)
  }
}

buildPredmat <- function(
    outlist,threshold,X,y,foldid,response.type,type.measure
){
  nfolds = max(foldid)
  predmat = matrix(NA, nfolds, length(threshold))
  for(i in 1:nfolds){ # predict for each fold
    which = foldid == i
    y.i = y[which]
    fitobj = outlist[[i]]
    x.i = X[which, , drop=FALSE]
    predy.i = sapply(
      fitobj, function(a) slr.predict(
        a,newdata=x.i,response.type=response.type))
    for(j in 1:length(threshold)){
      predy.ij = predy.i[, j]
      if (response.type == 'continuous'){ # mse, to be minimized
        if(type.measure != "mse"){
          stop("if response.type is continuous, then type.measure must be mse!!")
        }
        predmat[i, j] <- mean((as.numeric(y.i)-predy.ij)^2)
      } else if (response.type=='binary'){
        if(!(type.measure %in% c("accuracy", "auc"))){
          stop("if response.type is binary, then type.measure must be either accuracy or auc!!")
        }
        if(type.measure == "accuracy"){# accuracy, minimize the # that don't match
          predmat[i, j] <- mean((predy.ij > 0.5) != y.i)
        } else if(type.measure == "auc"){# auc, minimize 1 - auc
          predmat[i, j] = tryCatch({
            1 - pROC::auc(y.i,predy.ij, levels = c(0, 1), direction = "<", quiet = TRUE)
          }, error = function(e){return(NA)}
          )
        }
      }
    }
  }
  predmat
}

getOptcv <- function(threshold, cvm, cvsd){
  cvmin = min(cvm, na.rm = TRUE)
  idmin = cvm <= cvmin
  threshold.min = max(threshold[idmin], na.rm = TRUE)
  idmin = match(threshold.min, threshold)
  semin = (cvm + cvsd)[idmin]
  id1se = cvm <= semin
  threshold.1se = max(threshold[id1se], na.rm = TRUE)
  id1se = match(threshold.1se, threshold)
  index=matrix(c(idmin,id1se),2,1,dimnames=list(c("min","1se"),"threshold"))
  list(
    threshold.min = threshold.min,
    threshold.1se = threshold.1se,
    index = index
  )
}


#' Cross Validation for Supervised Log Ratio
#'
#' @param X A sample by variable matrix of strictly positive relative abundances (\code{n} by \code{p}). This matrix should not contain any zeros.
#' @param y The response vector of length \code{n}.
#' @param screen.method Method of variable screening: could be correlation based ("correlation") or wald score based ("wald"). Default is \code{wald}.
#' @param cluster.method Method of clustering: "spectral" or "hierarchical". Default is \code{spectral}.
#' @param response.type Type of the response variable: could be "continuous" or "binary".
#' @param threshold Optional user-supplied threshold sequence; default is \code{NULL} and \code{cv.slr} chooses its own sequence, which is recommended.
#' @param s0.perc The percentile of standard deviation values added to the denominator of the wald score statistic. Default is 0.
#' @param zeta A small positive value used to perturb the Aitchison similarity matrix when spectral clustering is used. Default is 0.
#' @param type.measure Loss used for cross-validation. If \code{response.type='continuous'}, then \code{type.measure="mse"}. For two-class logistic regression, \code{type.measure="auc"}.
#' @param nfolds Number of folds. Default is 10.
#' @param foldid An optional vector of values between 1 and \code{nfold} identifying which fold each observation is in. If supplied, \code{nfold} can be missing.
#' @param weights Observation weights. Default is 1 per observation.
#' @param trace.it If \code{trace.it=TRUE}, then progress bars are displayed.
#' @param plot If \code{TRUE}, then a visualization of the cross-validation results is displayed.
#'
#' @return An object of class \code{"cv.slr"} is returned, which is a list with the ingredients of cross-validation fit.
#' \item{threshold}{A vector of threshold values used. The range of threshold values depends on the marginal association between each variable and the response.}
#' \item{cvm}{The mean cross-validated error - a vector of length \code{length(threshold)}.}
#' \item{cvsd}{The estimate of standard error of \code{cvm}.}
#' \item{foldid}{The fold assignments used.}
#' \item{threshold.min}{Value of \code{threshold} that gives minimum \code{cvm}.}
#' \item{threshold.1se}{Largest value of \code{threshoold} such that the cross-validation error is within 1 standard error of the minimum. This choice usually leads to a sparser set of selected variables.}
#' \item{index}{A one column matrix with the indices of \code{threshold.min} and \code{threshold.1se} in the sequence of all threshold values.}
#'
#' @export
#'
#' @description Performs k-fold cross-validation for \code{slr} and returns an optimal value for \code{threshold}.
#'
#' @details Please see \code{slr} for details regarding the choice of \code{s0.perc} and \code{zeta}.
#'
#' @author Jing Ma and Kristyn Pantoja.
#'
#' Maintainer: Jing Ma (\url{jingma@fredhutch.org})
#'
#' @seealso \code{slr}
#'
#' @examples
#'
#' HIV <- load_data() # Load HIV data
#' X <- HIV[,1:60]
#' y <- ifelse(HIV[,62] == "Pos", 1, 0)
#' X.adjusted <- sweep(X+1,rowSums(X+1),MARGIN = 1, FUN='/')# zero handling
#'
#' cv.out <- cv.slr(X.adjusted, y, screen.method='wald', cluster.method ='spectral',
#'                  response.type = 'binary', threshold = NULL,type.measure = 'auc',
#'                  trace.it = TRUE, plot = TRUE)
#'
cv.slr <- function(
    X,
    y,
    screen.method = c('correlation', 'wald'),
    cluster.method = c('spectral', 'hierarchical'),
    response.type = c('continuous', 'binary'),
    threshold = NULL,
    s0.perc = 0,
    zeta = 0,
    type.measure = c("default", "mse", "deviance", "class", "auc", "mae", "C", "accuracy"),
    nfolds = 10,
    foldid = NULL,
    weights = NULL,
    trace.it = FALSE,
    plot = FALSE
){
  type.measure = match.arg(type.measure)
  N <- nrow(X)
  p <- ncol(X)

  if (is.null(weights)){
    weights = rep(1, nfolds)
  }

  if (is.null(threshold)) {
    xclr <- apply(X,1,function(a) log(a) - mean(log(a)))
    xclr.centered <- base::scale(t(xclr),center=TRUE, scale=TRUE)

    # determine threshold based on univariate score statistics or correlations
    threshold = sort(
      getFeatureScores(X, y, screen.method, response.type, s0.perc))
  }

  if (is.null(foldid)) {
    foldid = sample(rep(seq(nfolds), length = N))
  } else {
    nfolds = max(foldid)
  }
  if (nfolds < 3){
    stop("nfolds must be bigger than 3; nfolds=10 recommended")
  }

  if (trace.it){
    cat("Training\n")
  }

  outlist = as.list(seq(nfolds))
  for (i in seq(nfolds)) {
    if (trace.it){
      cat(sprintf("Fold: %d/%d\n", i, nfolds))
    }
    which.fold.i = foldid == i
    # x_in <- X[which.fold.i, ,drop=FALSE]
    x_sub <- X[!which.fold.i, ,drop=FALSE]
    y_sub <- y[!which.fold.i]
    outlist[[i]] <- lapply(threshold, function(l) slr(
      X = x_sub, y = y_sub,
      screen.method = screen.method, cluster.method = cluster.method,
      response.type = response.type,
      threshold = l,
      s0.perc = s0.perc, zeta = zeta))
  }

  # collect all out-of-sample predicted values
  #   with the updated code, this is more like a CV matrix
  predmat <- buildPredmat(
    outlist, threshold, X, y, foldid, response.type = response.type,
    type.measure = type.measure)

  cvm <- apply(predmat, 2, stats::weighted.mean, w=weights, na.rm = TRUE)
  cvsd <- apply(predmat, 2, stats::sd, na.rm = TRUE) / sqrt(nfolds)

  out <- list(
    threshold = threshold,
    cvm=cvm,cvsd = cvsd,
    # fit.preval = predmat,
    foldid = foldid
  )

  lamin <- with(out, getOptcv(threshold, cvm, cvsd))

  obj = c(out, as.list(lamin))
  class(obj) = "cv.slr"

  if (plot){
    df = with(obj,
              data.frame(threshold = threshold, l = cvm, cvsd = cvsd))
    if (response.type=='binary' && type.measure=='auc'){
      df = with(obj,
                data.frame(threshold = threshold, l = 1-cvm, cvsd = cvsd))
    }

    suppressWarnings({
      # This warning was about the length of the arrow, which can be zero.
      with(df, plot(threshold, l, ylim=c(min(l - cvsd),max(l+cvsd)), xlab="threshold", ylab=type.measure, col='red', type='b'))
      # Add error bars
      with(df, arrows(x0=threshold, y0=l - cvsd, x1=threshold, y1=l + cvsd, code=3, angle=90, length=0.05))
      # Add vertical lines
      with(obj, abline(v=c(threshold.1se,threshold.min), col=c("blue", "red"), lty=c(2,3), lwd=c(2,2)))
    })

  }
  return(obj)
}

