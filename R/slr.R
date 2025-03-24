# Function to obtain principal balances from principal components
#  C, a compositional data matrix with rows being samples and columns being features
#  ncomp, number of components to compute
#  angle, TRUE or FALSE, a flag that determines how the principal balance is calculated. See details. Defalut is \code{TRUE}.
# details
#  given a coda set C of n by p, this function returns the balance with D parts
#  that is closest in angle to the first PC.
#  If "angle=FALSE" then the balance that maximizes the variance of scores is
#  returned.
#
#  This function is modified from the fBalChipman.R function.
PBfromPC<-function(C, ncomp=2, angle=TRUE){
  # columns
  col <- dim(C)[2]
  nbal <- col-1
  # clr-transform
  clrC<-log(C) - rowMeans(log(C))

  # PCs
  pcClr <- stats::prcomp(clrC)
  out <- matrix(0,col,ncomp)
  colnames(out) <- paste0("z",1:ncomp)
  if (!is.null(colnames(C))){
    rownames(out) <- colnames(C)
  }
  vep <- rep(0,ncomp)
  for (k in 1:ncomp){
    #k-th PC
    pcClr1<-pcClr$rotation[,k]

    balsig<-sign(pcClr1)
    bal<-matrix(0,nbal,col)
    colnames(bal)<-colnames(C)

    # balances associated to the PCs
    # first balance
    bal[1,pcClr1==max(pcClr1)]<-1
    bal[1,pcClr1==min(pcClr1)]<--1
    numbal=1

    # other balance
    if (col>2){
      numbal=numbal+1
      while (numbal<col){
        bal[numbal,]<-bal[numbal-1,]
        useonly<-(bal[numbal-1,]==0)
        if (sum(useonly)==0){
          # no more zeros
          break
        } else {
          bal[numbal,abs(pcClr1)==max(abs(pcClr1[useonly]))]<-balsig[abs(pcClr1)==max(abs(pcClr1[useonly]))]
          numbal=numbal+1
        }
      }#end while
    }#end if

    # coefficients & angle
    VarSBP<-rep(0,nbal)
    for (f in 1:nbal) {
      den<-sum(bal[f,]==-1)
      num<-sum(bal[f,]==1)
      bal[f,bal[f,]==1] <- sqrt(den/((den+num)*num))
      bal[f,bal[f,]==-1] <- -sqrt(num/((den+num)*den))
      # angle between the balance and the PC
      VarSBP[f]<-abs(sum(bal[f,]*pcClr1))
    }
    # log-transform
    lC<-as.matrix(log(C))
    mvar=stats::var(as.vector(lC%*%bal[VarSBP==max(VarSBP),]))

    if (!angle) {
      # calculate variance in the balance direction
      VarSBP<-rep(0,nbal)

      for (i in 1:nbal){
        Proj<-as.vector(lC%*%(bal[i,]))
        VarSBP[i]<-stats::var(Proj)
      }# end for
      mvar=max(VarSBP)
    }# end if
    out[,k] <- bal[VarSBP==max(VarSBP),]
    vep[k] <- mvar
  } # end for

  # return results
  return(list(bal=out,var=vep))
}

# Function to compute the univariate feature scores.
getFeatureScores = function(x, y, covar=NULL, family){
  n = nrow(x)

  ## Compute univariate statistic
  xclr <- apply(x, 1, function(a) log(a) - mean(log(a)))
  xclr.centered <- base::scale(t(xclr),center=TRUE,scale=TRUE)

  if (family=='gaussian'){
    fs <- rep(0,ncol(xclr.centered))
    for (j in 1:ncol(xclr.centered)){
      df <- data.frame(cbind(y=y,x=xclr.centered[,j],covar))
      fit <- stats::lm(y~.,data=df)
      fs[j] <- stats::coef(summary(fit))[2,3]
    }
    # t-distribution with n-2 df
    fs <- stats::pt(abs(fs), df = n-2)
  } else if (family=='binomial'){
    fs <- rep(0,ncol(xclr.centered))
    for (j in 1:ncol(xclr.centered)){
      df <- data.frame(cbind(x=xclr.centered[,j],covar))
      fit <- glm(y~.,data=df,family=binomial(link='logit'))
      fs[j] <- coef(summary(fit))[2,3]
    }
    fs <- stats::pnorm(abs(fs))
  } else if (family=='cox'){
    if (!survival::is.Surv(y)){
      y <- survival::Surv(y)
    }
    fs <- rep(0,ncol(xclr.centered))
    for (j in 1:ncol(xclr.centered)){
      df <- data.frame(cbind(x=xclr.centered[,j],covar))
      fit <- survival::coxph(y~.,data=df)
      fs[j] <- coefficients(summary(fit))[1,4]
    }
    fs <- stats::pnorm(abs(fs))
  }
  return(fs)
}


#' Supervised Log Ratio for supervised principal balance analysis
#'
#' @param x A matrix of nobs by nvars strictly positive relative abundances (\code{n} by \code{p}). This matrix should not contain any zeros. See Details section.
#' @param y Response variable. Quantitative for \code{family="gaussian"}. For \code{family="binomial"} should be either a factor with two levels. For \code{family="cox"}, preferably a \code{Surv} object from the survival package: see Details section for more information.
#' @param family Type of the response variable: could be \code{gaussian}, \code{binomial} or \code{cox}.
#' @param covar A numeric matrix, containing additional covariates that you want to adjust for. If \code{NULL}, only feature abundances in \code{x} are used to predict the response variable.
#' @param threshold A nonnegative constant between 0 and 1. If \code{NULL}, then no variable screening is performed.
#' @param method Method used to compute the principal balance. Could be \code{spectral clustering}, \code{hierarchical clustering} or \code{constrained PC}. See Details.
#' @param zeta A small positive value used to perturb the similarity matrix when \code{method='spectral'}. Default is 0. See details.
#' @param feature.scores Optional vector of the feature-level score statistics. If \code{NULL}, the function will compute the feature scores internally.
#' @param Aitchison.var Optional matrix of the Aitchison variation matrix. If \code{NULL}, the function will compute the Aitchison variation internally. Only used if \code{method} is spectral clustering or hierarchical clustering.
#'
#' @description \code{slr} fits a balance regression model in which the balance is defined by a sparse set of input variables.
#'
#' @details \code{slr} is a computationally efficient method for selecting balance biomarkers.
#' \code{slr} consists of two main steps. First, it computes the univariate wald statistic for regressing \code{y} onto the centered log ratio transformed \code{x} and
#' forms a reduced data matrix with variables whose univariate wald statistic exceed a \code{threshold} in absolute value
#' (the threshold is chosen via cross-validation). Next, \code{slr} uses the leading principal balance (see details in Martín-Fernández et al. 2018) on the reduced data to
#'  define a balance biomarker. \code{slr} can be viewed as a semi-supervised approach: it leverages the response variable to filter out inactive features
#'  while using principal balance analysis to capture latent structures in the data. To obtain the leading principal balance, \code{slr} uses the constrained PC approach
#'  (\code{method='PC'}) or clustering (\code{'spectral'}, \code{'hierarchical'}). For clustering based methods, \code{slr} uses the Aitchison variation matrix
#'  on the reduced data as the distance matrix. When spectral clustering is used, it can be beneficial to perturb the Aitchison similarity matrix,
#'  obtained from the Aitchison variation, with a small positive value to improve the clustering performance. This leads to regularized spectral clustering. For more details on
#'  regularized spectral clustering of networks, see Amini et al. (2013). When \code{method='hierarchical'}, ward hierarchical clustering with the \code{ward.D2} agglomeration is used.
#'
#' By default, \code{slr} only uses the leading principal balance of the reduced data matrix, but more components can be used for more accurate prediction.
#' However, for easier interpretation and faster computation, we recommend only using the leading principal balance, which generally works sufficiently well.
#'
#' The design matrix \code{x} cannot contain any zeros. In general, zeros in raw data can be imputed using a pseudocount or
#' Bayesian methods prior to applying \code{slr}.
#'
#'
#' @return An object of class \code{"slr"} is returned, which is a list with
#'
#' \item{threshold}{The threshold used for variable screening.}
#' \item{zeta}{The parameter used to regularize spectral clustering.}
#' \item{family}{Type of the response variable.}
#' \item{sbp}{The binary partition of selected variables. A positive value of 1 indicates the variable is in the numerator, while -1 indicates a denominator variable.}
#' \item{log-contrast coefficients}{The coefficients in the linear log-contrast model.}
#' \item{feature.scores}{The feature scores for all variables.}
#' \item{glm}{The generalized linear model fit using balances as predictors.}
#' @export
#'
#' @author Jing Ma and Kristyn Pantoja.
#'
#' Maintainer: Jing Ma (\url{crystal.jing.ma@gmail.com})
#'
#' @seealso \code{cv.slr}
#'
#' @references
#' Martín-Fernández, J. A., Pawlowsky-Glahn, V., Egozcue, J. J., & Tolosona-Delgado, R. (2018). Advances in principal balances for compositional data. Mathematical Geosciences, 50, 273-298.
#'
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
#' fit <- slr(X.adjusted, y, method ='PC',
#'            family = 'binomial', threshold = 0.9)
#' fit$sbp
slr = function(
    x,
    y,
    covar = NULL,
    family = c('gaussian','binomial','cox'),
    threshold,
    method = c('constrained PC', 'spectral clustering', 'hierarchical clustering'),
    zeta = 0,
    feature.scores = NULL,
    Aitchison.var = NULL
){
  this.call <- match.call()
  family <- match.arg(family)
  ncomp <- 1 # number of balances to use

  if (family == 'cox'){
    if (!survival::is.Surv(y)){
      y <- survival::Surv(y)
    }
    if(!is.matrix(y)||!all(match(c("time","status"),dimnames(y)[[2]],0)))
      stop("Cox model requires a matrix with columns 'time' (>0)
           and 'status'  (binary) as a response; a 'Surv' object suffices",call.=FALSE)
  }
  if(!("data.frame" %in% class(x))) x = data.frame(x)

  n <- nrow(x)
  p <- ncol(x)

  if(is.null(feature.scores)){
    feature.scores = getFeatureScores(x, y, covar, family)
  }
  which.features <- (abs(feature.scores) > threshold)

  object <- list(threshold=threshold, zeta=zeta, family=family)

  if (sum(which.features)<2){
    # Fit an intercept or covariate only regression model
    if (is.null(covar)){
      df <- data.frame(y = y)
    } else {
      df <- data.frame(covar,y=y)
    }
    if (family=='gaussian'){
      model.train <- lm(y~.,data=df)
    } else if (family=='binomial'){
      model.train <- glm(y~.,data=df,family=binomial(link='logit'))
    } else if (family=='cox'){
      model.train <- NULL #'cannot fit a cox proportional hazards model with no covariate'
      if (is.null(covar)){
        model.train <- NULL
      } else {
        model.train <- survival::coxph(y~.,data=data.frame(covar))
      }
    }
    object$sbp <- NULL
    object$`log-contrast coefficients` <- rep(0,p)

  } else {
    x.reduced <- x[,which.features] # reduced data matrix
    # if (ncomp >= sum(which.features)){
    #   warning('not enough features to partition; reduce ncomp.')
    #   ncomp <- sum(which.features) - 1
    # }

    if (sum(which.features)==2){
      # hard code a two-part partition
      sbp.est <- matrix(c(1,-1),ncol=1)
      colnames(sbp.est) <- "z1"
      rownames(sbp.est) <- colnames(x)[which.features]
    } else {
      if (method %in% c("constrained PC","PC","C")){
        # use constrained PCA to obtain the principal balance
        estPB <- PBfromPC(x.reduced, ncomp = 1, angle = TRUE)
        sbp.est <- sign(estPB$bal[,1:ncomp,drop=FALSE])
      } else {
        if (is.null(Aitchison.var)){
          Aitchison.var = getAitchisonVar(x.reduced)
        } else {
          Aitchison.var = Aitchison.var[which.features,which.features]
        }

        if(method %in% c("spectral clustering","spectral",'S')){
          # if (ncomp>1){
          #   stop('cannot perform spectral clustering with ncomp greater than 1.')
          # }
          ## Perform spectral clustering
          Aitchison.sim <- max(Aitchison.var) - Aitchison.var
          sbp.est <- as.matrix(spectral.clust(Aitchison.sim, k=2, zeta = zeta))
          colnames(sbp.est) <- 'z1'
        } else if (method %in% c("hierarchical clustering","hierarchical","H")){
          ## Perform ward hierarchical clustering
          htree.est <- hclust(dist(Aitchison.var),method='ward.D2')
          sbp.est <- sbp.fromHclust(htree.est)[, 1:ncomp, drop=FALSE]
        } else{
          stop("invalid method arg was provided!!")
        }
      }
    }

    balance <- balance.fromSBP(x.reduced, sbp.est)
    df <- data.frame(cbind(balance,covar))
    colnames(df)[1:ncomp] <- paste0('z',1:ncomp)

    # model fitting
    if (family=='gaussian'){
      model.train <- lm(y~.,data=df)
      if(any(coefficients(model.train)[2:(ncomp+1)] < 0)){
        sbp.est =  t(t(sbp.est) * sign(coefficients(model.train)[2:(ncomp+1)]))
        balance <- balance.fromSBP(x.reduced, sbp.est)
        df <- data.frame(cbind(balance,covar))
        colnames(df)[1:ncomp] <- paste0('z',1:ncomp)
        model.train <- lm(y~.,data=df)
      }
      theta <- coef(model.train)[2:(ncomp+1)]
    } else if (family=='binomial'){
      model.train <- glm(y~.,data=df,family=binomial(link='logit'))
      if(any(coefficients(model.train)[2:(ncomp+1)] < 0)){
        sbp.est =  t(t(sbp.est) * sign(coefficients(model.train)[2:(ncomp+1)]))
        balance <- balance.fromSBP(x.reduced, sbp.est)
        df <- data.frame(cbind(balance,covar))
        colnames(df)[1:ncomp] <- paste0('z',1:ncomp)
        model.train <- glm(y~.,data=df,family=binomial(link='logit'))
      }
      theta <- coef(model.train)[2:(ncomp+1)]
    } else if (family=='cox'){
      model.train <- survival::coxph(y~.,data=df)
      if(any(coefficients(model.train)[1:ncomp] < 0)){
        sbp.est =  t(t(sbp.est) * sign(coefficients(model.train)[1:ncomp]))
        balance <- balance.fromSBP(x.reduced, sbp.est)
        df <- data.frame(cbind(balance,covar))
        colnames(df)[1:ncomp] <- paste0('z',1:ncomp)
        model.train <- survival::coxph(y~.,data=df)
      }
      theta <- coef(model.train)[1:ncomp]
    }
    object$sbp <- sbp.est
    # object$ncomp <- ncomp
    fullSBP <- matrix(0, nrow = p, ncol = ncomp)
    rownames(fullSBP) = colnames(x)
    for (k in 1:ncomp){
      fullSBP[match(names(sbp.est[,k]), rownames(fullSBP)),k] = sbp.est[,k]
    }
    colnames(fullSBP) <- colnames(sbp.est)

    object$`log-contrast coefficients` <- balReg2llc(theta,fullSBP)
  }
  object$feature.scores <- feature.scores
  object$glm <- model.train

  class(object) <- 'slr'
  return(object)
}

#' Predict Method for slr Fits
#' @description
#' Obtain predictions from a fitted generalized linear model object
#' @param object a fitted object of class inheriting from `glm`
#' @param newdata optionally, a data frame in which to look for variables with which to predict. If omitted, the fitted linear predictors are used.
#' @param family Type of the response variable: could be \code{gaussian}, \code{binomial} or \code{cox}. This helps determine the prediction model.
#' @return Predictions evaluated at \code{newdata}.
#' @details
#'
#' If \code{newdata} is omitted the predictions are based on the data used for the fit.
#'
#' @export predict.slr
#'
#' @importFrom stats binomial
#' @importFrom stats deviance
#' @importFrom stats coef
#' @importFrom stats lm
#' @importFrom stats glm
#' @importFrom stats coefficients
#' @importFrom stats prcomp
#' @importFrom stats predict
#' @importFrom stats var
#' @importFrom stats dist
#' @importFrom stats hclust
#'
#' @examples
#'
#' HIV <- load_data() # Load HIV data
#' X <- HIV[,1:60]
#' y <- ifelse(HIV[,62] == "Pos", 1, 0)
#' X.adjusted <- sweep(X+1,rowSums(X+1),MARGIN = 1, FUN='/') # zero handling
#'
#' # Run slr ----
#' fit <- slr(X.adjusted, y, method ='PC',
#'            family = 'binomial', threshold = 0.9)
#'
#' predict.slr(fit,family = 'binomial')
predict.slr <- function(
    object,newdata = NULL,family=c('cox','gaussian','binomial')
){
  # prediction will be based on the canonical space
  if (missing(newdata) || is.null(newdata)) {
    newdata = object$glm$data
    predictor <- switch (family,
                         cox = predict(object$glm,type='lp'),
                         gaussian = predict(object$glm),
                         binomial = predict(object$glm,type='response')
    )
    return(as.numeric(predictor))
  } else {
    if (family=='cox'){
      if (is.null(object$sbp) && is.null(object$glm)){
        stop('cannot perform prediction in a null cox model!')
      } else if (is.null(object$sbp) && !is.null(object$glm)){
        # covariates only model
        df <- data.frame(newdata[,colnames(newdata) %in% names(coef(object$glm)),drop=FALSE])
        predictor <- predict(object$glm,newdata=df,type='lp')
      } else {
        ncomp <- ncol(object$sbp)
        newdata.reduced <- newdata[,colnames(newdata) %in% rownames(object$sbp),drop=FALSE]
        new.balance <- balance.fromSBP(newdata.reduced,object$sbp)
        df <- data.frame(new.balance)
        if (length(coef(object$glm))>ncomp) {
          # also include covariates
          df <- data.frame(new.balance,newdata[,colnames(newdata) %in% names(coef(object$glm)),drop=FALSE])
        }
        predictor <- predict(object$glm,newdata=df,type='lp')

      }
    } else if (family=='binomial'){
      if (is.null(object$sbp) && object$glm$rank==1){
        # intercept only model
        predictor <- predict(object$glm,newdata=data.frame(rep(1,nrow(newdata))),type='response')
      } else if (is.null(object$sbp) && object$glm$rank>1){
        # covariates only model
        df <- data.frame(newdata[,colnames(newdata) %in% names(object$glm$model),drop=FALSE])
        predictor <- predict(object$glm,newdata=df,type='response')
      } else {
        ncomp <- ncol(object$sbp)
        newdata.reduced <- newdata[,colnames(newdata) %in% rownames(object$sbp)]
        new.balance <- balance.fromSBP(newdata.reduced,object$sbp)
        df <- data.frame(new.balance)
        if (length(coef(object$glm))>ncomp+1) {
          # also include covariates
          df <- data.frame(new.balance,newdata[,colnames(newdata) %in% names(coef(object$glm)),drop=FALSE])
        }
        predictor <- predict(object$glm,newdata=df,type='response')
      }
    } else if (family=='gaussian'){
      if (is.null(object$sbp) && object$glm$rank==1){
        # intercept only model
        predictor <- predict(object$glm,newdata=data.frame(rep(1,nrow(newdata))))
      } else if (is.null(object$sbp) && object$glm$rank>1){
        # covariates only model
        df <- data.frame(newdata[,colnames(newdata) %in% colnames(object$glm$model),drop=FALSE])
        predictor <- predict(object$glm,newdata=df)
      } else {
        ncomp <- ncol(object$sbp)
        newdata.reduced <- newdata[,colnames(newdata) %in% rownames(object$sbp)]
        new.balance <- balance.fromSBP(newdata.reduced,object$sbp)
        df <- data.frame(new.balance)
        if (length(coef(object$glm))>ncomp+1) {
          # also include covariates
          df <- data.frame(new.balance,newdata[,colnames(newdata) %in% colnames(object$glm$model),drop=FALSE])
        }
        predictor <- predict(object$glm,newdata=df)
      }
    }
    as.numeric(predictor)
  }
}

# this version is nfolds x n.thresh, like the other methods
buildPredmat <- function(
    outlist,threshold,x,y,covar,foldid,family,type.measure
){
  nfolds = max(foldid)
  predmat = matrix(NA, nfolds, length(threshold))
  for(i in 1:nfolds){ # predict for each fold
    which = foldid == i
    y.i = y[which,drop=FALSE]
    fitobj = outlist[[i]]
    x.i = x[which, , drop=FALSE]
    covar.i = covar[which, , drop=FALSE]
    df.i <- x.i
    if (!is.null(covar.i)){
      df.i <- cbind(x.i,covar.i)
    }
    predy.i = sapply(
      fitobj, function(a) predict.slr(
        a,newdata=df.i,family=family))

    if (family =='gaussian'){ # mse, to be minimized
      if(type.measure != "mse"){
        stop("if family is gaussian, then type.measure must be mse!!")
      }
      predmat[i, ] <- apply(predy.i, 2, function(j)  mean((as.numeric(y.i)-j)^2))
    }
    if (family=='binomial'){
      if(!(type.measure %in% c("accuracy", "auc"))){
        stop("if family is binomial, then type.measure must be either accuracy or auc!!")
      }
      predmat[i,] <- apply(predy.i, 2, function(j)
        if(type.measure == "accuracy"){# accuracy, minimize the # that don't match
          mean((j > 0.5) != y.i)
        } else if(type.measure == "auc"){# auc, minimize 1 - auc
          tryCatch({
            1 - pROC::auc(
              y.i,j, levels = c(0, 1), direction = "<", quiet = TRUE)
          }, error = function(e){return(NA)}
          )
        })
    }

    if (family=='cox'){
      # use C index or deviance
      if(!(type.measure %in% c("deviance", "C"))){
        stop("if family is cox, then type.measure must be either deviance or C!!")
      }
      predmat[i,] <- apply(predy.i, 2, function(j)
        if (type.measure=='deviance'){
          deviance(pred = j,y = y.i)
        } else if (type.measure =='C'){
          Cindex(pred = j, y = y.i)
        })
    }
  }
  predmat
}

#  pred linear predictions from a cox model
#  y a n by 2 response matrix with the first column recording survival time and the second column recording status
Cindex <- function(pred,y,weights=rep(1,nrow(y))){
  ###  This function links to the concordance function in the survival package
  if(!survival::is.Surv(y))y=survival::Surv(y[,"time"],y[,"status"])
  f=-pred
  if(missing(weights))
    survival::concordance(y~f)$concordance
  else
    survival::concordance(y~f,weights=weights)$concordance
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
  index = matrix(c(idmin,id1se),2,1,dimnames=list(c("min","1se"),"threshold"))

  list(
    threshold.min = threshold.min,
    threshold.1se = threshold.1se,
    index = index
  )
}


#' Cross Validation for Supervised Log Ratio
#'
#' @param x A sample by variable matrix of strictly positive relative abundances (\code{n} by \code{p}). This matrix should not contain any zeros.
#' @param y The response vector of length \code{n}.
#' @param covar A numeric matrix, containing additional covariates that you want to adjust for. If NULL, only the microbiome abundances are used to predict the response variable.
#' @param family Type of the response variable: could be \code{gaussian}, \code{binomial} or \code{cox}.
#' @param covar A numeric matrix, containing additional covariates that you want to adjust for. If \code{NULL}, only feature abundances in \code{x} are used to predict the response variable.
#' @param threshold Nonnegative values between 0 and 1 used to screen variables. If \code{NULL}, the function will compute the threshold values internally.
#' @param method Method used to compute the principal balance. Could be \code{spectral clustering}, \code{hierarchical clustering} or \code{constrained PC}. See Details.
#' @param zeta A small positive value used to perturb the Aitchison similarity matrix when spectral clustering is used. Default is 0. See Details.
#' @param type.measure Loss used for cross-validation. If \code{family='gaussin'}, then \code{type.measure="mse"}. For two-class logistic regression, \code{type.measure="auc"}. For cox regression, \code{type.measure="C"}.
#' @param nfolds Number of folds. Default is 5.
#' @param foldid An optional vector of values between 1 and \code{nfold} identifying which fold each observation is in. If supplied, \code{nfold} can be missing.
#' @param weights Observation weights. Default is 1 per observation.
#' @param trace.it If \code{trace.it=TRUE}, then progress bars are displayed.
#' @param plot If \code{TRUE}, then a visualization of the cross-validation results is displayed.
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
#' @details The function runs \code{slr} \code{nfolds} times to compute the fit with each of the folds omitted. The error is accumulated, and
#' the average error and standard deviation over the folds is computed. Note that the results of \code{cv.slr} are random, since the folds are selected at random.
#' Users can reduce this randomness by running \code{cv.slr} many times and averaging the error curves.
#'
#' Please see \code{slr} for details regarding the choice of method and choice of \code{zeta}.
#'
#' @author Jing Ma and Kristyn Pantoja.
#'
#' Maintainer: Jing Ma (\url{crystal.jing.ma@gmail.com})
#'
#' @seealso \code{slr}
#'
#' @examples
#' library(compositions)
#' HIV <- load_data() # Load HIV data
#' X <- HIV[,1:60]
#' y <- ifelse(HIV[,62] == "Pos", 1, 0)
#' X.adjusted <- sweep(X+1,rowSums(X+1),MARGIN = 1, FUN='/')# zero handling
#'
#' cv.out <- cv.slr(X.adjusted, y, method ='PC',
#'                  family = 'binomial',type.measure = 'auc',
#'                  trace.it = TRUE)
cv.slr <- function(
    x, y, covar = NULL,
    family=c('cox','gaussian','binomial'),
    threshold = NULL,
    method = c('spectral clustering', 'hierarchical clustering', 'constrained PC'),
    zeta = 0,
    type.measure = c(
      "default", "mse", "deviance", "class", "auc", "mae", "C", "accuracy"
    ),
    nfolds = 5,
    foldid = NULL,
    weights = NULL,
    trace.it = FALSE,
    plot = FALSE
){
  x <- data.frame(x)
  type.measure = match.arg(type.measure)
  n <- nrow(x)
  p <- ncol(x)
  family = match.arg(family)
  if (family=='cox' && !survival::is.Surv(y)){
    y = response.coxnet(y)
  }

  if (is.null(weights)){
    weights = rep(1, nfolds)
  }

  if (is.null(threshold)) {
    xclr <- apply(x,1,function(a) log(a) - mean(log(a)))
    xclr.centered <- base::scale(t(xclr),center=TRUE, scale=TRUE)

    # determine threshold based on univariate score statistics or correlations
    threshold = sort(getFeatureScores(x, y, covar, family))
  } else {
    threshold <- sort(threshold)
  }

  if (is.null(foldid)) {
    foldid = sample(rep(seq(nfolds), length = n))
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
    x_sub <- x[!which.fold.i, , drop=FALSE]
    y_sub <- y[!which.fold.i, drop=FALSE]
    covar_sub <- covar[!which.fold.i,,drop=FALSE]
    # x.unlab_sub = NULL

    feature.scores = getFeatureScores(x_sub, y_sub, covar_sub, family)
    Aitchison.var = getAitchisonVar(x_sub)
    outlist[[i]] <- lapply(threshold, function(l) {
      slr(x = x_sub, y = y_sub, covar = covar_sub,
          family = family,
          threshold = l,
          method = method,
          zeta = zeta,
          feature.scores = feature.scores,
          Aitchison.var = Aitchison.var)
    })
  }


  if (family=='cox'){
    # some models are NULL and need to be removed
    check.NULL <- sapply(outlist, function(b) sapply(b, function(a) is.null(a$glm)))
    if (any(rowSums(check.NULL)>0)){
      id2remove <- which(rowSums(check.NULL)>0)
      threshold <- threshold[-id2remove]
      outlist <- lapply(outlist, function(b) b[-id2remove])
    }
  }
  # collect all out-of-sample predicted values
  #   with the updated code, this is more like a CV matrix
  predmat <- buildPredmat(
    outlist, threshold, x, y, covar, foldid, family = family,
    type.measure = type.measure)

  cvm <- apply(predmat, 2, stats::weighted.mean, w=weights, na.rm = TRUE)
  cvsd <- apply(predmat, 2, stats::sd, na.rm = TRUE) / sqrt(nfolds)

  out <- list(
    threshold = threshold,
    cvm = cvm, cvsd = cvsd,
    fit.preval = predmat,
    foldid = foldid
  )

  lamin <- with(out, getOptcv(threshold, cvm, cvsd))

  obj = c(out, as.list(lamin))
  class(obj) = "cv.slr"

  if (plot){
    df = with(obj,
              data.frame(threshold = threshold, l = cvm, cvsd = cvsd))
    if (family=='binomial' && type.measure=='auc'){
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
  obj
}

