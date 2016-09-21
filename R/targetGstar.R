targetTm.glm <- function(formula, data, weights,
                         link,
                         subset, na.action, start = NULL,
                         etastart, mustart, offset,
                         control = list(...),
                         model = TRUE, 
                         x = FALSE, y = TRUE,
                         contrasts = NULL, ...) {
  ## - - - - - - - - - - - 
  ## Adaptation of 'glm.R'
  ## - - - - - - - - - - - 
  call <- match.call()
  ## method
  method <- "targetTm.glm.fit"
  ## extract x, y, etc from the model formula and frame
  if(missing(data)) data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action",
               "etastart", "mustart", "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  if(identical(method, "model.frame")) return(mf)

  if (!is.character(method) && !is.function(method))
      stop("targetTm.glm.fit: invalid 'method' argument")
  ## for back-compatibility in return result
  if (identical(method, "targetTm.glm.fit"))
      control <- do.call("glm.control", control)

  mt <- attr(mf, "terms") # allow model.frame to have updated it

  Y <- model.response(mf, "any") # e.g. factors are allowed
  ## avoid problems with 1D arrays, but keep names
  if(length(dim(Y)) == 1L) {
    nm <- rownames(Y)
    dim(Y) <- NULL
    if(!is.null(nm)) names(Y) <- nm
  }
  ## null model support
  X <- if (!is.empty.model(mt)) model.matrix(mt, mf, contrasts) else matrix(,NROW(Y), 0L)
  ## avoid any problems with 1D or nx1 arrays by as.vector.
  weights <- as.vector(model.weights(mf))
  if(!is.null(weights) && !is.numeric(weights))
      stop("targetTm.glm.fit: 'weights' must be a numeric vector")
  ## check weights and offset
  if( !is.null(weights) && any(weights < 0) )
      stop("targetTm.glm.fit: negative weights not allowed")

  offset <- as.vector(model.offset(mf))
  if(!is.null(offset)) {
    if(length(offset) != NROW(Y))
        stop(gettextf("targetTm.glm.fit: number of offsets is %d should equal %d (number of observations)", length(offset), NROW(Y)), domain = NA)
  }
  ## these allow starting values to be expressed in terms of other vars.
  mustart <- model.extract(mf, "mustart")
  etastart <- model.extract(mf, "etastart")

  ## We want to set the name on this call and the one below for the
  ## sake of messages from the fitter function
  fit <- eval(call(if(is.function(method)) "method" else method,
                   x = X, y = Y, weights = weights, start = start,
                   etastart = etastart, mustart = mustart,
                   offset = offset, control = control, link=link,
                   intercept = attr(mt, "intercept") > 0L))

  ## This calculated the null deviance from the intercept-only model
  ## if there is one, otherwise from the offset-only model.
  ## We need to recalculate by a proper fit if there is intercept and
  ## offset.
  ##
  ## The glm.fit calculation could be wrong if the link depends on the
  ## observations, so we allow the null deviance to be forced to be
  ## re-calculated by setting an offset (provided there is an intercept).
  ## Prior to 2.4.0 this was only done for non-zero offsets.
  if(length(offset) && attr(mt, "intercept") > 0L) {
    throw("Should not happen...")
    fit2 <-
        eval(call(if(is.function(method)) "method" else method,
                  x = X[, "(Intercept)", drop=FALSE], y = Y,
                  weights = weights, offset = offset, 
                  control = control, intercept = TRUE))
    ## That fit might not have converged ....
    if(!fit2$converged)
        warning("targetTm.glm.fit: fitting to calculate the null deviance did not converge -- increase 'maxit'?")
    fit$null.deviance <- fit2$deviance
  }
  if(model) fit$model <- mf
  fit$na.action <- attr(mf, "na.action")
  if(x) fit$x <- X
  if(!y) fit$y <- NULL
  fit <- c(fit, list(call = call, formula = formula,
                     terms = mt, 
                     offset = offset, control = control, method = method,
                     contrasts = attr(X, "contrasts"),
                     xlevels = .getXlevels(mt, mf)))
  class(fit) <- c(fit$class, "glm", "lm")
  fit
}


targetTm.glm.fit <- function (x, y, weights = rep(1, nobs), link,
                              start = NULL,
                              etastart = NULL, mustart = NULL, offset = rep(0, nobs),
                              control = list(), intercept = TRUE)
{
  ## - - - - - - - - - - - - 
  ## Adaptation of 'glm.fit'
  ## - - - - - - - - - - - -
  
  control <- do.call("glm.control", control)
  x <- as.matrix(x)
  xnames <- dimnames(x)[[2L]]
  ynames <- if(is.matrix(y)) rownames(y) else names(y)
  conv <- FALSE
  nobs <- NROW(y)
  nvars <- ncol(x)
  EMPTY <- nvars == 0
  ## define weights and offset if needed
  if (is.null(weights))
      weights <- rep.int(1, nobs)
  if (is.null(offset))
      offset <- rep.int(0, nobs)

  ## get family functions:
  
  ## - - - - 
  ## no more:
  ## - - - - 
  ##
  ## family <- binomial(link=link)
  ## family$initialize <- link$initialize
  ## variance <- family$variance
  ##
  ## - - - - 
  ## instead:
  ## - - - - 
  
  family <- link
  variance <- function(mu){mu*(1-mu)}
  linkinv  <- family$linkinv
  if (!is.function(variance) || !is.function(linkinv) )
      stop("targetTm.glm.fit: 'family' argument seems not to be a valid family object", call. = FALSE)
  dev.resids <- family$dev.resids
  aic <- family$aic
  mu.eta <- family$mu.eta
  unless.null <- function(x, if.null) if(is.null(x)) if.null else x
  valideta <- unless.null(family$valideta, function(eta) TRUE)
  validmu  <- unless.null(family$validmu,  function(mu) TRUE)

  if(is.null(mustart)) {
    ## calculates mustart and may change y and weights and set n (!)
    eval(family$initialize)
    start <- c(mustart[1], rep(0, nvars-1))
  } else {
    mukeep <- mustart
    eval(family$initialize)
    mustart <- mukeep
  }
  if(EMPTY) {
    throw("Should not happen: 'EMPTY' is TRUE in call to 'targetTm.glm.fit'")
    eta <- rep.int(0, nobs) + offset
    if (!valideta(eta))
        stop("targetTm.glm.fit: invalid linear predictor values in empty model", call. = FALSE)
    mu <- linkinv(eta)
    ## calculate initial deviance and coefficient
    if (!validmu(mu))
        stop("targetTm.glm.fit: invalid fitted means in empty model", call. = FALSE)
    dev <- sum(dev.resids(y, mu, weights))
    w <- ((weights * mu.eta(eta)^2)/variance(mu))^0.5
    residuals <- (y - mu)/mu.eta(eta)
    good <- rep_len(TRUE, length(residuals))
    boundary <- conv <- TRUE
    coef <- numeric()
    iter <- 0L
  } else {
    coefold <- NULL
    eta <- 
        if(!is.null(etastart)) etastart
        else if(!is.null(start))
            if (length(start) != nvars)
                stop(gettextf("targetTm.glm.fit: length of 'start' should equal %d and correspond to initial coefs for %s", nvars, paste(deparse(xnames), collapse=", ")),
                     domain = NA)
            else {
              coefold <- start
              offset + as.vector(if (NCOL(x) == 1L) x * start else x %*% start)
            }
        else family$linkfun(mustart)
    mu <- linkinv(eta)

    if (!(validmu(mu) && valideta(eta))) {
        stop("targetTm.glm.fit: cannot find valid starting values: please specify some", call. = FALSE)
      }
    ## calculate initial deviance and coefficient
    devold <- sum(dev.resids(y, mu, weights))
    conv <- FALSE
    
    ## - - - - - - - - - - - - - - - - 
    ## START of tailored optimization procedure
    ## - - - - - - - - - - - - - - - - 
    good <- weights > 0
    varmu <- variance(mu)[good]
    if (any(is.na(varmu))) {
      stop("targetTm.glm.fit: NAs in V(mu)")
    }
    if (any(varmu == 0)){
      stop("targetTm.glm.fit: 0s in V(mu)")
    }
    mu.eta.val <- mu.eta(eta)
    if (any(is.na(mu.eta.val[good]))) {
      stop("targetTm.glm.fit: NAs in d(mu)/d(eta)")
    }
    ## drop observations for which w will be zero
    good <- (weights > 0) & (mu.eta.val != 0)
    if (all(!good)) {
      conv <- FALSE
    }
    
    FUN <- function(par) {
      eta <- drop(x %*% par)
      sum(dev.resids(y, eta, weights))
    }
    GR <- function(par) {
      eta <- drop(x %*% par)
      mu <- linkinv(eta)
      out <- t(weights * ((1-y)/(1-mu)^2 - y/mu^2) * mu.eta(eta)) %*% x
      return(out)
    }
    opt <- optim(par=start, fn=FUN, gr=GR, method="BFGS",
                 control=list(maxit=control$maxit))
    conv <- opt$convergence==0
    nulldev <- FUN(start)
    dev <- opt$value
    coef <- opt$par
    names(coef) <- xnames
    eta <- drop(x %*% coef)
    mu <- linkinv(eta)

    ## - - - - - - - - - - - - - - - - - - - - 
    ## END of tailored optimization procedure
    ## - - - - - - - - - - - - - - - - - - - - 
        
    if (!conv) {
      warning("targetTm.glm.fit: algorithm did not converge", call. = FALSE)
    }
    eps <- 10*.Machine$double.eps
    ##    if (family$family == "binomial") {
    if (any(mu > 1 - eps) || any(mu < eps)) {
      warning("targetTm.glm.fit: fitted probabilities numerically 0 or 1 occurred", call. = FALSE)
    }
    ##}

    ## calculate null deviance -- corrected in glm() if offset and intercept
    wtdmu <-
        if (intercept){
          sum(weights * y)/sum(weights)
        } else {
          linkinv(offset)
        }
    
    list(coefficients = coef, fitted.values = mu,
         linear.predictors = eta, deviance = dev,
         family=family,
         null.deviance = nulldev, 
         prior.weights = weights, 
         y = y, converged = conv)
  }
}

predict.targetTm.glm <- function(object, newdata) {
  predict.lm <- function(object, newdata) {
    tt <- terms(object)
    Terms <- delete.response(tt)
    m <- model.frame(Terms, newdata, na.action=NULL,
                     xlev=object$xlevels)
    if(!is.null(cl <- attr(Terms, "dataClasses"))) .checkMFClasses(cl, m)
    X <- model.matrix(Terms, m, contrasts.arg=object$contrasts)
    offset <- rep(0, nrow(X))
    if (!is.null(off.num <- attr(tt, "offset"))) {
      for(i in off.num) {
        offset <- offset + eval(attr(tt, "variables")[[i+1]], newdata)
      }
    }
    if (!is.null(object$call$offset)) {
      offset <- offset + eval(object$call$offset, newdata)
    }
    beta <- object$coefficients
    predictor <- drop(X %*% beta)
    if (!is.null(offset)) {
      predictor <- predictor + offset
    }
    return(predictor)
  }
  pred <- family(object)$linkinv(predict.lm(object, newdata))
  return(pred)
}




targetGstar  <-  function(#Targets the Optimal Treatment Mechanism When Targeting the Average Treatment Effect (ATE)  Within a Parametric Model 
### Function  to  target  the  optimal  treatment  mechanism  within  a  given
### parametric model  when targeting the  Average Treatment Effect  (ATE), aka
### 'Gstar'.
                          obs,
### A \code{data.frame} of observations, as produced by \code{getSample}.
                          weights,
### The \code{vector} of weights upon which the estimation procedure relies.
                          tm.model=formula(A~1),
### A parametric  model \eqn{{\cal G}} of treatment  mechanisms. The procedure
### targets  the optimal treatment  mechanism within  this model.  Defaults to
### \code{formula(A~1)}.
                          targetLink,
### An  object of class  \code{link-glm} used  to carry  out the  targeting of
### 'Gstar'.
                          tm.control=glm.control(maxit=500),
### A  \code{list} of  options for  the targeting  of the  treatment mechanism
### within   the  model   defined   by  argument   'tm.model'.   Defaults   to
### \code{glm.control(maxit=500)}.
                          ...,
### Additional parameters.
                          verbose=FALSE
### A \code{logical}  or an \code{integer}  indicating the level  of verbosity
### (defaults to 'FALSE').
                          ) {
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ## Argument 'obs'
  obs <- validateArgumentObs(obs)

  ## Argument 'weights'
  weights <- Arguments$getNumerics(weights, c(0, Inf))

  ## Argument 'tm.model'
  if (!identical(class(tm.model), "formula") | tm.model[[2]]!="A") {
    throw("Argument 'tm.model' should be a formula with response 'A', not:", tm.model) 
  }

  ## Argument 'targetLink'
  if (!class(targetLink)=="link-glm") {
    throw("Argument 'targetLink' must be of class 'link-glm', not:", class(targetLink)) 
  }

  ## Argument 'verbose'
  verbose <- Arguments$getVerbose(verbose)

  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Core
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  obsF <- data.frame(A=factor(obs[, "A"], levels=0:1), 
                     extractW(obs))

  Gstar.fit <- do.call(targetTm.glm,
                       list(formula=tm.model, link=targetLink,
                            data=obsF, weights=weights,
                            control=tm.control))

  coeffs <- coef(Gstar.fit)
  Gstar <- function(W) {
    predict.targetTm.glm(Gstar.fit, W)
  }

  attr(Gstar, "model") <- coeffs
  
  return(Gstar)
### Returns \code {Gstar}, a \code{function} of the same form as \code{oneOne}
### coding  the   targeted  optimal  treatment  mechanism   within  the  given
### parametric model.   This function has  a "model" attribute  describing the
### latter optimal treatment mechanism.
}

############################################################################
## HISTORY:
## 2016-09-16
## o Created.
############################################################################

