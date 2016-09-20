estimateQ  <-  function(#Estimates  the  Conditional Expectation  of  Y  Given (A,W).
### Function  for the  estimation of  the conditional  expectation  of \eqn{Y}
### given \eqn{(A,W)}, aka 'Q'.
                        obs,
### A \code{data.frame} of observations, as produced by \code{getSample}.
                        flavor=c("parametric", "lasso"), 
### A  \code{character}  indicating  the   flavor  of  the  procedure,  either
### 'parametric' or 'lasso'.
                        weights,
### The \code{vector} of weights upon which the estimation procedure relies.
                        learnQ=NULL,
### A parametric  model \eqn{{\cal Q}} of conditional  expectations of \eqn{Y}
### given \eqn{(A,W)}  for both flavors  "parametric" and "lasso", given  as a
### \code{formula} or as a \code{function} outputing formulas.
                        Qmin,
### A \code{numeric},  used to  rescale the values  of the outcome  \eqn{Y} by
### using the scaling \code{function} \code{scaleY}.
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
  
  ## Argument 'flavor'
  flavor <- match.arg(flavor)
  learnMode <- switch(flavor,
                      parametric=c("formula", "function"),
                      lasso=c("formula", "function"))
  class <- class(learnQ)
  if (!class%in%learnMode) {
    throw("Argument 'learnQ' should be of class '", learnMode, "', not '", class, "' for flavor: ", flavor)
  }
  if (class=="formula") {
    if (learnQ[[2]]!="Y") {
      throw("Argument 'learnQ' should be a formula with response 'Y', not:", learnQ) 
    }
  }   
  ## Argument 'Qmin'
  Qmin <- Arguments$getNumeric(Qmin, c(0, 1/4))  
 
  ## Argument 'verbose'
  verbose <- Arguments$getVerbose(verbose)

  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Core
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if (flavor=="parametric") {
    Y <- scaleY(obs[, "Y"], thr=Qmin)
    A <- obs[, "A"]
    W <- extractW(obs)

    if (is.function(learnQ)) {
      form <- learnQ()
      learnQ <- form
    }
    
    dat <- data.frame(Y=Y, A=A, W)
    fit <- do.call(glm,
                   list(formula=learnQ,
                        data=dat,
                        family="binomial",
                        weights=weights))
    coeffs <- coef(fit)

    Q <- function(A, W) {
      dat <- data.frame(A=A, W)
      Qy <- predict(fit, newdata=dat, type="response")
      Qy <- pmin(1-2*Qmin, pmax(2*Qmin, Qy))
      attributes(Qy) <- attributes(Y)
      Qy <- scaleY(Qy, reverse=TRUE)
      Qy
    }
    attr(Q, "model") <- coeffs      
  }
  
  if (flavor=="lasso") {
    Y <- scaleY(obs[, "Y"], thr=Qmin)
    A <- obs[, "A"]
    W <- extractW(obs)

    if (is.function(learnQ)) {
      form <- learnQ()
      learnQ <- form
    }
    
    YY <- cbind(1-Y, Y)
    XX <- model.matrix(learnQ, data=data.frame(Y=Y, A=A, W))

    if (ncol(XX)<2) {
      throw("Flavor 'lasso' with such a 'learnQ' does not make any sense:",
            deparse(substitute(form, list(form=learnQ))))
    }

    cv.fits <- glmnet::cv.glmnet(x=XX, y=YY, weights=weights, family="binomial",
                                  intercept=FALSE)#, lambda.min.ratio=1e-3)
    fit <- cv.fits$glmnet.fit
    lambda <- cv.fits$lambda.min ## cv.fits$lambda.1se
    coeffs <- coef(fit, s=lambda)
    
    Q <- function(A, W) {
      XX <- model.matrix(learnQ, data=data.frame(Y=0.5, A=A, W))
      Qy <- glmnet::predict.glmnet(fit, newx=XX, s=lambda, type="response")
      Qy <- pmin(1-2*Qmin, pmax(2*Qmin, Qy))
      attributes(Qy) <- attributes(Y)
      Qy <- scaleY(Qy, reverse=TRUE)
      Qy
    }
    attr(Q, "model") <- coeffs
  }

  
  return(Q)
### Returns \code{Q}, a \code{function} estimating the conditional expectation
### of  \eqn{Y} given  \eqn{(A,W)}.   This function  has  a "model"  attribute
### describing the latter conditional expectation.
}


############################################################################
## HISTORY:
## 2016-02-05
## o Added pmin(., pmax(., ))
## 2014-06-13
## o Adapted to deal with option 'MOR'
## 2014-02-26
## o Created.
############################################################################

