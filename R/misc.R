logit <- qlogis

expit <- plogis

Qbar1 <- function#A Conditional Expectation of \eqn{Y} Given \eqn{(A,W)}
### A  conditional  expectation  of  \eqn{Y}   given  \eqn{(A,W)}  to  use  in
### \code{\link{getSample}}.
(AW,
### A \code{data.frame} of observations,  whose columns contain the components
### of \eqn{W} and the value of \eqn{A}.
 rho=1
### A non-negative \code{numeric}, with default value equal to one.
 ) {
  ##references<<  Chambaz,  van  der  Laan, Scand.  J.  Stat.,  41(1):104--140
  ##(2014).

  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ## Argument 'AW':
  ## AW <- Arguments$getNumerics(AW)
  if (!is.data.frame(AW)) {
    throw("Argument 'AW' should be a data.frame, not:", mode(AW))
  }

  ## Argument 'rho':
  rho <- Arguments$getNumeric(rho, c(0, Inf))
  
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Core
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  U <- AW[, "U"]
  V <- as.integer(AW[, "V"])
  A <- AW[, "A"]
  baseline <- 2*U^2 + 2*U + 1
  out1 <- rho*V + baseline
  out0 <- rho/(1 + V) + baseline
  out <- A*out1 + (1-A)*out0

  return(out)
### Returns a \code{vector} of conditional expectations given \eqn{(A,W)}.
}

Qbar2 <- function#A conditional Expectation of \eqn{Y} Given \eqn{(A,W)}
### A  conditional  expectation  of  \eqn{Y}   given  \eqn{(A,W)}  to  use  in
### \code{\link{getSample}}.
(AW,
### A \code{data.frame} of observations,  whose columns contain the components
### of \eqn{W} and the value of \eqn{A}.
 rho=1
### A non-negative \code{numeric}, with default value equal to 1.5.
 ) {
  ##references<< Chambaz, Zheng, van der Laan, https://hal.archives-ouvertes.fr/hal-01301297.

  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ## Argument 'AW':
  ## AW <- Arguments$getNumerics(AW)
  if (!is.data.frame(AW)) {
    throw("Argument 'AW' should be a data.frame, not:", mode(AW))
  }

  ## Argument 'rho':
  rho <- Arguments$getNumeric(rho, c(0, Inf))
  
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Core
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  U <- AW[, "U"]
  V <- as.integer(AW[, "V"])
  A <- AW[, "A"]
  out1 <- (1+0.75*cos(rho*V*pi*U))/2
  out0 <- (1+0.50*sin(3*rho/V*pi*U))/2
  out <- A*out1 + (1-A)*out0
  
  return(out)
### Returns a \code{vector} of conditional expectations given \eqn{(A,W)}.
}



Vbar1 <- function#A Conditional Variance of \eqn{Y} Given \eqn{(A,W)}
### A  conditional   variance  of   \eqn{Y}  given   \eqn{(A,W)}  to   use  in
### \code{\link{getSample}}.
(AW,
### A \code{data.frame} of observations,  whose columns contain the components
### of \eqn{W} and the value of \eqn{A}.
 rho=1
### A non-negative \code{numeric}, with default value equal to one.
 ) {
  ##references<<  Chambaz,  van der  Laan,  Scand.  J.  Stat.,  41(1):104--140
  ##(2014).

  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  ## Argument 'AW':
  ## AW <- Arguments$getNumerics(AW)
  if (!is.data.frame(AW)) {
    throw("Argument 'AW' should be a data.frame, not:", mode(AW), ".")
  }

  ## Argument 'rho':
  rho <- Arguments$getNumeric(rho, c(0, Inf))
  
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Core
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  U <- AW[, "U"]
  V <- as.integer(AW[, "V"])
  A <- AW[, "A"]
  out1 <- rho + U + V
  out0 <- 1/(rho+V) + U
  out <- (A*out1 + (1-A)*out0)^2

  ## if (!all(out>0)) {
  ##   throw("Conditional variances must be positive.")
  ## }
  
  return(out)
### Returns a \code{vector} of conditional variances given \eqn{(A,W)}.
}

Vbar2 <- function#A Conditional Variance of \eqn{Y} Given \eqn{(A,W)}
### A  conditional   variance  of   \eqn{Y}  given   \eqn{(A,W)}  to   use  in
### \code{\link{getSample}}. 
(AW,
### A \code{data.frame} of observations,  whose columns contain the components
### of \eqn{W} and the value of \eqn{A}.
 rho=0.1
### A non-negative \code{numeric}, with default value equal to one.
 ) {


  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  ## Argument 'AW':
  ## AW <- Arguments$getNumerics(AW)
  if (!is.data.frame(AW)) {
    throw("Argument 'AW' should be a data.frame, not:", mode(AW), ".")
  }

  ## Argument 'rho':
  rho <- Arguments$getNumeric(rho, c(0, Inf))
  
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Core
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  U <- AW[, "U"]
  V <- as.integer(AW[, "V"])
  A <- AW[, "A"]
  out1 <- rho #  + U + V
  out0 <- rho # 1/(rho+V) + U
  out <- (A*out1 + (1-A)*out0)^2

  ## if (!all(out>0)) {
  ##   throw("Conditional variances must be positive.")
  ## }
  
  return(out)
### Returns a \code{vector} of conditional variances given \eqn{(A,W)}.
}


oneOne <- function#Balanced Treatment Mechanism
### Balanced treatment mechanism.
(W
### A \code{data.frame} of covariates.
 ) {
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  ## Argument 'W':
  ## W <- Arguments$getNumerics(W)
  if (!is.data.frame(W)) {
    throw("Argument 'W' should be a data.frame, not:", mode(W))
  }

  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Core
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  return(rep(0.5, nrow(W)))
###Returns a \code{vector} of probabilities, all equal to 1/2.
}

scaleY <- function#Scaling Function.
### Scaling function.
(y,
### A \code{vector} of numerics.
 reverse=FALSE,
### A \code{logical} (default 'FALSE'). If  'TRUE' then the reverse scaling of
### a vector initially scaled using \code{scaleY} is performed.
 wrt=NULL,
### A \code{vector} of numerics from  which the characteristics of the scaling
### are drawn. Defaults to \code{NULL}, meaning that \code{y} is used.
 thr=1e-2
### A  small  positive  \code{numeric}  (default value  1e-2)  specifying  the
### minimum  (resp. 1 minus  the maximum)  value that  the scaled  vector must
### achieve. Only used if \code{reverse} is 'FALSE'.
 ) {

  ## - - - - - - - - - - - - -  - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - -  - - - - - - - - - - - - - - - -
  
  ## Argument 'y':
  y <- Arguments$getNumerics(y)
  
  ## Argument 'reverse':
  reverse <- Arguments$getLogical(reverse)

  test <-identical(names(attributes(y)), c("min", "max", "thr"))
  if (reverse & !test) {
    throw("Cannot unscale a vector which has not been scaled before!")
  }

  ## Argument 'wrt':
  if (!is.null(wrt)) {
    wrt <- Arguments$getNumerics(wrt)
  } else {
    wrt <- y
  }

  ## Argument 'thr':
  thr <- Arguments$getNumerics(thr, c(0,1))
  
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Core
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  if (!reverse) {
    m <- min(wrt)
    M <- max(wrt)
  } else {
    m <- attr(y, "min")
    M <- attr(y, "max")
    thr <- attr(y, "thr")
  }
  t <- (1-thr)/(1-2*thr)
  a <- t*m + (1-t)*M
  b <- t*M + (1-t)*m
  
  if (!reverse) {
    out <- (y - a)/(b-a)
    attr(out, "min") <- m
    attr(out, "max") <- M
    attr(out, "thr") <- thr
  } else {
    out <- a+(b-a)*y
    attributes(out) <- NULL
  }
  return(out)
### Returns a \code{vector} of (un)scaled values.
}

scaleQmat <- function#Scaling Function of a Matrix.
### Scaling function of a Matrix.
(Q,
### A \code{matrix} of numerics to be scaled columnwise based on the values of
### \code{wrt}.
 reverse=FALSE,
### A \code{logical} (default 'FALSE'). If  'TRUE' then the reverse scaling of
### a matrix initially scaled using \code{scaleQmat} is performed.
  wrt,
### A \code{vector} of numerics from  which the characteristics of the scaling
### are drawn.
 thr=1e-2
### A  small  positive  \code{numeric}  (default value  1e-2)  specifying  the
### minimum  (resp. 1 minus  the maximum)  value that  the scaled  matrix must
### achieve. Only used if \code{reverse} is 'FALSE'.
 ) {

  ## - - - - - - - - - - - - -  - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - -  - - - - - - - - - - - - - - - -

  ## Argument 'Q':
  Q <- Arguments$getNumerics(Q)
  if (!is.matrix(Q)) {
    throw("Argument 'Q' should be a matrix, not:", class(Q))
  }

  ## Argument 'reverse':
  reverse <- Arguments$getLogical(reverse)

  if (!reverse) {
    ## Argument 'wrt':
    wrt <- Arguments$getNumerics(wrt)
  }

  test <- all(c("min", "max", "thr") %in% names(attributes(Q)))
  if (reverse & !test) {
    throw("Cannot unscale a matrix which has not been scaled before!")
  }

  
  ## Argument 'thr':
  thr <- Arguments$getNumerics(thr, c(0,1/2))
  
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Core
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  if (!reverse) {
    m <- min(wrt)
    M <- max(wrt)
  } else {
    m <- attr(Q, "min")
    M <- attr(Q, "max")
    thr <- attr(Q, "thr")
  }
  t <- (1-thr)/(1-2*thr)
  a <- t*m + (1-t)*M
  b <- t*M + (1-t)*m
  
  if (!reverse) {
    out <- (Q - a)/(b-a)
    attr(out, "min") <- m
    attr(out, "max") <- M
    attr(out, "thr") <- thr
  } else {
    out <- a+(b-a)*Q
    attr(out, "min") <- NULL
    attr(out, "max") <- NULL
    attr(out, "thr") <- NULL
  }
  return(out)
### Returns a \code{matrix} with (un)scaled columns.
}


## checkCont <- function#Checks Continuity
## ### Checks whether a covariate is continuous or not.
## (W
## ### A \code{numeric} vector of covariates.
##  ) {
##   ## details<<  Somewhat approximate  method to check  whether a  covariate is
##   ##  continuous   or  not  based   on  independent  copies  drawn   from  its
##   ##  distribution.   The  covariate  is  declared continuous  iff  there  are
##   ## \code{n} different values among the \code{n} copies.

##   ## - - - - - - - - - - - - -  - - - - - - - - - - - - - - - - -
##   ## Validate arguments
##   ## - - - - - - - - - - - - - -  - - - - - - - - - - - - - - - -
  
##   ## Argument 'W'
##   W <- Arguments$getNumerics(W)

##   ## - - - - - - - - - - - - -  - - - - - - - - - - - - - - - - -
##   ## Core
##   ## - - - - - - - - - - - - - -  - - - - - - - - - - - - - - - -
##   out <- TRUE
##   test <- all(diff(sort(W, method="quick"))>0)
##   if (!test) {
##     out <- FALSE
##   }
##   return(out)
## ### Returns  a  \code{logical}, 'TRUE'  if  the  covariate  is continuous  and
## ### 'FALSE' otherwise.
## }


blik <- function(x, y)
{
  x*y+(1-x)*(1-y)
}

getOptTm <- function#Computes the Optimal Treatment Mechanism Within a Parametric Model
### Computes the optimal  treatment mechanism within a  given parametric model
### of  treatment  mechanisms.   The  computation  is  based  on  Monte  Carlo
### simulation.
(tm.model=formula(A~1),
### A parametric  model \eqn{{\cal G}} of treatment  mechanisms. The procedure
### targets  the optimal treatment  mechanism within  this model.  Defaults to
### \code{formula(A~1)}.
 Gmin=1e-2,
### A  small positive  \code{numeric}, the  minimum value  of elements  of the
### parametric  model \eqn{{\cal  G}}  of treatment  mechanisms (see  argument
### \code{tm.model}). The maximal value is \code{1-Gmin}.
 piV=c(1/2, 1/3, 1/6),
### Marginal distribution of \eqn{V}. Defaults to \code{c(1/2, 1/3, 1/6)}.
 Qbar=Qbar1,
### A   \code{function},  the   conditional  expectation   of   \eqn{Y}  given
### \eqn{(A,W)}. Defaults to \code{Qbar1}.
 Vbar=Vbar1,
### A   \code{function},   the   conditional   variance   of   \eqn{Y}   given
### \eqn{(A,W)}. Defaults to \code{Vbar1}.
 n=1e3,
### An  \code{integer}, the  sample size  of the  simulated data  set  used to
### compute the optimal treatment mechanism.
 tm.control=glm.control(maxit=500),
### A \code{list}  of options to pass  to \code{glm} for the  targeting of the
### treatment  mechanism  within the  model  defined  by argument  'tm.model'.
### Defaults to \code{glm.control(maxit=500)}.
 verbose=FALSE
### A \code{logical}  or an \code{integer}  indicating the level  of verbosity
### (defaults to 'FALSE').
 ) {
  ##seealso<< getSample, getOptVar
  
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  ## Argument 'tm.model'
  if (!identical(class(tm.model), "formula") | tm.model[[2]]!="A") {
    throw("Argument 'tm.model' should be a formula with response 'A', not:",
          deparse(substitute(form, list(form=tm.model)))) 
  }

  ## Argument 'Gmin'
  Gmin <- Arguments$getNumeric(Gmin, c(0, 1/2))

  ## TRICK
  ## Gmin <- Gmin/2
  
  ## Argument 'Qbar':
  mode <- mode(Qbar);
  if (mode != "function") {
    throw("Argument 'Qbar' should be of mode 'function', not '", mode)
  }
  
  ## Argument 'Vbar':
  mode <- mode(Vbar);
  if (mode != "function") {
    throw("Argument 'Vbar' should be of mode 'function', not '", mode)
  }
  
  ## Argument 'n':
  n <- Arguments$getNumeric(n)

  ## Argument 'verbose'
  verbose <- Arguments$getVerbose(verbose)
  verbose <- less(verbose, 10)

  
  ##
  ## preparing 'targetLink' (just like in 'TSMLCARA')
  ##

  ## link
  linkfun <- function(mu) {# maps 'mu' (probs) to 'eta' (reals)
    logit((mu-Gmin)/(1-2*Gmin))
  }
  ## inverse link
  linkinv <- function(eta) {# maps 'eta' (reals) to 'mu' (probs)
    Gmin + (1-2*Gmin)*expit(eta)
  }
  ## derivative of inverse link wrt eta
  mu.eta <- function(eta) {
    expit.eta <- expit(eta)
    (1-2*Gmin)*expit.eta*(1-expit.eta)
  }
  ## test of validity for eta 
  valideta <- function(eta) {
    TRUE
  }
  ## "deviance residuals" as a function of eta
  dev.resids <- function(y, eta, wt) {
    mu <- linkinv(eta)
    wt*(y/mu + (1-y)/(1-mu))
  }
  ## computing 'mustart'
  initialize <- expression({
    if (is.factor(y)) {
      y <- y != levels(y)[1L]
    }
    n <- rep.int(1, nobs)
    y[weights == 0] <- 0
    if (any(y < 0 | y > 1)) {
      stop("y values must be 0 <= y <= 1")
    }
    mustart <- rep(1/(1+sqrt( sum(weights*(1-y))/sum(weights*y) )),
                   times=nobs)
  })
  link <- "targetLink"
  targetLink <- structure(list(linkfun=linkfun, linkinv=linkinv,
                             mu.eta=mu.eta, valideta=valideta, 
                             dev.resids=dev.resids,
                             initialize=initialize, name=link),
                        class = "link-glm")

  ##  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Core
  ##  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ##
  ## Computing Gstar
  ##

  obs <- getSample(n, tm=oneOne, piV=piV, Qbar=Qbar, Vbar=Vbar)
  AW <- extractAW(obs)
  
  ## targeting 'Gstar'
  verbose && enter(verbose, "Targeting Gstar")

  weights <- (obs[, "Y"]-Qbar(AW))^2

  Gstar <- targetGstar(obs=obs, weights=weights, tm.model=tm.model,
                       targetLink=targetLink, tm.control=tm.control,
                       verbose=verbose)

  
  return(Gstar)
### Returns  the  optimal  treatment  mechanism within  the  parametric  model
### \eqn{{\cal G}}. 
}

getOptVar <- function#Computes the Optimal Variance Given a Parametric Model
### Computes the optimal variance, defined as the variance associated with the
### optimal treatment  mechanism within a given parametric  model of treatment
### mechanisms. The computation is based on Monte Carlo simulation.
(tm.model=formula(A~1),
### A parametric  model \eqn{{\cal G}} of treatment  mechanisms. The procedure
### targets  the optimal treatment  mechanism within  this model.  Defaults to
### \code{formula(A~1)}.
 Gmin=1e-2,
### A  small positive  \code{numeric}, the  minimum value  of elements  of the
### parametric  model \eqn{{\cal  G}}  of treatment  mechanisms (see  argument
### \code{tm.model}). The maximal value is \code{1-Gmin}.
 piV=c(1/2, 1/3, 1/6),
### Marginal distribution of \eqn{V}. Defaults to \code{c(1/2, 1/3, 1/6)}.
 Qbar=Qbar1,
### A   \code{function},  the   conditional  expectation   of   \eqn{Y}  given
### \eqn{(A,W)}. Defaults to \code{Qbar1}.
 Vbar=Vbar1,
### A   \code{function},   the   conditional   variance   of   \eqn{Y}   given
### \eqn{(A,W)}. Defaults to \code{Vbar1}.
 n=1e3,
### An \code{integer}, the  common sample size of the  two simulated data sets
### used  to compute the  optimal treatment  mechanism and  associated optimal
### variance.
 tm.control=glm.control(maxit=500),
### A \code{list}  of options to pass  to \code{glm} for the  targeting of the
### treatment  mechanism  within the  model  defined  by argument  'tm.model'.
### Defaults to \code{glm.control(maxit=500)}.
 verbose=FALSE
### A \code{logical}  or an \code{integer}  indicating the level  of verbosity
### (defaults to 'FALSE').
 ) {
  ##seealso<< getSample, getOptTm
  
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  ## Argument 'tm.model'
  if (!identical(class(tm.model), "formula") | tm.model[[2]]!="A") {
    throw("Argument 'tm.model' should be a formula with response 'A', not:",
          deparse(substitute(form, list(form=tm.model)))) 
  }

  
  ## Argument 'Qbar':
  mode <- mode(Qbar);
  if (mode != "function") {
    throw("Argument 'Qbar' should be of mode 'function', not '", mode)
  }
  
  ## Argument 'Vbar':
  mode <- mode(Vbar);
  if (mode != "function") {
    throw("Argument 'Vbar' should be of mode 'function', not '", mode)
  }
  
  ## Argument 'n':
  n <- Arguments$getNumeric(n)

  ## Argument 'verbose'
  verbose <- Arguments$getVerbose(verbose)
  verbose <- less(verbose, 10)

  

  ##  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Core
  ##  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ##
  ## Computing the optimal treatment mechanism
  ##

  optTm <- getOptTm(tm.model=tm.model,   Gmin=Gmin,
                    piV=piV,  Qbar=Qbar, Vbar=Vbar,
                    n=n, verbose=verbose)
  optVar <- getSample(n=n, tm=optTm,
                      piV=piV,  Qbar=Qbar, Vbar=Vbar,
                      what="ATE")[2]^2
  names(optVar) <- NULL
  
  return(optVar)
### Returns  the  optimal  variance   associated  with  the  parametric  model
### \eqn{{\cal G}}.
}


ruleQbar <- function#Computes the Treatment Rule Associated with Qbar
### Computes the Treatment Rule Associated with Qbar
(Qbar,
### A \code{function}, a conditional expectation of Y given (A,W)
  ...
### Arguments to be passed to 'Qbar'  
) {
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  ## Argument 'Qbar':
  mode <- mode(Qbar)
  if (mode != "function") {
    throw("Argument 'Qbar' should be of mode 'function', not '", mode)
  }
  
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Core
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  rule <- function(W) {
    ## Argument 'W':
    if (!is.data.frame(W)) {
      throw("Argument 'W' should be a data.frame, not:", mode(W))
    }
    Qbar(1, W, ...) - Qbar(0, W, ...) > 0
  }
  return(rule)
###Returns a \code{function}, the rule associated with 'Qbar'.
}

makeLearnQ <- function#Builds a Parametric Model Based on Sample Size
### Builds a Parametric Model Based on Sample Size
() {
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ## No argument!
  
  if(exists("obs", envir=sys.parent())){
    nobs <- nrow(get("obs", envir=sys.parent()))
    deg <- 3 + floor(nobs/500)
    nlev <- ceiling(nobs/250)
  } else{
    deg <- 3
    nlev <- 1
  }
  
  
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Core
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if (nlev==1) {
    learnQ <- as.formula(sprintf("Y ~ 
                                I(A==0)*V*poly(U, %s, coefs=attr(poly(seq(0, 1, 1e-3), %s), \"coefs\")) +
                                I(A==1)*V*poly(U, %s, coefs=attr(poly(seq(0, 1, 1e-3), %s), \"coefs\"))",
                                as.character(deg),
                                as.character(deg),
                                as.character(deg),
                                as.character(deg)))
  } else {
    learnQ <- as.formula(sprintf("Y ~ 
                                I(A==0)*V*(poly(U, %s, coefs=attr(poly(seq(0, 1, 1e-3), %s), \"coefs\"))+
                                I(factor(ceiling(%s*U), levels=1:%s))) +
                                I(A==1)*V*(poly(U, %s, coefs=attr(poly(seq(0, 1, 1e-3), %s), \"coefs\"))+
                                I(factor(ceiling(%s*U), levels=1:%s)))",
                                as.character(deg),
                                as.character(deg),
                                as.character(nlev),
                                as.character(nlev),
                                as.character(deg),
                                as.character(deg),
                                as.character(nlev),
                                as.character(nlev)))
  }
  
  
  return(learnQ)
###Returns a \code{formula} describing the parametric model.
}

makeLearnQ.piecewise <- function#Builds a Parametric Model Based on Sample Size
### Builds a Parametric Model Based on Sample Size
() {
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ## No argument!
  
  if(exists("obs", envir=sys.parent())){
    nobs <- nrow(get("obs", envir=sys.parent()))
    deg <- min(ceiling(nobs/100), 6)
    nlev <- min(ceiling(nobs/100), 10)
  } else{
    deg <- 3
    nlev <- 1
  }
  
  
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Core
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if (nlev==1) {
    learnQ <- as.formula(sprintf("Y ~ 
                                I(A==0)*V*poly(U, %s, coefs=attr(poly(seq(0, 1, 1e-3), %s), \"coefs\")) +
                                I(A==1)*V*poly(U, %s, coefs=attr(poly(seq(0, 1, 1e-3), %s), \"coefs\"))",
                                as.character(deg),
                                as.character(deg),
                                as.character(deg),
                                as.character(deg)))
  } else {
    piecewise <- paste("paste(\"poly(I(U*as.integer(ceiling(",
                       sprintf("%s*U)==\",", nlev),
                       sprintf("1:%s, sep=\"\")", nlev), sep="")
    piecewise <- eval(parse(text=piecewise))
    piecewise <- paste(piecewise, ")), DEG, coefs=attr(poly(seq(0, 1, 1e-3), DEG), \"coefs\"))",
                       sep="", collapse="+")
    replace <- paste("gsub(\"DEG\", ",
                       sprintf("%i, piecewise)", deg), sep="")
    piecewise <- eval(parse(text=replace))
    
    learnQ <- as.formula(sprintf(paste("Y ~ 
                                I(A==0)*V*(poly(U, %s, coefs=attr(poly(seq(0, 1, 1e-3), %s), \"coefs\"))+",
                                piecewise, ")+",
                                "I(A==1)*V*(poly(U, %s, coefs=attr(poly(seq(0, 1, 1e-3), %s), \"coefs\"))+",
                                piecewise, ")", sep=""),
                                as.character(deg),
                                as.character(deg),
                                as.character(deg),
                                as.character(deg)))
  }
  
  
  return(learnQ)
###Returns a \code{formula} describing the parametric model.
}



############################################################################
## HISTORY:
## 2016-02-05
## o Added 'ruleQbar', 'makeLearnQ'
## 2014-04-01
## o Added 'getOptTm' and 'getOptVar'
## 2014-03-25
## o Removed 'LegendreExpandCov'
## 2014-02-26
## o Created.
############################################################################

