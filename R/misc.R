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
  ##seealso<< Qbar2
  
  ##references<<  Chambaz, van  der  Laan, Scand.   J.  Stat.,  41(1):104--140
  ##(2014).

   ##details<< This conditional expectation of \eqn{Y} given \eqn{(A,W)} where
   ##\eqn{W=(U,V)} is  \deqn{\rho * V  + 2*U^2 + 2*U  + 1} when  \eqn{A=1} and
   ##\deqn{\rho/(1+V) + 2*U^2 + 2*U  + 1} when \eqn{A=0}.
   
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

Qbar3 <- function#A conditional Expectation of \eqn{Y} Given \eqn{(A,W)}
### A  conditional  expectation  of  \eqn{Y}   given  \eqn{(A,W)}  to  use  in
### \code{\link{getSample}}.
(AW,
### A \code{data.frame} of observations,  whose columns contain the components
### of \eqn{W} and the value of \eqn{A}.
 rho=1,
### A non-negative \code{numeric}, with default value equal to 1.
  tau=0.1
### A non-negative \code{numeric}, with default value equal to 0.1.

) {
  ##seealso<< Qbar1, Qbar2

  ##references<<       Chambaz,        Zheng,       van        der       Laan,
  ##https://hal.archives-ouvertes.fr/hal-01301297.

   ##details<< This conditional expectation of \eqn{Y} given \eqn{(A,W)} where
   ##\eqn{W=(U,V)} is \deqn{\frac{1}{2} + \frac{3}{8}*\cos(\rho*V*\pi*U)} when
   ##\eqn{A=1} and \deqn{\frac{1}{2}  + \frac{1}{4}*\sin(3*\rho/V*\pi*U)} when
   ##\eqn{A=0} PROVIDED THAT  the absolute value of the  corresponding blip is
   ##larger than 'tau'. If the blip is between 0 and tau, then the conditional
   ##expectation  of \eqn{Y}  given \eqn{(A=0,W)}  is set  to the  conditional
   ##expectation of  \eqn{Y} given  \eqn{(A=1,W)} minus tau.   If the  blip is
   ##between -tau  and 0,  then the conditional  expectation of  \eqn{Y} given
   ##\eqn{(A=1,W)}  is set  to the  conditional expectation  of \eqn{Y}  given
   ##\eqn{(A=0,W)} minus tau. This is a somewhat brutal way to ensure that the
   ##margin assumption is met.

   
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

  ## Argument 'tau':
  tau <- Arguments$getNumeric(tau, c(0, 1/2))

   
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Core
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  U <- AW[, "U"]
  V <- as.integer(AW[, "V"])
  A <- AW[, "A"]
  out1 <- (1+0.75*cos(rho*V*pi*U))/2
  out0 <- (1+0.50*sin(3*rho/V*pi*U))/2

  blip <- out1 - out0
  close.top <- which(0<=blip&blip<=tau)
  close.bottom <- which(-tau<=blip&blip<0)
  out0[close.top] <- out1[close.top]-tau
  out1[close.bottom] <- out0[close.bottom]-tau
   
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
   ##seealso<< Qbar1

  ##references<<       Chambaz,        Zheng,       van        der       Laan,
  ##https://hal.archives-ouvertes.fr/hal-01301297.

   ##details<< This conditional expectation of \eqn{Y} given \eqn{(A,W)} where
   ##\eqn{W=(U,V)} is \deqn{\frac{1}{2} + \frac{3}{8}*\cos(\rho*V*\pi*U)} when
   ##\eqn{A=1}  and \deqn{\frac{1}{2}  + \frac{1}{4}*\sin(3*\rho/V*\pi*U)} when
   ##\eqn{A=0}.

   
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

   ##details<< This  conditional variance  of \eqn{Y} given  \eqn{(A,W)} where
   ##\eqn{W=(U,V)}     is    \deqn{(\rho+U+V)^2}     when    \eqn{A=1}     and
   ##\deqn{(\frac{1}{\rho+V}+U)^2} when \eqn{A=0}.

   
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

   ##details<< This  conditional expectation  of \eqn{Y} given  \eqn{(A,W)} is
   ##always equal to \eqn{\rho}.


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
### \code{tm.model}). The maximum value is \code{1-Gmin}.
 piV=c(1/2, 1/3, 1/6),
### Marginal distribution of \eqn{V}. Defaults to \code{c(1/2, 1/3, 1/6)}.
  family=c("beta", "gamma"),
### A \code{character}, either "beta" (default)  or "gamma", the nature of the
### law of outcome.
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

  ## Argument 'piV'
  piV <- Arguments$getNumerics(piV, range=c(0,1))
  if (sum(piV)!=1) {
    throw("Argument 'piV' should consist of non-negative weights summing to one.") 
  }

  ## Argument 'family':
  family <- match.arg(family)
   
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

   obs <- getSample(n, tm=oneOne, piV=piV, Qbar=Qbar, Vbar=Vbar, family=family)
   ## ## not  'what="ATE"' because otherwise no data is output!
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

getOptVar <- function#Computes the Optimal Variance Given a Parametric Model When Targeting the Average Treatment Effect
### Computes the optimal variance when  targeting the Average Treatment Effect
### (ATE). The optimal variance is defined as the variance associated with the
### optimal treatment mechanism  within a given parametric  model of treatment
### mechanisms. The computation is based on Monte Carlo simulation.
(tm.model=formula(A~1),
### A parametric  model \eqn{{\cal G}} of treatment  mechanisms. The procedure
### targets  the optimal treatment  mechanism within  this model.  Defaults to
### \code{formula(A~1)}.
 Gmin=1e-2,
### A  small positive  \code{numeric}, the  minimum value  of elements  of the
### parametric  model \eqn{{\cal  G}}  of treatment  mechanisms (see  argument
### \code{tm.model}). The maximum value is \code{1-Gmin}.
 piV=c(1/2, 1/3, 1/6),
### Marginal distribution of \eqn{V}. Defaults to \code{c(1/2, 1/3, 1/6)}.
  family=c("beta", "gamma"),
### A \code{character}, either "beta" (default)  or "gamma", the nature of the
### law of outcome.
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

  ## Argument 'piV'
  piV <- Arguments$getNumerics(piV, range=c(0,1))
  if (sum(piV)!=1) {
    throw("Argument 'piV' should consist of non-negative weights summing to one.") 
  }

  ## Argument 'family':
  family <- match.arg(family)

   
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
                    piV=piV,  family=family, Qbar=Qbar, Vbar=Vbar,
                    n=n, verbose=verbose)
  optVar <- getSample(n=n, tm=optTm,
                      piV=piV,  family=family, Qbar=Qbar, Vbar=Vbar,
                      what="ATE")[2]^2
  names(optVar) <- NULL
  
  return(optVar)
### Returns  the  optimal  variance   associated  with  the  parametric  model
### \eqn{{\cal G}} when targeting the Average Treatment Effect (ATE).
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

makeLearnQ <- function#Builds a Parametric Working Model Based on Sample Size
### Builds a Parametric Working Model Based on Sample Size
() {
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ##details<< This functions builds a sample-size-dependent parametric working
  ##model.   Two fine-tune  parameters  are determined  based  on sample  size
  ##\eqn{n}:   \eqn{deg=3  +   \lfloor   n/500\rfloor}  and   \eqn{nlev=\lceil
  ##n/250\rceil}.  If  \eqn{nlev} equals one, then  the model is given  by the
  ##\code{formula}       \deqn{Y~I(A=0)*V*poly(U)+I(A=1)*V*poly(U)}      where
  ##\eqn{poly(U)} consists  of \eqn{deg} orthogonal  polynoms of degrees  1 to
  ##\eqn{deg}. If \eqn{nlev} is larger than or equal to two, then the model is
  ##given                 by                the                 \code{formula}
  ##\deqn{Y~I(A=0)*V*(poly(U)+step(U))+I(A=1)*V*(poly(U)+step(U))}       where
  ##\eqn{poly(U)} consists  of \eqn{deg} orthogonal  polynoms of degrees  1 to
  ##\eqn{deg} and  \eqn{step(U)} consits of \eqn{nlev}  indicator functions of
  ##subsets of \eqn{[0,1]} of the form \eqn{\{x:x<\frac{k}{nlev}\}}.
  
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
### Builds a Parametric Working Model Based on Sample Size
() {
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ##details<< This functions builds a sample-size-dependent parametric working
  ##model.   Two fine-tune  parameters  are determined  based  on sample  size
  ##\eqn{n}: \eqn{deg=\min(\lceil  n/100\rceil, 6)}  and \eqn{nlev=\min(\lceil
  ##n/100\rceil, 10)}.  If  \eqn{nlev} equals one, then the model  is given by
  ##the   \code{formula}    \deqn{Y~I(A=0)*V*poly(U)+I(A=1)*V*poly(U)}   where
  ##\eqn{poly(U)} consists  of \eqn{deg} orthogonal  polynoms of degrees  1 to
  ##\eqn{deg}. If \eqn{nlev} is larger than or equal to two, then the model is
  ##given  by   the  \code{formula}   \deqn{Y~I(A=0)*V*(poly(U)+\sum_k  step_k
  ##(U))+I(A=1)*V*(poly(U)+\sum_k step_k (U))} where \eqn{poly(U)} consists of
  ##\eqn{deg} orthogonal  polynoms of degrees  1 to \eqn{deg}  and \eqn{step_k
  ##(U)} consits of \eqn{nlev} ortogonal polynoms of degrees 1 to \eqn{deg} in
  ##\eqn{U} times the  indicator function of the subset of  \eqn{[0,1]} of the
  ##form \eqn{\{x:\frac{k}{nlev}\le x <\frac{k}{nlev}\}}.

  
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

smoothIndicator <- function#Smooth Approximation of \eqn{1\{x>0\}} Over \eqn{[-1,+1]}
### Smooth approximation of \eqn{1\{x>0\}} over  \eqn{[-1,+1]}, used to map an
### estimated blip function to a stochastic treatment mechanism.
(xx,
### A \code{vector} of \code{numerics} between -1 and +1.
  exploit,
### A small positive  \code{numeric}.
  explore
### A small positive  \code{numeric}.
) {
  ##references<<       Chambaz,        Zheng,       van        der       Laan,
  ##https://hal.archives-ouvertes.fr/hal-01301297.
 
  ##details<< This  function is  a non-decreasing, Lipschitz  approximation of
  ##\eqn{1\{x>0\}} over \eqn{[-1,+1]}.  It takes  its values in \eqn{[t, 1-t]}
  ##where \eqn{t}  equals \code{exploit},  equals \eqn{t}  if its  argument is
  ##smaller  than \eqn{-\xi}  and \eqn{1-t}  if  its argument  is larger  than
  ##\eqn{\xi}, and 0.5 if its argument equals 0.

  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 
  ## Argument 'xx':
  xx <- Arguments$getNumerics(xx, c(-1, 1))

  ## Argument 'exploit':
  exploit <- Arguments$getNumeric(exploit, c(0, 1/2))

  ## Argument 'explore':
  explore <- Arguments$getNumeric(explore, c(0, 1/2))

  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Core
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  aa <- -(1/2 - exploit)/(2*explore^3)
  bb <- (1/2 - exploit)/(2*explore/3)
  cc <- 1/2
  
  out <- rep(exploit, length(xx))
  pos <- (xx > explore)
  btwn <- (-explore <= xx & xx <= explore)
  out[pos] <- 1-exploit
  out[btwn] <- aa*xx[btwn]^3 + bb*xx[btwn] + cc

  return(out)
}

############################################################################
## HISTORY:
## 2016-09-16
## o Created.
############################################################################

