setConstructorS3("TSMLCARA", function(#Creates a TSMLCARA Object.
### Creates a TSMLCARA object.
                                     what=c("ATE", "MOR"),
### A  \code{character}  indicating the  parameter  of  interest to  estimate.
### Either "ATE" for the Average  Treatment Effect, the difference between the
### means under  '\eqn{do(A=1)}' and  '\eqn{do(A=0)}', or  "MOR" for  the Mean
### under the Optimal treatment Rule '\eqn{do(A=r(W))}'.
                                     obs=data.frame(matrix(nrow=0, ncol=4, dimnames=list(NULL, c("W", "G", "A", "Y")))),
### A  \code{data.frame}  of  observations,  as  produced  by  \code{function}
### \code{getSample}.
                                     tm.ref=oneOne,
### A \code{function}, the reference treatment mechanism of the \code{TSMLCARA}
### object. Defaults to 'oneOne', the balanced treatment mechanism. 
                                     Gmin=1e-2,
### A  small positive  \code{numeric}, with  default value  \code{1e-2}.  When
### \code{what}  equals "ATE",  it is  the minimum  value of  elements of  the
### parametric  model \eqn{{\cal  G}}  of treatment  mechanisms (see  argument
### \code{tm.model}).  The  maximum value is \code{1-Gmin}.   When \code{what}
### equals "MOR",  it is the minimum  value of the conditional  probability of
### \eqn{do(A=r_n(W))} given \eqn{W}.
                                     Gexpl=1e-2,
### A small positive \code{numeric}, with default value \code{1e-2}, only used
### when  \code{what}  equals  "MOR",  in   which  case  it  lower-bounds  the
### conditional probability of \eqn{do(A=1-r_n(W))} given \eqn{W}.
                                     threxpl=1e-2,
### Either a small positive \code{numeric}, with default value \code{1e-2}, or
### a  function of  sample  size giving  such small  numbers,  only used  when
### \code{what}  equals  "MOR".    If  \eqn{0\in[Q_n-\theta,Q_n+\theta]}  with
### \eqn{\theta} equal  to \code{threxpl} then \eqn{r_n(W)}  is the proportion
### of the interval which lies above  0, thresholded at levels \code{Gmin} and
### \code{1-Gmin}.
                                     Qmin=0,
### A  small positive  \code{numeric}, the  minimum value  of  scaled outcomes
### \eqn{Y}. The maximum value is \code{1-Qmin}.
                                     flavor=c("parametric", "lasso"), 
### A  \code{character}  indicating  the   flavor  of  the  procedure,  either
### 'parametric' or 'lasso'.
                                     learnQ,
### A  model  \eqn{{\cal Q}}  of  conditional  expectations of  \eqn{Y}  given
### \eqn{(A,W)}  for  both  flavors  'parametric'  and  'lasso',  given  as  a
### \code{formula}  or  a  \code{function}  outputing  formulas.  Defaults  to
### \code{formula(Y~1)} for flavors  'parametric' and 'lasso'.
                                     tm.model=formula(A~1),
### A parametric model \eqn{{\cal G}}  of treatment mechanisms, only used when
### \code{what} equals "ATE". The procedure targets the optimal treatment
### mechanism within this model.  Defaults to \code{formula(A~1)}.
                                     tm.control=glm.control(maxit=500), 
### A \code{list}  of options to pass  to \code{glm} for the  targeting of the
### treatment mechanism within the model  defined by argument 'tm.model', only
### used     when      \code{what}     equals     "ATE".       Defaults     to
### \code{glm.control(maxit=500)}.
                                     conf.level=0.95,
###  A  \code{numeric},  the  confidence  level of  the  resulting  confidence
###  interval.
                                     ...,
### Additional parameters.
                                     verbose=FALSE
### A \code{logical}  or an \code{integer}  indicating the level  of verbosity
### (defaults to 'FALSE').
                                     ) {
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ## Argument 'what':
  what <- match.arg(what)
  
  ## Argument 'obs':
  obs <- validateArgumentObs(obs)

  ## Argument 'tm.ref'
  mode <- mode(tm.ref)
  if (mode != "function") {
    throw("Argument 'tm.ref' should be of mode 'function', not '", mode)
  }
  
  ## Argument 'Gmin'
  Gmin <- Arguments$getNumeric(Gmin, c(0, 1/2))

  ## Argument 'Gexpl'
  Gexpl <- Arguments$getNumeric(Gexpl, c(0, 1/2))

  ## Argument 'threxpl'
  mode <- mode(threxpl)
  if (!(mode %in% c("numeric", "function"))) {
    throw("Argument 'threxpl' should be of mode either 'numeric' or 'function', not ", mode) 
  } else if (mode=="numeric") {
    threxpl <- Arguments$getNumeric(threxpl, c(0, 1/2))
  } 
  
  ## Argument 'Qmin'
  Qmin <- Arguments$getNumeric(Qmin, c(0, 1/4))

  ## Arguments 'flavor' and 'learnQ'
  flavor <- match.arg(flavor)
  if (missing(learnQ)) {
    learnQ <- switch(flavor,
                     parametric=formula(Y~1),
                     lasso=formula(Y~1))
  } else {
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
  }
    
  ## Argument 'tm.model'
  if (!identical(class(tm.model), "formula") | tm.model[[2]]!="A") {
    throw("Argument 'tm.model' should be a formula with response 'A', not:",
          deparse(substitute(form, list(form=tm.model)))) 
  }

  
  ## Argument 'conf.level'
  conf.level <- Arguments$getNumeric(conf.level, c(0, 1))

  ## Argument 'verbose'
  verbose <- Arguments$getVerbose(verbose)
  verbose <- less(verbose, 10)
  
  ## Arguments '...':
  args <- list(...);
  if (length(args) > 0) {
    argsStr <- paste(names(args), collapse=", ")
    throw("Unknown arguments: ", argsStr)
  } 

  ##
  ## preparing 'targetLink'
  ##

  if (what=="ATE") {
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
  } else if (what=="MOR") {
    ## link
    linkfun <- function(mu) {# maps 'mu' (probs) to 'eta' (reals)
      logit((mu-Qmin)/(1-2*Qmin))
    }
    ## inverse link
    linkinv <- function(eta) {# maps 'eta' (reals) to 'mu' (probs)
      Qmin + (1-2*Qmin)*expit(eta)
    }
    ## derivative of inverse link wrt eta
    mu.eta <- function(eta) {
      expit.eta <- expit(eta)
      (1-2*Qmin)*expit.eta*(1-expit.eta)
    }
    ## test of validity for eta 
    valideta <- function(eta) {
      TRUE
    }
    ## "deviance residuals" as a function of eta
    dev.resids <- binomial()$dev.resids
    ## computing 'mustart'
    initialize <- binomial()$initialize
  }
  link <- "targetLink"
  targetLink <- structure(list(linkfun=linkfun, linkinv=linkinv,
                               mu.eta=mu.eta, valideta=valideta, 
                               dev.resids=dev.resids,
                               initialize=initialize, name=link),
                          class = "link-glm")
  
  
  ## setting 'G'
  Gtab <- extractG(obs)

  if (what=="ATE") {
    nms <- c("nobs", "psi", "psi.sd")
  } else if (what=="MOR") {
    nms <- c("nobs", "psi", "psi.sd", "regret", "regret.sd")
  }
  hist <- matrix(c(nrow(obs), rep(NA, length(nms)-1)), 1, length(nms))
  rownames(hist) <- "create"
  colnames(hist) <- nms
  history <- list(hist=hist, Gstar.model=list(), learnedQ=list())
  class(history) <- ".history"

  if (what=="ATE") {
    extend(Object(), "TSMLCARA",
           .what=what,
           .obs=obs,
           .tm.ref=tm.ref,
           .flavor=flavor,
           .learnQ=learnQ,
           .tm.model=tm.model,
           .tm.control=tm.control,
           .targetLink=targetLink,
           .weights=NULL,
           .Gmin=Gmin,
           .Gexpl=Gexpl,
           .threxpl=threxpl,
           .Qmin=Qmin,
           .Gstar=oneOne,
           .learnedQ=NULL,
           .Gtab=Gtab, .Qtab=matrix(NA, 0, 3), .QEsptab=matrix(NA, 0, 3),
           .psi=NA, .psi.sd=NA,
           .history=history, .step=c(update=0, target=0),
           .conf.level=conf.level
           )
  } else if (what=="MOR") {
    extend(Object(), "TSMLCARA",
           .what=what,
           .obs=obs,
           .tm.ref=tm.ref,
           .flavor=flavor,
           .learnQ=learnQ,
           .tm.model=tm.model,
           .tm.control=tm.control,
           .targetLink=targetLink,
           .weights=NULL,
           .Gmin=Gmin,
           .Gexpl=Gexpl,
           .threxpl=threxpl,
           .Qmin=Qmin,
           .Gstar=oneOne,
           .learnedQ=NULL,
           .Gtab=Gtab, .Qtab=matrix(NA, 0, 3), .QEsptab=matrix(NA, 0, 3),
           .psi=NA, .psi.sd=NA,
           .regret=NA, .regret.sd=NA,
           .history=history, .step=c(update=0, target=0),
           .conf.level=conf.level
           )
  }
})

setMethodS3("getWhat", "TSMLCARA", function(#Retrieves What Is the Parameter of Interest
### Retrieves the \code{character} indicating what is the parameter of interest.
    this,
### An object of class \code{TSMLCARA}.
    ...
### Not used.
    ) {
  ##alias<< getWhat
  ##seealso<< tsml.cara.rct

  this$.what
### The parameter of interest.
})


setMethodS3("getObs", "TSMLCARA", function(#Retrieves the Current Data Set
### Retrieves the \code{data.frame} of observations currently involved in the procedure.
    this,
### An object of class \code{TSMLCARA}.
    ...
### Not used.
    ) {
  ##alias<< getObs
  ##seealso<< tsml.cara.rct

  this$.obs
### The data set involved in  the procedure.
})

setMethodS3("addNewSample", "TSMLCARA", function(#Adds the Newly Sampled Observations
### Adds the newly sampled \code{data.frame} of observations and updates 
    this,
### An object of class \code{TSMLCARA}.
    newobs,
### The  newly  sampled  \code{data.frame}  of  observations  to  add  to  the
### \code{TSMLCARA} object.
    ...
### Not used.
    ) {
  ##alias<<   addNewSample
  ##seealso<<    tsml.cara.rct
  ##details<<  The  output  \code{data.frame}  has a  "from"  attribute  which
  ##indicates the  row in the  concatenated data set  of the first row  of the
  ##added data set.
  
  ## Argument 'newobs':
  newobs <- validateArgumentObs(newobs)

  ## retrieving 'obs'
  obs <- getObs(this)
  nobs <- nrow(obs)

  ## updating 'obs' 
  newobs <- rbind(obs, newobs)
  rownames(newobs) <- 1:nrow(newobs)
  this$.obs <- newobs
})


setMethodS3("getTm.ref", "TSMLCARA", function(this, ...) {
  this$.tm.ref
})

setMethodS3("setTm.ref", "TSMLCARA", function(#Sets Reference Treatment Mechanism.
### Sets the reference treatment mechanism of a \code{TSMLCARA} object.
    this,
### An object of class \code{TSMLCARA}.
    tm=oneOne,
### A \code{function} of \eqn{W}, a treatment mechanism, with \code{oneOne} as
### default value.
    ...
### Not used.
    ) {
  ##alias<< setTm.ref
  ##seealso<< as.character
  
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  ## Argument 'tm':
  mode <- mode(tm)
  if (mode != "function") {
    throw("Argument 'tm' should be of mode 'function', not '", mode)
  }
    
  this$.tm.ref <- tm
})

setMethodS3("getTm.model", "TSMLCARA", function(this, ...) {
  this$.tm.model
})

setMethodS3("setTm.model", "TSMLCARA", function(#Sets Model of  Treatment Mechanisms.
### Sets the parametric model of treatment mechanisms of a \code{TSMLCARA} object.
    this,
### An object of class \code{TSMLCARA}.
    tm.model=formula(A~1),
### A \code{formula}  (default \code{formula(A~1)}) describing  the parametric
### model. The formula must be of the form 'A~predictors'.
    ...
### Not used.
    ) {
  ##alias<< setTm.model
  ##seealso<< as.character
  
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  ## Argument 'tm.model':
  if (!identical(class(tm.model), "formula") | tm.model[[2]]!="A") {
    throw("Argument 'tm.model' should be a formula with response 'A', not:",
          deparse(substitute(form, list(form=tm.model)))) 
  }

  this$.tm.model <- tm.model
})

setMethodS3("getTm.control", "TSMLCARA", function(this, ...) {
  this$.tm.control
})

setMethodS3("setTm.control",  "TSMLCARA", function(#Sets Control  Arguments
### Sets the 'control' arguments of a \code{TSMLCARA} object 
    this,
### An object of class \code{TSMLCARA}.
    tm.control=glm.control(maxit=500),
### A \code{list} of control  arguments produced by \code{glm.control} for the
### targeting of the treatment mechanism  within the model defined by argument
### 'tm.model' of the object.  Defaults to \code{glm.control(maxit=500)}.
    ...
### Not used.
    ) {
  ##alias<< setTm.control
  ##seealso<< as.character
  
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  ## Argument 'tm.control':

  this$.tm.control <- tm.control
})

setMethodS3("getTargetLink", "TSMLCARA", function(this, ...) {
  this$.targetLink
})



setMethodS3("getGstar", "TSMLCARA", function(this, ...) {
  this$.Gstar
})

setMethodS3("getQ", "TSMLCARA", function(this, ...) {
  get("Q", envir=environment(this$.Gstar))
})


setMethodS3("setGstar", "TSMLCARA", function(#Sets the Targeted Treatment Mechanism
### Sets the targeted treatment mechanism of a \code{TSMLCARA} object.
    this,
### An object of class \code{TSMLCARA}.
    tm,
### A \code{function}  of \eqn{W}, a targeted treatment  mechanism as obtained
### by applying \code{estimateGstar}.
    ...
### Not used.
    ) {
  ##alias<< setGstar
  ##seealso<< as.character, estimateGstar
  
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  ## Argument 'tm':
  mode <- mode(tm)
  if (mode != "function") {
    throw("Argument 'tm' should be of mode 'function', not '", mode)
  }
    
  this$.Gstar <- tm
})


setMethodS3("getGmin", "TSMLCARA", function(this, ...) {
  this$.Gmin
})

setMethodS3("setGmin", "TSMLCARA", function(#Sets Value of 'Gmin'.
### Sets the value of 'Gmin' of a \code{TSMLCARA} object.
    this,
### An object of class \code{TSMLCARA}.
    Gmin,
### A \code{numeric},  the minimal value  of elements of the  parametric model
### \eqn{{\cal   G}}   of  treatment   mechanisms.   The   maximum  value   is
### \code{1-Gmin}.
    ...
### Not used.
    ) {
  ##alias<< setGmin
  ##seealso<< as.character
  
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  ## Argument 'Gmin':
  Gmin <- Arguments$getNumeric(Gmin, range=c(0, 1/2))
  
  this$.Gmin <- Gmin
})

setMethodS3("getGexpl", "TSMLCARA", function(this, ...) {
  this$.Gexpl
})

setMethodS3("setGexpl", "TSMLCARA", function(#Sets Value of 'Gexpl'.
### Sets the value of 'Gexpl' of a \code{TSMLCARA} object.
    this,
### An object of class \code{TSMLCARA}.
    Gexpl,
### A  \code{numeric},  the  minimal  value  the  conditional  probability  of
### \eqn{do(A=1-r_n(W))} given \eqn{W}.
    ...
### Not used.
    ) {
  ##alias<< setGexpl
  ##seealso<< as.character
  
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  ## Argument 'Gexpl':
  Gexpl <- Arguments$getNumeric(Gexpl, range=c(0, 1/2))
  
  this$.Gexpl <- Gexpl
})

setMethodS3("getThrexpl", "TSMLCARA", function(this, ...) {
  this$.threxpl
})

setMethodS3("setThrexpl", "TSMLCARA", function(#Sets Value of 'threxpl'.
### Sets the value of 'threxpl' of a \code{TSMLCARA} object.
    this,
### An object of class \code{TSMLCARA}.
    threxpl=1e-2,
### Either a small positive \code{numeric}, with default value \code{1e-2}, or
### a  \code{function} of  sample size  giving  such numbers,  only used  when
### \code{what}  equals 'MOR'.  If  \eqn{0\in[Q_n-\theta,Q_n+\theta]} with
### \eqn{\theta} equal  to \code{threxpl} then \eqn{r_n(W)}  is the proportion
### of the interval which lies  above 0, thresholded at levels \code{Gmin} and
### \code{1-Gmin}.
   ...
### Not used.
    ) {
  ##alias<< setThrexpl
  ##seealso<< as.character
  
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  ## Argument 'threxpl'
  mode <- mode(threxpl)
  if (!(mode %in% c("numeric", "function"))) {
    throw("Argument 'threxpl' should be of mode either 'numeric' or 'function', not ", mode) 
  } else if (mode=="numeric") {
    threxpl <- Arguments$getNumeric(threxpl, c(0, 1/2))
  } 
  
  this$.threxpl <- threxpl
})



setMethodS3("getQmin", "TSMLCARA", function(this, ...) {
  this$.Qmin
})

setMethodS3("setQmin", "TSMLCARA", function(#Sets Value of 'Qmin'.
### Sets the value of 'Qmin' of a \code{TSMLCARA} object.
    this,
### An object of class \code{TSMLCARA}.
    Qmin,
### A \code{numeric},  used to  rescale the values  of the outcome  \eqn{Y} by
### using the scaling \code{function} \code{scaleY}.
    ...
### Not used.
    ) {
  ##alias<< setQmin
  ##seealso<< scaleY
  
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  ## Argument 'Qmin':
  Qmin <- Arguments$getNumeric(Qmin, range=c(0, 1/4))
  
  this$.Qmin <- Qmin
})


setMethodS3("getConf.level", "TSMLCARA", function(this, ...) {
  this$.conf.level
})

setMethodS3("setConfLevel", "TSMLCARA", function(#Sets a Confidence Level
### Sets the confidence level of a \code{TSMLCARA} object.
                                                  this,
### An object of class \code{TSMLCARA}.
                                                  confLevel,
### A \code{numeric}, confidence interval level.
                                                  ...
### Not used.
                                                  ) {
  ##alias<< setConfLevel
  ##seealso<< as.character
  
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  ## Argument 'confLevel':
  conf.level <- Arguments$getNumeric(confLevel, range=c(0, 1))
  
  this$.conf.level <- conf.level
})

setMethodS3("as.character", "TSMLCARA", function(#Returns a Description of a TSMLCARA Object
### Returns a short string describing the TSMLCARA object.
                                                 x,
### An object of class \code{TSMLCARA}.
                                                 ...
### Not used.
                                                 ) {
  ##alias<< as.character
  ##seealso<< tsml.cara.rct, setConfLevel
  
  this <- x  ## To please R CMD check
  s <- sprintf("%s object:", class(this)[1])
  s <- c(s, "")
  ## what
  if (getWhat(this)=="ATE") {
    what <- "risk difference"
  } else if (getWhat(this)=="MOR") {
    what <- "optimal treatment rule"
  }
  s <- c(s, "Parameter of interest: ", what)
  s <- c(s, "")

  ## sample size
  s <- c(s, sprintf("Sample size: %s", nrow(getObs(this))))
  s <- c(s, "")

  ## updates and fluctuations
  step <- getStep(this)
  s <- c(s, sprintf("Number of updating steps:\t%s", step["update"]))
  s <- c(s, sprintf("Number of targeting steps:\t%s", step["target"]))
  s <- c(s, "")
  
  ## psi
  s <- c(s, sprintf("Estimator of psi:\t\t%s", signif(getPsi(this), 3)))

  ## std error
  psi.sd <- getPsiSd(this)
  s <- c(s, sprintf("Estimated standard error:\t%s", signif(psi.sd, 3)))
  
  ## confidence intervals
  psi <- getPsi(this)
  alpha <- 1-getConf.level(this)
  n <- nrow(getObs(this))
  CI <- psi+c(-1, 1)*psi.sd*qnorm(1-alpha/2)/sqrt(n)
  CI <- signif(CI, 3)
  s <- c(s, sprintf("%s-confidence interval:\t[%s, %s]", 1-alpha, CI[1], CI[2]))

  class(s) <- "GenericSummary"
  s
### A character string summarizing the content of the object. The summary contains:
### \itemize{
### \item{The name of the parameter of interest.}
### \item{The sample size of the data set involved in the procedure.}
### \item{The value of the TSMLE and its estimated standard error.}
### \item{A confidence interval with default level of 95% (the level can be changed by using \code{\link{setConfLevel}}).}
### }
}, private=TRUE)

setMethodS3("getPsi", "TSMLCARA", function(#Returns the Current Estimator
### Returns the current value of the estimator.
                                           this,
### An object of class \code{TSMLCARA}.
                                           ...
### Not used.
                                           ) {
  ##alias<< getPsi
  ##seealso<< tsml.lcara.rct, getHistory, getPsiSd
  this$.psi
### Retrieves  the current value  of the estimator. 
})

setMethodS3("getPsiSd", "TSMLCARA", function(#Returns the Current Estimated Standard Deviation of the Current Estimator 
### Returns  the current  value of  the  estimated standard  deviation of  the
### current estimator.
    this,
### An object of class \code{TSMLCARA}.
    ...
### Not used.
    ) {
  ##alias<< getPsiSd
  ##seealso<< tsml.lcara.rct, getHistory, getPsi
  this$.psi.sd
  ### Retrieves the  estimated standard deviation of the  current estimator of
  ###  the parameter  of interest.
})

setMethodS3("getRegret", "TSMLCARA", function(#Returns the Current Estimator of the  Empirical Regret
### Returns the current value of the estimator of the empirical regret when 'what' equals 'MOR'. 
    this,
### An object of class \code{TSMLCARA}.
    ...
### Not used.
    ) {
  ##alias<< getRegret
  ##seealso<< tsml.lcara.rct, getHistory, getRegretSd
  what <- this$.what
  if (what!="MOR") {
    throw("Function 'getRegret' can be used when 'what' equals 'MOR', not ", what)
  }
  this$.regret
  ### The current estimate of the empirical regret. 
})

setMethodS3("getRegretSd", "TSMLCARA", function(#Returns the Current Estimated Standard Deviation of the Current Estimator of the Empirical Regret
### Returns  the current  value of  the  estimated standard  deviation of  the
### current estimator of the empirical regret when 'what' equals 'MOR'.
    this,
### An object of class \code{TSMLCARA}.
    ...
### Not used.
    ) {
  ##alias<< getRegretSd
  ##seealso<< tsml.lcara.rct, getHistory, getRegret
  what <- this$.what
  if (what!="MOR") {
    throw("Function 'getRegretSd' can be used when 'what' equals 'MOR', not ", what)
  }
  this$.regret.sd
  ### The  estimated  standard  deviation  of  the current  estimator  of  the
  ### empirical regret.
})



setMethodS3("getGtab", "TSMLCARA", function(this, ...) {
    this$.Gtab
})

setMethodS3("setGtab", "TSMLCARA", function(this, Gtab, ...) {
  ## Argument 'Gtab':
  Gtab <- Arguments$getNumerics(Gtab, c(0, 1))
  
  ## retrieving 'obs'
  obs <- getObs(this)
  
  if (length(Gtab) != nrow(obs)) {
    throw("Vector 'Gtab' has fewer lines than the number of rows of 'obs'...")
  }
  this$.Gtab <- Gtab
})


setMethodS3("getQtab", "TSMLCARA", function(this, ...) {
    this$.Qtab
})

setMethodS3("setQtab", "TSMLCARA", function(this, Qtab, ...) {
  ## Argument 'Qtab':
  Qtab <- Arguments$getNumerics(Qtab)
  
  ## retrieving 'obs'
  obs <- getObs(this)
  
  if (nrow(Qtab) != nrow(obs)) {
    throw("Matrix 'Qtab' has fewer rows than 'obs'...")
  }
  this$.Qtab <- Qtab
})

setMethodS3("getQEpstab", "TSMLCARA", function(this, ...) {
    this$.QEpstab
})

setMethodS3("setQEpstab", "TSMLCARA", function(this, Qtab, ...) {
  ## Argument 'Qtab':
  Qtab <- Arguments$getNumerics(Qtab)
  
  ## retrieving 'obs'
  obs <- getObs(this)
  
  if (nrow(Qtab) != nrow(obs)) {
    throw("Matrix 'Qtab' has fewer rows than 'obs'...")
  }
  this$.QEpstab <- Qtab
})

setMethodS3("getLearnedQ", "TSMLCARA", function(this, ...) {
    this$.learnedQ
})

setMethodS3("setLearnedQ", "TSMLCARA", function(this, learnedQ, ...) {
  this$.learnedQ <- learnedQ
})

setMethodS3("getLearnQ", "TSMLCARA", function(this, ...) {
    this$.learnQ
})

setMethodS3("setLearnQ", "TSMLCARA", function(this, learnQ, ...) {
  this$.learnQ <- learnQ
})


printHistory <- function(#Prints a Summary of an History of a TSMLCARA Object
### Prints a summary of an history of a TSMLCARA object.
                          x,
### The \code{history} of a TSMLCARA object.)
                          ...)
### Not used.
{
  ##seealso<< getHistory
  cat("$hist:")
  print(x$hist)
  cat("\n$Gstar.model: a list with", length(x$Gstar.model), "entries\n\n")
  cat("$learnedQ: a list with", length(x$learnedQ), "entries\n")
  invisible(x)
}

setMethodS3("getHistory", "TSMLCARA", function(#Retrieves the History of a TSMLCARA Object
### Retrieves the \code{history} of a TSMLCARA object.
                                              this,
### A \code{TSMLCARA} object.
                                              ...
### Not used.
                                              ) {
  ##alias<< getHistory
  ##seealso<< as.character
  
  this$.history
### The history of the given TSMLCARA object.
})

setMethodS3("updateHistory", "TSMLCARA", function(#Updates History of a TSMLECARA Object.
### Updates the "history" of a \code{TSMLECARA} object.
                                                 this,
### A \code{TSMLCARA} object.
                                                 when=c("update", "target"),
### A \code{character} indicating  at which step of the  procedure the history
### updating takes place. Either "update" or "target".
                                                 ...
### Not used.
                                                 ) {
  
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  when <- match.arg(when)
  what <- getWhat(this)
  history <- tsml.cara.rct::getHistory.TSMLCARA(this)
  step <- getStep(this)[when]
  nm <- paste(when, as.character(step), sep="")
  hh <- history$hist
  rownames.hh <- rownames(hh)
  if (when=="update") {
    if (what=="ATE") {
      cols <- c("psi", "psi.sd")
    } else if (what=="MOR") {
      cols <- c("psi", "psi.sd", "regret", "regret.sd")
    }
    hh <- rbind(hh, c(nobs=nrow(getObs(this)), hh[nrow(hh), cols]))
    rownames(hh) <- c(rownames.hh, nm) 
    history$hist <- hh
    Gstar <- getGstar(this)
    learnedQ <- getLearnedQ(this)
    history$Gstar.model[[nm]] <- attr(Gstar, "model")
    history$learnedQ[[nm]] <- learnedQ
  } else {
    if (what=="ATE") {
      hh <- rbind(hh, c(nobs=nrow(getObs(this)), psi=getPsi(this), psi.sd=getPsiSd(this)))
    } else if (what=="MOR") {
      hh <- rbind(hh, c(nobs=nrow(getObs(this)), psi=getPsi(this), psi.sd=getPsiSd(this),
                        regret=getRegret(this), regret.sd=getRegretSd(this)))      
    }
    rownames(hh) <- c(rownames.hh, nm)
    history$hist <- hh
  }
  this$.history <- history
})



setMethodS3("getFlavor", "TSMLCARA", function(this, ...) {
    this$.flavor
})

setMethodS3("getStep", "TSMLCARA", function(this, name, ...) {
  this$.step
})

setMethodS3("plot", "TSMLCARA", function#Plots a TSMLCARA Object
### Plots a  TSMLCARA object.
(x,
### An object of class \code{TSMLCARA}.
 truth=NULL,
### Defaults  to \code{NULL}.   When estimating  the Average  Treatment Effect
### ('what' equal to "ATE"), can be set to a \code{vector}, the true values of
### both  the  targeted parameter  and  standard  deviation of  the  efficient
### influence  curve under  the optimal  treatment mechanism  given a  working
### model,  as output  by 'getOptVar'.   When  estimating the  Mean under  the
### Optimal treatment  Rule ('what' equal  to 'MOR'), can  be set to  either a
### \code{vector} or  a list of two  vectors.  In the former  case, the vector
### gives  the  true  values  of  both the  targeted  parameter  and  standard
### deviation of  the efficient  influence curve  under the  optimal treatment
### rule. In the latter case, the list  contains the same vector as before and
### a  second  vector   giving  the  values  of   the  targeted  data-adaptive
### parameters.
 regret=NULL,
### Defauls to \code{NULL}. When estimating the optimal treatment rule ('what'
### equal to 'MOR'), can  be set to a vector or a list  of two vectors. In the
### former case,  the vectors gives  the true values  of the empirical  or the
### counterfactual regrets. In the latter case, the list contains both.
  regret.mirror=FALSE,
### Should the regret(s) be mirrored? If  \code{regret} is given, should it be
### multiplied by (-1)? (defaults to \code{FALSE}).
  lower.bound=TRUE,
### Defaults to \code{TRUE}. Should the  confidence lower bound on the regrets
### be printed or a confidence interval?
  xlog=TRUE,
### Defaults  to \code{TRUE}.  Should x-axis be on a log-scale?
  colors=c("tomato", "olivedrab4"),
### A length 2 vector of  colors. Defauls to \code{c("tomato", "olivedrab4")},
### colors used for plotting.
 file=NULL,
### A \code{character}.   Defaults to NULL,  and if a \code{character}  then a
### PDF plot is produced, with name "'file'.pdf".
    ...
### Not used.
    ) {
  ##alias<< plot
  ##seealso<< tsml.cara.rct

  this <- x;
      
  ## Argument 'file'  
  if (!is.null(file)) {
    file <- Arguments$getCharacters(file)
    file <- paste(file, ".pdf", sep="")
  }

  ## Retrieving 'what'
  what <- getWhat(this)
  what.expand <- switch(what,
                        ATE="Average Treatment Effect (ATE)",
                        MOR="Mean of the Optimal treatment Rule (MOR)")

  ## Argument 'regret'
  if (what!="MOR" && !is.null(regret)) {
    throw("Regrets are only handled when 'what' equals 'MOR', not: ", what)
  }

  ## Argument 'regret.mirror'
  regret.mirror <- Arguments$getLogical(regret.mirror)
  
  
  ## Argument 'lower.bound'
  lower.bound <- Arguments$getLogical(lower.bound)

  ## Argument 'xlog'
  xlog <- Arguments$getLogical(xlog)

  ## Argument 'colors'
  colors <- Arguments$getCharacters(colors)
  
  ## Core
  history <- tsml.cara.rct::getHistory.TSMLCARA(this)
  Gstar.model <- history$Gstar.model
  history <- history$hist
  alpha <- 1-getConf.level(this)

  idx.target <- which(sapply(strsplit(rownames(history), "target"), length)==2)
  if (length(idx.target)==0) {
    throw("Nothing to plot since no fluctuation has been performed yet.")
  } else {
    hp <- history[idx.target, "psi"]
    hs <- history[idx.target, "psi.sd"]
    hn <- history[idx.target, "nobs"]
    
    ics <-  qnorm(alpha/2) * c(-1,1) %*% t(hs/sqrt(hn))
                                                        
    pch <- 20
    if (!is.null(truth)) {
      if (!is.list(truth)) {
        usd <- truth["psi"]+qnorm(alpha/2)*truth["psi.sd"]/sqrt(hn)
        lsd <- truth["psi"]-qnorm(alpha/2)*truth["psi.sd"]/sqrt(hn)
        if (what=="MOR") {
          uopt <- truth["psi"]+qnorm(alpha/2)*truth["psi.sd.opt"]/sqrt(hn)
          lopt <- truth["psi"]-qnorm(alpha/2)*truth["psi.sd.opt"]/sqrt(hn)
        } else {
          uopt <- usd
          lopt <- lsd
        }
        ylim <- range(ics+rep(1,2) %*% t(hp),
                      truth["psi"],
                      uopt, lopt, usd, lsd)
      } else {
        uopt <- truth[[1]]["psi"]+qnorm(alpha/2)*truth[[1]]["psi.sd.opt"]/sqrt(hn)
        lopt <- truth[[1]]["psi"]-qnorm(alpha/2)*truth[[1]]["psi.sd.opt"]/sqrt(hn)
        usd <- truth[[1]]["psi"]+qnorm(alpha/2)*truth[[1]]["psi.sd"]/sqrt(hn)
        lsd <- truth[[1]]["psi"]-qnorm(alpha/2)*truth[[1]]["psi.sd"]/sqrt(hn)        
        ylim <- range(ics+rep(1,2) %*% t(hp),
                      truth[[1]]["psi"],
                      uopt, lopt, usd, lsd,
                      truth[[2]][, "psin"])
      }
    } else {
      ylim <- range(ics+rep(1,2) %*% t(hp))
    }

    if (!is.null(file)) {
      pdf(file=file)
    }
    par(mfrow=c(2, 1), oma=c(0,0,0,0), mar=c(4,6,3,2)+.1)

    ##
    ## first plot
    ##
    plot(hn, hp, ylim=ylim, pch=pch, cex=2,
         xlab="Sample size at targeting steps",
         ylab=expression(paste(psi[n], "*")), log=ifelse(xlog, "x", ""))
    ## title(main=paste("Representing the '", deparse(substitute(this)), "' TSMLCARA object\n(option '", getWhat(this), "')", sep="")) 
    mtext(paste(sprintf("[%s-confidence interval]", 1-alpha)), 2, line=2)
    mtext(what.expand, 3, line=2, font=2)
    mtext("Targeted estimators", 3, line=1, font=4)
    dummy <- sapply(seq(along=hn), function(x) {lines(c(hn[x],hn[x]), hp[x]+ics[, x])})
    if (!is.null(truth)) {
      if (!is.list(truth)) {
        abline(a=truth["psi"], b=0, col=4, lwd=2)
      } else {
        abline(a=truth[[1]]["psi"], b=0, col=4, lwd=2)
        points(hn, truth[[2]][, "psin"], col=colors[1], pch=4, cex=3, lwd=2)
      }
    }
    lines(hn, uopt, col="grey48", lwd=2)
    lines(hn, lopt, col="grey48", lwd=2)
    lines(hn, usd, col="grey64", lwd=2)
    lines(hn, lsd, col="grey64", lwd=2)
    
    if (what=="ATE") {
      idx.update <- which(sapply(strsplit(rownames(history), "update"), length)==2)
      hn <- history[idx.update, "nobs"]
      nms <- names(Gstar.model[[1]])
      Gstar.param <- matrix((unlist(Gstar.model)), nrow=length(idx.update), byrow=TRUE)
      ##
      ## second plot
      ##
      matplot(hn, Gstar.param,
              xlab="Sample size at updating steps",
              ylab="", type="b", pch=pch, log=ifelse(xlog, "x", ""))
      legend("left", legend=nms, text.col=1:4)
      mtext("Targeted treatment mechanisms", 3, line=1, font=4)
    } else if (what=="MOR" & !is.null(truth)) {
      ##
      ## second plot
      ##
      hr <- history[idx.target, "regret"]
      hrs <- history[idx.target, "regret.sd"]
      if (regret.mirror) {
        hr <- -hr
        if (!is.null(regret)) {
          if (!is.list(regret)) {
            regret <- -regret
          } else {
            regret[[1]] <- -regret[[1]]
            regret[[2]] <- -regret[[2]]
          }
        }
      }
          
      if (lower.bound) {
        ics <-  qnorm(alpha) * c(1,1) %*% t(hrs/sqrt(hn))
      } else {
        ics <-  qnorm(alpha/2) * c(-1,1) %*% t(hrs/sqrt(hn))
      }
                                                        
      pch <- 20
      if (!is.null(regret)) {
        if (!is.list(regret)) {
          ylim <- range(c(ics+rep(1,2) %*% t(hr), regret), 0)
        } else {
          ylim <- range(c(ics+rep(1,2) %*% t(hr),
                          regret[[1]],
                          regret[[2]]), 0)
        }
      } else {
        ylim <- range(ics+rep(1,2) %*% t(hr), 0)
      }
      if (!regret.mirror) { ## regret defined as: n^-1 sum Y_i - Y_i(r(W_i))
        plot(hn, hr, ylim=ylim, pch=pch, cex=2,
             xlab="Sample size at updating steps",
             ylab=expression(paste(psi[n], "*-",
                                   frac(1,n),
                                   sum(Y[i], i==1, n))),
             log=ifelse(xlog, "x", ""))
      } else { ## regret defined as: n^-1 sum Y_i(r(W_i)) - Y_i
        plot(hn, hr, ylim=ylim, pch=pch, cex=2,
             xlab="Sample size at updating steps",
             ylab=expression(paste(psi[n], "*-",
                                   frac(1,n),
                                   sum(Y[i], i==1, n))),
             log=ifelse(xlog, "x", ""))        
      }
      if (!is.null(regret)) {
        if (is.list(regret)) {
          mtext("Regrets", 3, line=1, font=2)
        } else {
          mtext("Regret", 3, line=1, font=2)
        }
      }
      if (lower.bound) {
        dummy <- sapply(seq(along=hn), function(x) {
          lines(hn[x]+ c(-1,1)*20, rep(hr[x]+ics[1, x], 2))
          lines(c(hn[x],hn[x]), c(hr[x]+ics[1, x], max(ylim))) 
        })
        mtext(paste(sprintf("[%s-lower confidence-bound]", 1-alpha)), 2, line=2)
      } else {
        dummy <- sapply(seq(along=hn), function(x) {lines(c(hn[x],hn[x]), hr[x]+ics[, x])})
        mtext(paste(sprintf("[%s-confidence interval]", 1-alpha)), 2, line=2)
      }
      if (!is.null(regret)) {
        if (is.list(regret)) {
          points(hn, regret[[2]], col=colors[2], pch=1, lwd=2, cex=3)
          points(hn, regret[[1]], col=colors[1], pch=4, lwd=2, cex=3)
        } else {
          points(hn, regret, col=colors[2], pch=4, lwd=2, cex=3)
        }
      }
     }
    if(!is.null(file)) {
      dev.off()
    }
  }
}, private=TRUE)



############################################################################
## HISTORY:
## 2016-09-16
## o Created.
############################################################################

