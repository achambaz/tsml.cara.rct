setMethodS3("update", "TSMLCARA", function(#Updates a TSMLCARA Object
### Updates a TSMLCARA object.
                                          this,
### A \code{TSMLCARA} object, as created by \code{TSMLCARA}.
                                          flavor=c("parametric", "lasso"),  
### A \code{character}  indicating the  'flavor' of the  procedure. 
                                          ...,
### Additional parameters.
                                          verbose=FALSE
### A \code{logical}  or an \code{integer}  indicating the level  of verbosity
### (defaults to 'FALSE').
                                        ) {
  ##alias<<  update
  ##references<<  Chambaz, van  der Laan,  Zheng, Chapter  16, Modern Adaptive Randomized Clinical  Trials: Statistical, Operational, and Regulatory    Aspects,    by    A.    Sverdlov    (CRC    Press,    2015).
  ##seealso<<tsml.cara.rct, targetPsi
  
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ## Arguments 'flavor'
  if (missing(flavor)) {
    flavor <- getFlavor(this)
  }
  
 
  ## Argument 'verbose'
  verbose <- Arguments$getVerbose(verbose)
  verbose <- less(verbose, 0)

  ## retrieving 'what'
  what <- getWhat(this)
  
  ## retrieving 'obs', preparing 'A'
  obs <- getObs(this)
  A <- obs[, "A"]
                     
  ## retrieving 'Qmin'
  Qmin <- getQmin(this)
    
  ## retrieving 'tm.ref'
  tm.ref <- getTm.ref(this)

  ## retrieving 'learnQ'
  learnQ <- getLearnQ(this)

  if (what=="ATE") {
    ## retrieving 'tm.model' and 'tm.control'
    tm.model <- getTm.model(this)
    tm.control <- getTm.control(this)
    
    ## retrieving 'targetLink'
    targetLink <- getTargetLink(this)
  } 
  
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## learning
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
  
  verbose && ruler(verbose)
  verbose && enter(verbose, "UPDATING STEP")

  ## estimating 'Q'
  verbose && enter(verbose, "Estimating Q")
  
  W <- extractW(obs)
  GA <- blik(A, obs[, "G"])
  weights <- blik(A, tm.ref(W))
  weights <- weights/GA

  Q <- estimateQ(obs, learnQ=learnQ, weights=weights, Qmin=Qmin,
                 flavor=flavor, ..., verbose=verbose)

  
  
  Qtab <- matrix(NA, nrow=nrow(obs), ncol=3)
  colnames(Qtab) <- c("A", "A=0", "A=1")
  Qtab[, "A=1"] <- Q(1, W)
  Qtab[, "A=0"] <- Q(0, W)
  Aone <- (A==1)
  Qtab[Aone, "A"] <- Qtab[Aone, "A=1"]
  Qtab[!Aone, "A"] <- Qtab[!Aone, "A=0"]

  setQtab(this, Qtab)
  setLearnedQ(this, attr(Q, "model"))

  verbose && str(verbose, Qtab)
  verbose && exit(verbose)
  
  ## targeting 'Gstar'
  verbose && enter(verbose, "Targeting Gstar")

  if (what=="ATE") {
    weights <- (obs[, "Y"]-Qtab[, "A"])^2/GA
    
    Gstar <- targetGstar(obs=obs, weights=weights, tm.model=tm.model,
                         targetLink=targetLink, tm.control=tm.control,
                         ...,
                         verbose=verbose)
  } else if (what=="MOR") {
    Gstar <- targetOptRule(this, Q, ..., verbose=verbose)
  } 
  
  setGstar(this, Gstar)

  verbose && str(verbose, Gstar)
  verbose && exit(verbose)

  ## updating 'Gtab'
  verbose && enter(verbose, "Updating Gtab")

  setGtab(this, Gstar(W))
  
  verbose && str(verbose, Gstar)
  verbose && exit(verbose)
  
  ## updating 'history'
  verbose && enter(verbose, "Updating history")
  step <- getStep(this)
  step["update"] <- step["update"]+1
  this$.step <- step
  updateHistory(this, "update")
  verbose && str(verbose, tsml.cara.rct::getHistory.TSMLCARA(this))
  verbose && exit(verbose)

  verbose && exit(verbose)
  verbose && ruler(verbose)

})

setMethodS3("targetPsi", "TSMLCARA", function(#Targets a TSMLCARA Object Toward the Parameter Psi
### Targets a TSMLCARA object toward the parameter Psi.
                                        this,
### A \code{TSMLCARA} object, as created by \code{TSMLCARA}.
                                        ...,
### Additional parameters.
                                        verbose=FALSE
### A \code{logical}  or an \code{integer}  indicating the level  of verbosity
### (defaults to 'FALSE').
                                        ) {
  ##alias<< targetPsi

  ##references<<  Chambaz, van der Laan,  Zheng, Chapter 16, Modern Adaptive Randomized Clinical  Trials: Statistical, Operational, and Regulatory  Aspects,  by  A.    Sverdlov  (CRC  Press,  2015).

  ##seealso<< tsml.cara.rct, update
  
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ## Argument 'verbose'
  verbose <- Arguments$getVerbose(verbose)
  verbose <- less(verbose, 0)

  ## retrieving 'what'
  what <- getWhat(this)

  if (what=="ATE") {
    targetLink <- NULL
  }
  else if (what=="MOR") {
    targetLink <- getTargetLink(this)
  }
  
  ## retrieving 'obs', 'Gtab', 'Qtab', preparing 'Y', 'A', 'G' and 'GA'
  obs <- getObs(this)
  Y <- obs[, "Y"]
  A <- obs[, "A"]
  G <- obs[, "G"]
  GA <- blik(A, G)
  GstarW <- getGtab(this)
  Qtab <- getQtab(this)
  
  ## retrieving 'Qmin'
  Qmin <- getQmin(this)
  
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## targeting Psi
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && ruler(verbose)
  verbose && enter(verbose, "TARGETING PSI")

  W <- extractW(obs)
  GAstarW <- blik(A, GstarW)
  weights <- GAstarW/GA

  verbose && enter(verbose, "Fluctuating")
    
  eps <- fluctuateQ(obs=obs, Qtab=Qtab, GAstarW=GAstarW,
                    what=what, weights=weights, Qmin=Qmin,
                    targetLink=targetLink,
                    ..., verbose)
  verbose && str(verbose, eps)
  
  ## QEpstab  
  epsHtab <- matrix(NA, nrow=nrow(obs), ncol=3)
  colnames(epsHtab) <- colnames(Qtab)

  if (what=="ATE") {
    epsHtab[, "A=1"] <- eps/GstarW
    epsHtab[, "A=0"] <- -eps/(1-GstarW)
  } else if (what=="MOR") {
    rA <- ((Qtab[, "A=1"]-Qtab[, "A=0"])>0)
    epsHtab[, "A=1"] <- eps/GstarW
    epsHtab[!rA, "A=1"] <- 0
    epsHtab[, "A=0"] <- eps/(1-GstarW)
    epsHtab[rA, "A=0"] <- 0
  }
  Aone <- (A==1)
  epsHtab[Aone, "A"] <- epsHtab[Aone, "A=1"]
  epsHtab[!Aone, "A"] <- epsHtab[!Aone, "A=0"]
  
    
  scaledQtab <- scaleQmat(Qtab, wrt=Y, thr=Qmin)
  if (what=="ATE") {
    scaledQEpstab <- binomial()$linkinv(binomial()$linkfun(scaledQtab) + epsHtab) ##ie, expit(logit(...))
  } else if (what=="MOR") {
    scaledQEpstab <- targetLink$linkinv(targetLink$linkfun(scaledQtab) + epsHtab)
  }
  QEpstab <- scaleQmat(scaledQEpstab, reverse=TRUE)
  attr(QEpstab, "eps") <- eps
  
  setQEpstab(this, QEpstab)
  
  verbose && str(verbose, QEpstab)
  verbose && exit(verbose)
  
           
  ## estimating Psi
  verbose && enter(verbose, "estimating psi")
            
  ## CAUTION: same form in both cases! no need to fork
  
  if (what=="ATE") {
    ## epsHtab <- epsHtab/abs(eps) * matrix(weights, nrow=length(weights), ncol=3)
    epsHtab <- epsHtab/eps * matrix(weights, nrow=length(weights), ncol=3)
    psis <- estimatePsi(obs=obs, what=what, Qtab=NULL, QEpstab=QEpstab, epsHtab=epsHtab, verbose=verbose)
  } else if (what=="MOR") {
    epsHtab <- epsHtab/eps * matrix(weights, nrow=length(weights), ncol=3)
    psis <- estimatePsi(obs=obs, what=what, Qtab=Qtab, QEpstab=QEpstab, epsHtab=epsHtab, verbose=verbose)
  }
  
  this$.psi <- psis["psi"]
  
  this$.psi.sd <- psis["psi.sd"]
  
  verbose && str(verbose, psis)
  verbose && exit(verbose)

  if (what=="MOR") {
    ## if  targeting mean  under  optimal rule,  then  estimating the  empirical
    ##  regret as well
    verbose && enter(verbose, "estimating the regret")
    
    regret <- estimateRegret(obs=obs, what=what, Qtab=Qtab, QEpstab=QEpstab,
                             epsHtab=epsHtab, psi=psis["psi"], verbose=verbose)
    
    this$.regret <- regret["regret"]
    
    this$.regret.sd <- regret["regret.sd"]

    verbose && str(verbose, regret)
    verbose && exit(verbose)
  }
  
  
  ## updating 'history'
  verbose && enter(verbose, "Updating history")
  step <- getStep(this)
  step["target"] <- step["target"]+1
  this$.step <- step
  updateHistory(this, "target")
  verbose && str(verbose, tsml.cara.rct::getHistory.TSMLCARA(this))
  verbose && exit(verbose)

  verbose && exit(verbose)
  verbose && ruler(verbose)

})



tsml.cara.rct <- structure(function#Targeted  Minimum Loss  Covariate-Adjusted  Response-Adaptive RCT Design and Statistical Analysis
### Simulates a targeted minimum loss covariate-adjusted response-adaptive RCT
### design and statistical analysis.
(what=c("ATE", "MOR"),
### A  \code{character}  indicating the  parameter  of  interest to  estimate.
### Either "ATE" for the Average  Treatment Effect, the difference between the
### means under  '\eqn{do(A=1)}' and  '\eqn{do(A=0)}', or  "MOR" for  the Mean
### under the Optimal treatment Rule '\eqn{do(A=r(W))}'.
  flavor=c("parametric", "lasso"),
### A \code{character}  indicating the  'flavor' of the  procedure. 
  ninit=50,
### An     \code{integer},    number    of     subjects    to     sample    at
### initialization. Defaults to 50.
  by=25,
### An  \code{integer},  number of  subjects  to  sample  at each  step  after
### initialization. Defaults to 25.
  nmax=500,
### An  \code{integer},  maximum  number  of  subjects to  sample  during  the
### trial. Must be larger than 'ninit+by'. Defaults to 500.
  tm.init=oneOne,
### A  \code{function}  describing  the  initial  treatment  mechanism  to  be
### employed.   Defaults  to  the  balanced  (1:1)  treatment  mechanism,  ie,
### \code{function} \code{\link{oneOne}}.
  tm.ref=oneOne,
### A  \code{function}  describing the  reference  treatment  mechanism to  be
### employed.   Defaults  to  the  balanced  (1:1)  treatment  mechanism,  ie,
### \code{function} \code{\link{oneOne}}.
  learnQ,
### A  model  \eqn{{\cal Q}}  of  conditional  expectations of  \eqn{Y}  given
### \eqn{(A,W)}  for  both  flavors  'parametric'  and  'lasso',  given  as  a
### \code{formula}  or  a  \code{function}  outputing  formulas.  Defaults  to
### \code{formula(Y~1)} for flavors  'parametric' and 'lasso'.
  tm.model=formula(A~1),
### A parametric model \eqn{{\cal G}}  of treatment mechanisms, used only when
### 'what'  equals  "ATE".   The   procedure  targets  the  optimal  treatment
### mechanism within this model.  Defaults to \code{formula(A~1)}.
  tm.control=glm.control(maxit=500),
### A  \code{list} of  options for  the targeting  of the  treatment mechanism
### within the  model defined by  argument 'tm.model'.  Used only  when 'what'
### equals "ATE", it defaults to \code{glm.control(maxit=500)}.
  Gmin=1e-2,
### A  small positive  \code{numeric}, with  default value  \code{1e-2}.  When
### \code{what} equals 'ATE', it is  the minimum value of elements of the
### parametric  model \eqn{{\cal  G}}  of treatment  mechanisms (see  argument
### \code{tm.model}).  The  maximum value is  \code{1-Gmin}.  When \code{what}
### equals 'MOR', it  is the minimum value of  the conditional probability
### of \eqn{A=r_n(W)} given \eqn{W}.
  Gexploit=Gmin,
### A small positive  \code{numeric}, with default value  that of \code{Gmin},
### or a  function of sample  size giving such  small numbers, only  used when
### \code{what} equals "MOR", in conjunction with \code{Gexplore}.
  Gexplore=1e-2,
### Either a small positive \code{numeric}, with default value \code{1e-2}, or
### a  function of  sample  size giving  such small  numbers,  only used  when
### \code{what} equals "MOR", in conjunction with \code{Gexploit}.
  Qmin=1e-2,
### A  small positive  \code{numeric}, the  minimum value  of  scaled outcomes
### \eqn{Y}. The maximum value is \code{1-Qmin}.
  conf.level=0.95,
###  A  \code{numeric},  the  confidence  level of  the  resulting  confidence
###  interval.
  verbose=FALSE,
### A \code{logical}  or an \code{integer}  indicating the level  of verbosity
### (defaults to 'FALSE').
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
  Bn=1e5,
### An \code{integer}, the sample size used  to estimate the true value of the
### data-adaptive parameter at  each step of the procedure  when 'what' equals
### 'MOR'. Defaults to 1e5.
  slice.by=1e5
### An \code{integer}. If it is smaller than argument 'n' of 'getSample', then
### the simulation  is decomposed  into 'n%/%slice.by' smaller  simulations of
### 'slice.by' observations and one of 'n%%slice.by' observations. Defaults to
### 1e5 (hence, no decomposition if 'n' smaller than 4e5). Mainly for internal
### use.
) {
  ##alias<< tsml.cara.rct
  ##references<< Chambaz, van der Laan,  Zheng, Chapter  16, Modern Adaptive Randomized Clinical  Trials: Statistical, Operational, and Regulatory    Aspects,    by    A.    Sverdlov    (CRC    Press,    2015). 
  ##seealso<< update, targetPsi, getSample
  
  ##details<<  Defines  a  lower-bound   on  the  conditional  probability  of
  ## \eqn{do(A=1-r_n(W))} given \eqn{W}.
  
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  ## Argument 'what'
  what <- match.arg(what)
  
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
  
  ## Argument 'ninit'
  ninit <- Arguments$getNumeric(ninit)
  
  ## Argument 'by'
  by <- Arguments$getNumeric(by)
  
  ## Argument 'nmax'
  nmax <- Arguments$getNumeric(nmax, c(ninit+by, Inf))
  
  ## Argument 'tm.init'
  mode <- mode(tm.init)
  if (mode != "function") {
    throw("Argument 'tm.init' should be of mode 'function', not '", mode)
  }
  
  ## Argument 'tm.ref'
  mode <- mode(tm.ref)
  if (mode != "function") {
    throw("Argument 'tm.ref' should be of mode 'function', not '", mode)
  }
  
  ## Argument 'tm.model'
  if (!identical(class(tm.model), "formula") | tm.model[[2]]!="A") {
    throw("Argument 'tm.model' should be a formula with response 'A', not:",
          deparse(substitute(form, list(form=tm.model)))) 
  }
  
  ## Argument 'tm.control'
  
  ## Argument 'Gmin'
  Gmin <- Arguments$getNumeric(Gmin, c(0, 1/2))
  
  ## Argument 'Gexploit'
  mode <- mode(Gexploit)
  if (!(mode %in% c("numeric", "function"))) {
    throw("Argument 'Gexploit' should be of mode either 'numeric' or 'function', not ", mode) 
  } else if (mode=="numeric") {
    Gexploit <- Arguments$getNumeric(Gexploit, c(0, 1/2))
  } 

  
  ## Argument 'Gexplore'
  mode <- mode(Gexplore)
  if (!(mode %in% c("numeric", "function"))) {
    throw("Argument 'Gexplore' should be of mode either 'numeric' or 'function', not ", mode) 
  } else if (mode=="numeric") {
    Gexplore <- Arguments$getNumeric(Gexplore, c(0, 1/2))
  } 
  
  ## Argument 'verbose'
  verbose <- Arguments$getVerbose(verbose)
  verbose <- less(verbose, 10)
  
  ## Argument 'piV':
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
  
  ## Argument 'Bn':
  Bn <- Arguments$getInteger(Bn);
  if (Bn <= 1e5) {
    verbose && str(verbose, "isn't 'Bn' too small?")
  }
  
  ## Argument 'slice.by'
  slice.by <- Arguments$getInteger(slice.by, c(1, Inf));
  
  
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Core
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  ## initial sampling (under 1:1 treatment mechanism)
  ## (a) sampling
  obs <- getSample(ninit, tm=tm.ref,
                   piV=piV, Qbar=Qbar, Vbar=Vbar,
                   family=family)
  nobs <- nrow(obs)
  ## (b) declaring the 'TSMLCARA' object
  tsml.cara <- TSMLCARA(what=what,
                        flavor=flavor,
                        obs=obs,
                        tm.ref=tm.ref,
                        tm.model=tm.model,
                        tm.control=tm.control,
                        learnQ=learnQ,
                        Gmin=Gmin,
                        Gexploit=Gexploit,
                        Gexplore=Gexplore,
                        Qmin=Qmin,
                        conf.level=conf.level,
                        verbose=verbose)
  
  ## update of 'tsml.cara'
  update(tsml.cara, verbose=verbose)
  targetPsi(tsml.cara, verbose=verbose)
  
  ## if 'what=="MOR"', then compute and store the targeted data-adaptive parameters
  if (what=="MOR") {
    Qn <- getQ(tsml.cara)
    ruleQn <- ruleQbar(Qn)
    ##
    ## first data-adaptive parameter:
    ## mean reward under the current best estimate of the optimal rule
    ##
    truthn <- getSample(Bn, tm=oneOne, rule=ruleQn,
                        piV=piV, Qbar=Qbar, Vbar=Vbar, what,
                        family=family, slice.by=slice.by)[c("psi", "psi.sd")]
    ##
    ## second data-adaptive parameter:
    ## empirical cumulated regret
    ##
    W <- extractW(obs)
    rA <- as.numeric(ruleQn(W))
    eregretn <- mean(obs[, "Y"] - Qbar(cbind(A=rA, W)))
    ##
    ## third data-adaptive parameter:
    ## counterfactual cumulated regret
    ##
    col <- match(c("Y0", "Y1"), colnames(obs))
    rA <- col[rA+1]
    counterfactual <- as.numeric(obs[cbind(1:nrow(obs), rA)])
    cregretn <- mean(obs[, "Y"] - counterfactual)
  }  
  ## successive updates of 'tsml.cara'
  while (nobs < nmax) {
    if (nobs+by>nmax) {
      by <- (nmax-nobs)
    }
    newObs <- getSample(by, tm=getGstar(tsml.cara),
                        piV=piV, Qbar=Qbar, Vbar=Vbar,
                        family=family)
    addNewSample(tsml.cara, newObs)
    nobs <- nrow(getObs(tsml.cara))
    update(tsml.cara, verbose=verbose)
    targetPsi(tsml.cara, verbose=verbose)
    if (what=="MOR") {
      Qn <- getQ(tsml.cara)
      ruleQn <- ruleQbar(Qn)
      ##
      ## first data-adaptive parameter:
      ## mean reward under the current best estimate of the optimal rule
      ##
      truthn <- rbind(truthn,
                      getSample(Bn, tm=oneOne, rule=ruleQn,
                                piV=piV, Qbar=Qbar, Vbar=Vbar, what,
                                family=family, slice.by=slice.by)[c("psi", "psi.sd")])
      ##
      ## second data-adaptive parameter:
      ## empirical cumulated regret
      ##
      obs <- getObs(tsml.cara)
      W <- extractW(obs)
      rA <- as.numeric(ruleQn(W))
      eregretn <- c(eregretn,
                    mean(obs[, "Y"] - Qbar(cbind(A=rA, W))))
      ##
      ## third data-adaptive parameter:
      ## counterfactual cumulated regret
      ##
      rA <- col[rA+1]
      counterfactual <- as.numeric(obs[cbind(1:nrow(obs), rA)])
      cregretn <- c(cregretn,
                    mean(obs[, "Y"] - counterfactual))
    }  
  }
  if (what=="MOR") {
    truthn[, "psi.sd"] <- truthn[, "psi.sd"]/sqrt(Bn)
    colnames(truthn) <- c("psin", "psin.sd")
    rownames(truthn) <- NULL
    attr(tsml.cara, "truthn") <- truthn
    attr(tsml.cara, "eregretn") <- eregretn
    attr(tsml.cara, "cregretn") <- cregretn
  }
  
  return(tsml.cara)
### Returns a  \code{TSMLCARA} object  which summarizes the  TSMLCARA undertaken
### procedure.
}, ex=function() {
  ##
  log <- Arguments$getVerbose(-1, timestamp=TRUE)
  set.seed(12345)
  
  ## ########################
  ## AVERAGE TREATMENT EFFECT
  ## ########################
  tm.model <- formula(A~.)
  psi.sd <- sqrt(getOptVar(n=1e5,
                           tm.model=tm.model,
                           piV=c(1/2, 1/3, 1/6),
                           family="gamma",
                           Qbar=Qbar1,
                           Vbar=Vbar1))
  truth <- c(psi=91/72, psi.sd=psi.sd)
  
  ## parametric example
  learnQ <- formula(Y~I(as.integer(A)):(U+V)+I(as.integer(1-A)):(U+V))
  ATE.param <- tsml.cara.rct(what="ATE",
                             flavor="parametric",
                             ninit=200,
                             by=100,
                             nmax=400,
                             tm.init=oneOne,
                             tm.ref=oneOne,
                             learnQ=learnQ,
                             tm.model=tm.model,
                             conf.level=0.95,
                             piV=c(1/2, 1/3, 1/6),
                             family="gamma",
                             Qbar=Qbar1,
                             Vbar=Vbar1)
  ATE.param
  ## Not run:
  plot(ATE.param, truth=truth)
  ## End(**Not run**)

  ## See the vignette for more examples...
})

############################################################################
## HISTORY:
## 2016-09-16
## o Created.
############################################################################

