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
### A \code{TSMLECARA} object, as created by \code{TSMLCARA}.
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

  ## if  targeting mean  under  optimal rule,  then  estimating the  empirical
  ##  regret as well
  verbose && enter(verbose, "estimating the regret")

  if (what=="MOR") {
    regret <- estimateRegret(obs=obs, what=what, Qtab=Qtab, QEpstab=QEpstab,
                             epsHtab=epsHtab, psi=psis["psi"], verbose=verbose)
    
    this$.regret <- regret["regret"]
    
    this$.regret.sd <- regret["regret.sd"]
  }
  
  verbose && str(verbose, regret)
  verbose && exit(verbose)
  
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
### Either 'ATE' for  the difference between the means  under '\eqn{A<-1}' and
### '\eqn{A<-0}'  or 'MOR'  for  the  mean under  the  optimal treatment  rule
### '\eqn{A<-r(W)}'.
  flavor=c("parametric", "lasso"),
### A \code{character}  indicating the  'flavor' of the  procedure. 
  ninit=50,
### An     \code{integer},    number    of     subjects    to     sample    at
### initialization. Defaults to 50.
  by=25,
### An  \code{integer},  number of  subjects  to  sample  at each  step  after
### initialization. Defaults to 25.
  nmax=500,
### An  \code{integer},  maximal  number  of  subjects to  sample  during  the
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
### 'what'  equals 'ATE'.  The  procedure targets  the optimal  treatment
### mechanism within this model.  Defaults to \code{formula(A~1)}.
  tm.control=glm.control(maxit=500),
### A  \code{list} of  options for  the targeting  of the  treatment mechanism
### within the  model defined  by argument 'tm.model'.  Used only  when 'what'
### equals 'ATE', it defaults to \code{glm.control(maxit=500)}.
  Gmin=1e-2,
### A  small positive  \code{numeric}, with  default value  \code{1e-2}.  When
### \code{what} equals 'ATE', it is  the minimum value of elements of the
### parametric  model \eqn{{\cal  G}}  of treatment  mechanisms (see  argument
### \code{tm.model}).  The  maximal value is  \code{1-Gmin}.  When \code{what}
### equals 'MOR', it  is the minimum value of  the conditional probability
### of \eqn{A=r_n(W)} given \eqn{W}.
  Gexpl=1e-2,
### A small positive \code{numeric}, with default value \code{1e-2}, only used
### when  \code{what} equals  'MOR',  in which  case  it lower-bounds  the
### conditional  probability  of  \eqn{A=1-r_n(W)}  given  \eqn{W}.
  threxpl=1e-2,
### Either a small positive \code{numeric}, with default value \code{1e-2}, or
### a  \code{function} of  sample size  giving  such numbers,  only used  when
### \code{what}  equals 'MOR'.  If  \eqn{0\in[Q_n-\theta,Q_n+\theta]} with
### \eqn{\theta} equal  to \code{threxpl} then \eqn{r_n(W)}  is the proportion
### of the interval which lies  above 0, thresholded at levels \code{Gmin} and
### \code{1-Gmin}.
  Qmin=1e-2,
### A  small positive  \code{numeric}, the  minimum value  of  scaled outcomes
### \eqn{Y}. The maximal value is \code{1-Qmin}.
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
  
  ## Argument 'Gexpl'
  Gexpl <- Arguments$getNumeric(Gexpl, c(0, 1/2))
  
  ## Argument 'threxpl'
  mode <- mode(threxpl)
  if (!(mode %in% c("numeric", "function"))) {
    throw("Argument 'threxpl' should be of mode either 'numeric' or 'function', not ", mode) 
  } else if (mode=="numeric") {
    threxpl <- Arguments$getNumeric(threxpl, c(0, 1/2))
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
                        Gexpl=Gexpl,
                        threxpl=threxpl,
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
  
  ## ## ########################
  ## ## AVERAGE TREATMENT EFFECT
  ## ## ########################
  
  ## psi.sd <- sqrt(getOptVar(n=1e5))
  ## truth <- c(psi=91/72, psi.sd=psi.sd)
  
  ## ## parametric example
  ## tm.model <- formula(A~.)
  
  ## learnQ <- formula(Y~I(as.integer(A)):(U+V)+I(as.integer(1-A)):(U+V))
  ## rd.param <- tsml.cara.rct(what="ATE",
  ##                          flavor="parametric",
  ##                          ninit=250,
  ##                          by=100,
  ##                          nmax=2500,
  ##                          tm.init=oneOne,
  ##                          tm.ref=oneOne,
  ##                          learnQ=learnQ,
  ##                          tm.model=tm.model,
  ##                          conf.level=0.95,
  ##                          piV=c(1/2, 1/3, 1/6),
  ##                          Qbar=Qbar1,
  ##                          Vbar=Vbar1)
  
  ## rd.param
  ## plot(rd.param, truth=truth)
  
  ## ## lasso example
  ## learnQ <- formula(Y~I(as.integer(A)):(poly(U, 3)*V)+I(as.integer(1-A)):(poly(U, 3)*V))
  ## rd.lasso <- tsml.cara.rct(what="ATE",
  ##                          flavor="lasso",
  ##                          ninit=250,
  ##                          by=100,
  ##                          nmax=2500,
  ##                          learnQ=learnQ,
  ##                          tm.init=oneOne,
  ##                          tm.ref=oneOne,
  ##                          tm.model=tm.model,
  ##                          conf.level=0.95,
  ##                          piV=c(1/2, 1/3, 1/6),
  ##                          Qbar=Qbar1,
  ##                          Vbar=Vbar1)
  
  ## rd.lasso
  ## x11()
  ## plot(rd.lasso, truth=truth)
  
  ## ## manually, another parametric example
  ## learnQ <- formula(Y~A*V)
  ## tm.model <- formula(A~V)
  
  ## obs <- getSample(200, 
  ##                  piV=c(1/2, 1/3, 1/6),
  ##                  Qbar=Qbar1,
  ##                  Vbar=Vbar1)
  ## nobs <- nrow(obs)
  ## rd.param2 <- TSMLCARA(what="ATE",
  ##                      obs=obs,
  ##                      learnQ=learnQ,
  ##                      tm.model=tm.model,
  ##                      Gmin=1e-2,
  ##                      flavor="parametric",
  ##                      verbose=log)
  ## while (nobs<=900) {
  ##   update(rd.param2, verbose=log)
  ##   targetPsi(rd.param2, verbose=log)
  ##   newObs <- getSample(100, tm=getGstar(rd.param2),
  ##                       piV=c(1/2, 1/3, 1/6),
  ##                       Qbar=Qbar1,
  ##                       Vbar=Vbar1)
  ##   addNewSample(rd.param2, newObs)
  ##   nobs <- nrow(getObs(rd.param2))
  ## }
  ## update(rd.param2, verbose=log)
  ## targetPsi(rd.param2, verbose=log)
  
  ## rd.param2
  ## x11()
  ## plot(rd.param2, truth=truth)
  
  ## ##################################
  ## MEAN OF THE OPTIMAL TREATMENT RULE
  ## ##################################
  
  set.seed(54321)
  family <- "beta"
  truth <- getSample(1e5,
                     tm=oneOne,
                     rule=NULL,
                     piV=c(1/2, 1/3, 1/6),
                     Qbar=Qbar2,
                     Vbar=Vbar2,
                     what="MOR",
                     family=family)
  ## truth <- c(psi=2.58445, psi.sd=0.48729) ## CHECK!!!
  
  ## parametric example
  learnQ <- formula(Y~I(A==1):(U+V)+I(A==0):(U+V))
  or.param <- tsml.cara.rct(what="MOR",
                            flavor="parametric",
                            ninit=200,
                            by=100,
                            nmax=500,
                            tm.init=oneOne,
                            tm.ref=oneOne,
                            learnQ=learnQ,
                            conf.level=0.95,
                            piV=c(1/2, 1/3, 1/6),
                            Qbar=Qbar2,
                            Vbar=Vbar2,
                            family=family)
  or.param
  ## plot(or.param, truth=truth)
  
  ## lasso example
  learnQ <- formula(Y~I(A==0)*poly(U, 3)*V+I(A==1)*poly(U, 3)*V)
  or.lasso <- tsml.cara.rct(what="MOR",
                            flavor="lasso",
                            ninit=200,
                            by=100,
                            nmax=500,
                            learnQ=learnQ,
                            tm.init=oneOne,
                            tm.ref=oneOne,
                            conf.level=0.95,
                            piV=c(1/2, 1/3, 1/6),
                            Qbar=Qbar2,
                            Vbar=Vbar2,
                            family=family)
  or.lasso
  ## x11()
  ## plot(or.lasso, truth=truth)
  
  ## manually, another parametric example
  learnQ <- formula(Y~I(A==1):U+I(A==0):U)
  
  obs <- getSample(200,
                   piV=c(1/2, 1/3, 1/6),
                   Qbar=Qbar2,
                   Vbar=Vbar2,
                   family=family)
  nobs <- nrow(obs)
  or.param2 <- tsml.cara.rct:::TSMLCARA(what="MOR",
                                        obs=obs,
                                        learnQ=learnQ,
                                        Gmin=1e-2,
                                        Qmin=1e-2,
                                        flavor="parametric", #"lasso",
                                        verbose=log)
  while (nobs<=400) {
    tsml.cara.rct:::update.TSMLCARA(or.param2, verbose=log)
    tsml.cara.rct:::targetPsi.TSMLCARA(or.param2, verbose=log)
    newObs <- getSample(100, tm=tsml.cara.rct:::getGstar.TSMLCARA(or.param2),
                        piV=c(1/2, 1/3, 1/6),
                        Qbar=Qbar2,
                        Vbar=Vbar2,
                        family=family)
    tsml.cara.rct:::addNewSample.TSMLCARA(or.param2, newObs)
    nobs <- nrow(getObs(or.param2))
  }
  tsml.cara.rct:::update.TSMLCARA(or.param2, verbose=log)
  tsml.cara.rct:::targetPsi.TSMLCARA(or.param2, verbose=log)
  
  or.param2
  ## x11()
  ## plot(or.param2, truth=truth)
  
  
})

############################################################################
## HISTORY:
## 2014-06-13
## o Adapted to deal with option 'MOR'
## 2014-03-25
## o Removed 'wExpand', 'covLegendre', 'degLeg' and related stuff
## 2014-02-26
## o Created.
############################################################################

