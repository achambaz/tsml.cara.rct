targetOptRule  <-  function(#Targets the Optimal Treatment Rule.
### Function to target the optimal treatment rule.
                            this,
### A \code{TSMLCARA} object, as created by \code{TSMLCARA}.
                            Q,
### A \code{function},  the fitted  estimating the conditional  expectation of
### \eqn{Y} given \eqn{(A,W)}.
                            ...,
### Additional parameters.
                          verbose=FALSE
### A \code{logical}  or an \code{integer}  indicating the level  of verbosity
### (defaults to 'FALSE').
                          ) {
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ## Argument 'Gmin'
  Gmin <- getGmin(this)

  ## Argument 'Gexploit'
  GEXPLOIT <- getGexploit(this)
  mode <- mode(GEXPLOIT)
  if (mode=="numeric") {
    Gexploit <- GEXPLOIT
  } else {
    Gexploit <- Arguments$getNumerics(GEXPLOIT(nrow(getObs(this))), c(0, 1/2))
  }
  Gexploit <- max(Gexploit, Gmin)

  ## Argument 'Gexplore'
  GEXPLORE <- getGexplore(this)
  mode <- mode(GEXPLORE)
  if (mode=="numeric") {
    Gexplore <- GEXPLORE
  } else {
    Gexplore <- Arguments$getNumerics(GEXPLORE(nrow(getObs(this))), c(0, 1/2))
  }
  
  ## Argument 'Q'
  if (!mode(Q)=="function"|is.null(attr(Q, "model"))) {
    throw("Argument 'Q' must be a function as output by 'estimateQ'.")
  }
  
  ## Argument 'verbose'
  verbose <- Arguments$getVerbose(verbose)

  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Core
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  coeffs <- attr(Q, "model")
  ## --
  ## old version, kept for documenting
  ## --
  ## Gstar <- function(W) {
  ##   q <- Q(1,W) - Q(0,W)
  ##   pos <- (q > Gexplore)
  ##   neg <- (-q > Gexplore)
  ##   btwn <- !(pos|neg)
  ##   ##
  ##   out <- q
  ##   out[pos] <- (1-Gexploit)
  ##   out[neg] <- Gexploit
  ##   out[btwn] <- pmin(1-Gmin, pmax(Gmin, (q[btwn]+Gexplore)/(2*Gexplore)))
  ##   return(out)
  ## }

  ## --
  ## new version, based on 'smoothIndicator'
  ## --
  Gstar <- function(W) {
    q <- Q(1,W) - Q(0,W)
    out <- smoothIndicator(q, Gexploit, Gexplore)
    return(out)
  }
  attr(Gstar, "model") <- coeffs
  attr(Gstar, "expl") <- c(Gexploit=Gexploit, Gexplore=Gexplore)
  
  return(Gstar)
### Returns \code {Gstar}, a \code{function} of the same form as \code{oneOne}
### coding  the  treatment  rule.   This  function  has  a  "model"  attribute
### describing the fitted conditional expectation of \eqn{Y} given \eqn{(A,W)}
### which  is used  for its  construction and  a "expl"  attribute, describing
### which smooth approximation of \eqn{1\{x>0\}} over \eqn{[-1,+1]} is used.
}

############################################################################
## HISTORY:
## 2016-09-16
## o Created.
############################################################################

