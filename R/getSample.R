getSample <- structure(
    function#Generates  Data
### Generates  a run  of simulated  observations of  the  form \eqn{(W,G,A,Y,Y1,Y0)}
### under  the  specified treatment  mechanism.  
    (n,
### An \code{integer}, the number of observations to be generated.
     tm=oneOne,
### A  \code{function}  describing the  treatment  mechanism  to be  employed.
### Defaults to  the balanced  (1:1) treatment mechanism,  ie, \code{function}
### \code{\link{oneOne}}. The \code{G} column equals the vector \code{tm(W)}.
     rule=NULL,
### Either 'NULL' (default value) or  a \code{function} describing the rule to
### be employed  when \code{what} equals  "MOR".  In that  case, it must  be a
### function of  \eqn{W} with values  in \eqn{\{0,1\}}. If  \code{what} equals
### "MOR" and \code{rule} is 'NULL', then the optimal treatment rule is used.
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
     what,
### A \code{character}. If it  is not missing then it must  be equal to either
### "ATE"  (Average  Treatment  Effect)  or  "MOR"  (Mean  under  the  Optimal
### treatment Rule).  In  that case, ONLY estimates of the  true parameter and
### optimal  standard deviation  under the  specified simulation  scheme (and,
### possibly, treatment mechanism) are computed and returned.
      slice.by=n
### An \code{integer}. If it is smaller than argument 'n', then the simulation
### is  decomposed  into  'n%/%slice.by'  smaller  simulations  of  'slice.by'
### observations and one  of 'n%%slice.by' observations. Defaults  to 'n' (ie,
### no decomposition). Mainly for internal use.
###
     ) {
      ##references<< Chambaz,  van der Laan, Scand.  J.  Stat., 41(1):104--140
      ## (2014).

      ##details<<  By default,  implements the  same simulation  scheme as  in
      ## Chambaz, van der Laan, Scand.  J. Stat., 41(1):104--140 (2014).
      
      ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ## Validate arguments
      ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      ## Argument 'n'
      n <- Arguments$getInteger(n, c(1, Inf))
      
      ## Argument 'tm'
      mode <- mode(tm)
      if (mode != "function") {
        throw("Argument 'tm' should be of mode 'function', not '", mode)
      }

      ## Argument 'piV'
      piV <- Arguments$getNumerics(piV, range=c(0,1))
      if (sum(piV)!=1) {
        throw("Argument 'piV' should consist of non-negative weights summing to one.") 
      }

      ## Argument 'family'
      family <- match.arg(family)
      
      ## Argument 'Qbar'
      mode <- mode(Qbar);
      if (mode != "function") {
        throw("Argument 'Qbar' should be of mode 'function', not '", mode)
      }

      ## Argument 'Vbar':
      mode <- mode(Vbar);
      if (mode != "function") {
        throw("Argument 'Vbar' should be of mode 'function', not '", mode)
      }

      ## Argument 'what':
      if (missing(what)) {
        truth <- FALSE
      } else {
        truth <- TRUE
        if (!(what %in% c("ATE", "MOR"))) {
          throw("Argument 'what' must be equal to either 'ATE' or 'MOR', not ", what)
        }
      }

      ## Argument 'rule'
      if (!is.null(NULL)) {
        if (what=="ATE") {
          warning("Argument 'rule' meaningless when 'what' equals 'ATE'!")
        } else if (mode(rule)!="function") {
          throw("Argument 'rule' should be of mode 'function', not '", mode(rule))
        }
      }

      
      ## Argument 'slice.by'
      slice.by <- Arguments$getInteger(slice.by, c(1, Inf));
      slice <- ifelse (slice.by<n, TRUE, FALSE)
      if (slice) { ## n=pp*qq+rr
        qq <- slice.by
        pp <- n%/%slice.by
        rr <- n%%slice.by
      }
      
      ##  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ## Core
      ##  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      ## sampling W
      V <- 1 + findInterval(runif(n), cumsum(piV))
      V <- factor(V, levels=1:length(piV))
      U <- runif(n)
      W <- data.frame(V=V, U=U)

      ## computing G
      G <- tm(W)
      
      ## sampling A|W
      A <- rbinom(n, size=1, prob=G)
      AW <- cbind(A=A, W)

      ##
      ## sampling Y|A,W
      ##

      ## intervention 'A=1'
      AW1 <- cbind(A=1, W)
      if (slice) {
        avg1 <- rep(NA, n)
        sig1 <- rep(NA, n)
        for (ii in 1:pp) {
          idx <- ((ii-1)*qq+1):(ii*qq)
          avg1[idx] <- Qbar(AW1[idx, , drop=FALSE])
          sig1[idx] <- Vbar(AW1[idx, , drop=FALSE])
        }
        if (rr>0) {
          idx <- (pp*qq+1):(pp*qq+rr)
          avg1[idx] <- Qbar(AW1[idx, , drop=FALSE])
          sig1[idx] <- Vbar(AW1[idx, , drop=FALSE])
        }
      } else {
        avg1 <- Qbar(AW1)
        sig1 <- Vbar(AW1)
      }
      if (family=="gamma") {
        ##
        ## Gamma distribution
        ##
        scale1 <- sig1/avg1
        Y1 <- rgamma(n, shape=avg1/scale1, scale=scale1)
                                        # yes indeed, 'rgamma' can take
                                        # vectors of shapes and scales
      } else if (family=="beta") {
        ##
        ## Beta distribution
        ##
        alpha <- ((1-avg1)/sig1 - 1/avg1)*avg1^2
        beta <- alpha*(1/avg1-1)
        Y1 <- rbeta(n, shape1=alpha, shape2=beta)
      }
      ## intervention 'A=0'
      AW0 <- cbind(A=0, W)
      if (slice) {
        avg0 <- rep(NA, n)
        sig0 <- rep(NA, n)
        for (ii in 1:pp) {
          idx <- ((ii-1)*qq+1):(ii*qq)
          avg0[idx] <- Qbar(AW0[idx, , drop=FALSE])
          sig0[idx] <- Vbar(AW0[idx, , drop=FALSE])
        }
        if (rr>0) {
          idx <- (pp*qq+1):(pp*qq+rr)
          avg0[idx] <- Qbar(AW0[idx, , drop=FALSE])
          sig0[idx] <- Vbar(AW0[idx, , drop=FALSE])
        }
        rm(idx)
      } else {
        avg0 <- Qbar(AW0)
        sig0 <- Vbar(AW0)
      }
      if (family=="gamma") {
        ##
        ## Gamma distribution
        ##
        scale0 <- sig0/avg0
        Y0 <- rgamma(n, shape=avg0/scale0, scale=scale0)
                                        # yes indeed, 'rgamma' can take
                                        # vectors of shapes and scales
      } else if (family=="beta") {
        ##
        ## Beta distribution
        ##
        alpha <- ((1-avg0)/sig0 - 1/avg0)*avg0^2
        beta <- alpha*(1/avg0-1)
        Y0 <- rbeta(n, shape1=alpha, shape2=beta)
      } 

      ## Y and avg
      Y <- A*Y1+(1-A)*Y0
      avg <- A*avg1 + (1-A)*avg0
      
      if (truth) {
        ## option 'truth'
        if (what=="ATE") {
          psi <- mean(Y1-Y0)
          psi.sd <- sqrt(mean((avg1-avg0-psi + (2*A-1)/blik(A,G)*(Y-avg))^2))
          truth <- c(psi=psi, psi.sd=psi.sd)
        } else if (what=="MOR") {
          if (is.null(rule)) {
            ## optimal treatment assignment
            rA <- (avg1>avg0)
          } else {
            if (slice) {
              rA <- rep(NA, n)
              for (ii in 1:pp) {
                idx <- ((ii-1)*qq+1):(ii*qq)
                rA[idx] <- rule(W[idx, , drop=FALSE])
              }
              if (rr>0) {
                idx <- (pp*qq+1):(pp*qq+rr)
                rA[idx] <- rule(W[idx, , drop=FALSE])
              }
              rm(idx)
            } else {
              rA <- rule(W)
            }
          }
          rY <- rA*Y1+(1-rA)*Y0
          rAvg <- rA*avg1 + (1-rA)*avg0
          ##
          psi <- mean(rY)
          psi.sd.opt <- sqrt( mean( (rAvg-psi + (rY-rAvg))^2 ) )
          psi.sd <- sqrt( mean( (rAvg-psi + (A==rA)/blik(A,G)*(Y-avg))^2 ) )
          truth <- c(psi=psi, psi.sd=psi.sd, psi.sd.opt=psi.sd.opt)
        }
        out <- truth
      } else {
        out <- cbind(W, G=G, A=A, Y=Y, Y1=Y1, Y0=Y0)
      }
      
      return(out)
### If \code{what} is 'NULL' then  returns a \code{data.frame} of observations
### with columns  \code{Y1} and \code{Y0} (counterfactual  outcomes), \code{Y}
### (actual outcome), \code{A} (assigned treatment), \code{G} (the conditional
### probability that A=1  given W), the rest being  interpreted as covariates.
### Otherwise,  returns  the  estimated  true  parameter  value  and  standard
### deviation of the efficient influence  curve under the specified simulation
### scheme and  treatment mechanism, if  'what' equals "ATE", or  the standard
### deviationS of the efficient influence curve under the specified simulation
### scheme  and (i)  treatment mechanism  and (ii)  either the  treatment rule
### given by  argument \code{rule}  or the optimal  treatment rule,  if 'what'
### equals "MOR".
    }, ex=function() {
      ## Setting the verbosity parameter
      library(R.utils)
      log <- Arguments$getVerbose(-8, timestamp=TRUE)

      ##
      obs <- getSample(10)
    })



############################################################################
## HISTORY:
## 2016-09-16
## o Created.
############################################################################

