estimateRegret  <-  function(#Estimates the Empirical Regret
### Computes the estimator of the empirical regret and its standard deviation.
                          obs,
### A \code{data.frame} of observations, as produced by \code{getSample}.
                          what="MOR",
### A \code{character} indicating the parameter of interest to estimate.  Must
### be 'MOR'.
                          Qtab=NULL,
### A \code{matrix} of  conditional means of \eqn{Y} given  \eqn{(A,W)} at the
### observations as they are currently estimated. 
                          QEpstab,
### A  \code{matrix} of  estimated conditional  expectations of  \eqn{Y} given
### \eqn{(A,W)} with columns \code{c("A", "A=0", "A=1")}.
                          epsHtab,
### A  \code{matrix} of  clever  covariates with  columns \code{c("A",  "A=0",
### "A=1")}.
                          psi,
### A \code{numeric}, the targeted estimator of the mean under the optimal rule.
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
  Y <- obs[, "Y"]

  ## Arguments 'what'
  if (what!="MOR") {
    throw("Argument 'what' must be  'MOR', not: ", what)
  }

  ## Argument 'Qtab', only used if 'what' equals 'MOR'
  if (what=="MOR") {
    Qtab <- Arguments$getNumerics(Qtab)
    test <- !is.matrix(Qtab) ||
        !identical(colnames(Qtab), c("A", "A=0", "A=1")) ||
            !(nrow(obs)==nrow(Qtab))
    if (test) {
      throw("Argument 'Qtab' does not meet the constraints.")
    }
  }
  
  ## Argument 'QEpstab'
  QEpstab <- Arguments$getNumerics(QEpstab)
  test <- !is.matrix(QEpstab) ||
      !identical(colnames(QEpstab), c("A", "A=0", "A=1")) ||
          !(nrow(obs)==nrow(QEpstab))
  if (test) {
    throw("Argument 'QEpstab' does not meet the constraints.")
  }

  ## Argument 'epsHtab'
  epsHtab <- Arguments$getNumerics(epsHtab)
  test <- !is.matrix(epsHtab) ||
      !identical(colnames(epsHtab), c("A", "A=0", "A=1")) ||
          !(nrow(obs)==nrow(epsHtab))
  if (test) {
    throw("Argument 'epsHtab' does not meet the constraints.")
  }

  ## Argument 'psi'
  psi <- Arguments$getNumeric(psi)
  
  ## Argument 'verbose'
  verbose <- Arguments$getVerbose(verbose)

  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Core
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ##
  ## regret
  ##

  regret <- mean(Y) - psi
  names(regret) <- NULL

  A <- obs[, "A"]
  rA <- as.numeric((Qtab[, "A=1"]-Qtab[, "A=0"])>0)
  rAequalsA <- (rA==A)
  col <- match(c("A=0", "A=1"), colnames(QEpstab))
  rA <- col[rA+1]

  if (FALSE) {
    ##
    ## original version
    ##
    psin <- mean(Qtab[cbind(1:nrow(Qtab), rA)])
    
    regret.sd <- sqrt(mean(
    (epsHtab[, "A"]*(Y-QEpstab[, "A"])*rAequalsA +
     QEpstab[cbind(1:nrow(QEpstab), rA)] - psi -
     Qtab[cbind(1:nrow(Qtab), rA)] + psin
    )^2))
  } else if (FALSE) {
    ##
    ## alternative version
    ##
    regret.sd <- sqrt(mean(
    (epsHtab[, "A"]*(Y-QEpstab[, "A"])*rAequalsA +
     QEpstab[cbind(1:nrow(QEpstab), rA)] - ## psi -
     Qtab[cbind(1:nrow(Qtab), rA)] ## + psin
    )^2))
  } else {
    ##
    ## yet another alternative version
    ##
    regret.sd <- sqrt(mean(
    (epsHtab[, "A"]*(Y-QEpstab[, "A"])*rAequalsA ## +
      ## QEpstab[cbind(1:nrow(QEpstab), rA)] - psi -
      ## Qtab[cbind(1:nrow(Qtab), rA)] ## + psin
    )^2))
  }
  out <- c(regret=regret, regret.sd=regret.sd)
    
  return(out)
### Returns a \code{vector}  giving the estimator of the  empirical and regret
### and its estimated standard deviation.
}



############################################################################
## HISTORY:
## 2016-09-16
## o Created.
############################################################################

