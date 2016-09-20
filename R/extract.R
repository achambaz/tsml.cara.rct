extract <- function(dat, exclude) {
  theVar <- setdiff(colnames(dat), exclude)
  dat[, theVar, drop=FALSE]
}

extractW <- function(#Extracts W Columns from Data.Frame of Observations
### Extracts the \code{W} column(s) from a \code{data.frame} of observations.
    dat
### A  \code{data.frame}  of  observations,  as  the  \code{obs}  argument  of
### \code{function} \code{\link{tsmle.npvi}}.
    ) {
  ##details<< Mainly for internal use.
  extract(dat, exclude=c("G", "A", "Y", "Y1", "Y0"))
### The  \code{data.frame}  extracted from  \code{dat}  by  removing the  five
### \code{G}, \code{A}, \code{Y}, \code{Y1} and \code{Y0} columns.
}

extractG <- function(#Extracts G Column from Data.Frame of Observations
### Extracts the \code{G} column from a \code{data.frame} of observations.
                     dat
### A  \code{data.frame}  of  observations,  as  the  \code{obs}  argument  of
### \code{function} \code{\link{tsmle.npvi}}.
                     ) {
  ##details<< Mainly for internal use.
  dat[, "G", drop=FALSE]
### The \code{data.frame} extracted from  \code{dat} by removing the \code{W},
### \code{A} and \code{Y} columns.
}


extractAW <- function(#Removes the G,  Y, Y1 and Y0 Columns from Data.Frame of Observations
### Removes  the \code{G}, \code{Y},  \code{Y1} and  \code{Y0} columns  from a
### \code{data.frame} of observations.
    dat
### A \code{data.frame}  of observations,  as produced by  the \code{function}
### \code{\link{getSample}}.
    ) {
  ##details<< Mainly for internal use.
  dat <- extract(dat, exclude=c("G", "Y", "Y1", "Y0"))

  ## enforcing the order to be A, W
  cn <- colnames(dat)
  wA <- which(cn=="A")
  cbind(dat[, wA, drop=FALSE], dat[, -wA, drop=FALSE])
### The \code{data.frame} extracted from  \code{dat} by removing the \code{G},
###  \code{Y}, \code{Y1}  and \code{Y0} columns in  such a way  that the first
###  column is \code{A}.
}

############################################################################
## HISTORY:
## 2016-09-16
## o Created.
############################################################################

