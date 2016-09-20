validateArgumentObs <- function(obs) {
  if (!is.data.frame(obs)) {
    throw("Argument 'obs' should be a matrix")
  }
  
  ## Ensuring that columns 'G', 'A', 'Y', 'Y1' and 'Y0' are present
  varNames <- c("G", "A", "Y", "Y1", "Y0")
  nms <- colnames(obs)
  m <- match(varNames, nms)
  if (nrow(obs)>0) {
    idxs <- which(is.na(m))
    if (length(idxs)) {
      throw("Missing column:", varNames[idxs])
    }
  }
  
  obs
}

############################################################################
## HISTORY:
## 2016-09-16
## o Created.
############################################################################

