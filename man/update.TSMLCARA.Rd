\name{update.TSMLCARA}
\alias{update.TSMLCARA}
\alias{update}
\title{Updates a TSMLCARA Object}
\description{Updates a TSMLCARA object.}
\usage{\method{update}{TSMLCARA}(this, flavor = c("parametric", "lasso"), ..., 
    verbose = FALSE)}
\arguments{
  \item{this}{A \code{TSMLCARA} object, as created by \code{TSMLCARA}.}
  \item{flavor}{A \code{character}  indicating the  'flavor' of the  procedure. }
  \item{\dots}{Additional parameters.}
  \item{verbose}{A \code{logical}  or an \code{integer}  indicating the level  of verbosity
(defaults to 'FALSE').}
}


\references{Chambaz, van  der Laan,  Zheng, Chapter  16, Modern Adaptive Randomized Clinical  Trials: Statistical, Operational, and Regulatory    Aspects,    by    A.    Sverdlov    (CRC    Press,    2015).}
\author{Antoine Chambaz [aut, cre]}



\seealso{tsml.cara.rct, targetPsi}
\examples{FALSE}
