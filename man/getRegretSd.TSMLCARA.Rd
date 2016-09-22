\name{getRegretSd.TSMLCARA}
\alias{getRegretSd.TSMLCARA}
\alias{getRegretSd}
\title{Returns the Current Estimated Standard Deviation of the Current Estimator of the Empirical Regret}
\description{Returns  the current  value of  the  estimated standard  deviation of  the
current estimator of the empirical regret when 'what' equals 'MOR'.}
\usage{\method{getRegretSd}{TSMLCARA}(this, ...)}
\arguments{
  \item{this}{An object of class \code{TSMLCARA}.}
  \item{\dots}{Not used.}
}

\value{The  estimated  standard  deviation  of  the current  estimator  of  the
empirical regret.}

\author{Antoine Chambaz [aut, cre]}



\seealso{tsml.lcara.rct, getHistory, getRegret}
\examples{FALSE}
