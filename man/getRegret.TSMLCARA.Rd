\name{getRegret.TSMLCARA}
\alias{getRegret.TSMLCARA}
\alias{getRegret}
\title{Returns the Current Estimator of the  Empirical Regret}
\description{Returns the current value of the estimator of the empirical regret when 'what' equals 'MOR'. }
\usage{\method{getRegret}{TSMLCARA}(this, ...)}
\arguments{
  \item{this}{An object of class \code{TSMLCARA}.}
  \item{\dots}{Not used.}
}

\value{The current estimate of the empirical regret. }

\author{Antoine Chambaz [aut, cre]}



\seealso{tsml.lcara.rct, getHistory, getRegretSd}
\examples{FALSE}
