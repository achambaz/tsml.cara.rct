\name{as.character.TSMLCARA}
\alias{as.character.TSMLCARA}
\alias{as.character}
\title{Returns a Description of a TSMLCARA Object}
\description{Returns a short string describing the TSMLCARA object.}
\usage{\method{as.character}{TSMLCARA}(x, ...)}
\arguments{
  \item{x}{An object of class \code{TSMLCARA}.}
  \item{\dots}{Not used.}
}

\value{A character string summarizing the content of the object. The summary contains:
\itemize{
\item{The name of the parameter of interest.}
\item{The sample size of the data set involved in the procedure.}
\item{The value of the TSMLE and its estimated standard error.}
\item{A confidence interval with default level of 95\% (the level can be changed by using \code{\link{setConfLevel}}).}
}}

\author{Antoine Chambaz [aut, cre]}



\seealso{tsml.cara.rct, setConfLevel}
\examples{FALSE}
