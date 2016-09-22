\name{addNewSample.TSMLCARA}
\alias{addNewSample.TSMLCARA}
\alias{addNewSample}
\title{Adds the Newly Sampled Observations}
\description{Adds the newly sampled \code{data.frame} of observations and updates }
\usage{\method{addNewSample}{TSMLCARA}(this, newobs, ...)}
\arguments{
  \item{this}{An object of class \code{TSMLCARA}.}
  \item{newobs}{The  newly  sampled  \code{data.frame}  of  observations  to  add  to  the
\code{TSMLCARA} object.}
  \item{\dots}{Not used.}
}
\details{The  output  \code{data.frame}  has a  "from"  attribute  which
indicates the  row in the  concatenated data set  of the first row  of the
added data set.}


\author{Antoine Chambaz [aut, cre]}



\seealso{tsml.cara.rct}
\examples{FALSE}
