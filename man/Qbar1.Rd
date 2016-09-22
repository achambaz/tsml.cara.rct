\name{Qbar1}
\alias{Qbar1}
\title{A Conditional Expectation of \eqn{Y} Given \eqn{(A,W)}}
\description{A  conditional  expectation  of  \eqn{Y}   given  \eqn{(A,W)}  to  use  in
\code{\link{getSample}}.}
\usage{Qbar1(AW, rho = 1)}
\arguments{
  \item{AW}{A \code{data.frame} of observations,  whose columns contain the components
of \eqn{W} and the value of \eqn{A}.}
  \item{rho}{A non-negative \code{numeric}, with default value equal to one.}
}
\details{This conditional expectation of \eqn{Y} given \eqn{(A,W)} where
\eqn{W=(U,V)} is  \deqn{\rho * V  + 2*U^2 + 2*U  + 1} when  \eqn{A=1} and
\deqn{\rho/(1+V) + 2*U^2 + 2*U  + 1} when \eqn{A=0}.}
\value{Returns a \code{vector} of conditional expectations given \eqn{(A,W)}.}
\references{Chambaz, van  der  Laan, Scand.   J.  Stat.,  41(1):104--140
(2014).}
\author{Antoine Chambaz [aut, cre]}





