\name{Vbar1}
\alias{Vbar1}
\title{A Conditional Variance of \eqn{Y} Given \eqn{(A,W)}}
\description{A  conditional   variance  of   \eqn{Y}  given   \eqn{(A,W)}  to   use  in
\code{\link{getSample}}.}
\usage{Vbar1(AW, rho = 1)}
\arguments{
  \item{AW}{A \code{data.frame} of observations,  whose columns contain the components
of \eqn{W} and the value of \eqn{A}.}
  \item{rho}{A non-negative \code{numeric}, with default value equal to one.}
}
\details{This  conditional variance  of \eqn{Y} given  \eqn{(A,W)} where
\eqn{W=(U,V)}     is    \deqn{(\rho+U+V)^2}     when    \eqn{A=1}     and
\deqn{(\frac{1}{\rho+V}+U)^2} when \eqn{A=0}.}
\value{Returns a \code{vector} of conditional variances given \eqn{(A,W)}.}
\references{Chambaz,  van der  Laan,  Scand.  J.  Stat.,  41(1):104--140
(2014).}
\author{Antoine Chambaz [aut, cre]}





