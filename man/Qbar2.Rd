\name{Qbar2}
\alias{Qbar2}
\title{A conditional Expectation of \eqn{Y} Given \eqn{(A,W)}}
\description{A  conditional  expectation  of  \eqn{Y}   given  \eqn{(A,W)}  to  use  in
\code{\link{getSample}}.}
\usage{Qbar2(AW, rho = 1)}
\arguments{
  \item{AW}{A \code{data.frame} of observations,  whose columns contain the components
of \eqn{W} and the value of \eqn{A}.}
  \item{rho}{A non-negative \code{numeric}, with default value equal to 1.5.}
}
\details{This conditional expectation of \eqn{Y} given \eqn{(A,W)} where
\eqn{W=(U,V)} is \deqn{\frac{1}{2} + \frac{3}{8}*\cos(\rho*V*\pi*U)} when
\eqn{A=1}  and \deqn{\frac{1}{2}  + \frac{1}{4}*\sin(3*\rho/V*\pi*U)} when
\eqn{A=0}.}
\value{Returns a \code{vector} of conditional expectations given \eqn{(A,W)}.}
\references{Chambaz,        Zheng,       van        der       Laan,
https://hal.archives-ouvertes.fr/hal-01301297.}
\author{Antoine Chambaz [aut, cre]}





