\name{Vbar2}
\alias{Vbar2}
\title{A Conditional Variance of \eqn{Y} Given \eqn{(A,W)}}
\description{A  conditional   variance  of   \eqn{Y}  given   \eqn{(A,W)}  to   use  in
\code{\link{getSample}}. }
\usage{Vbar2(AW, rho = 0.1)}
\arguments{
  \item{AW}{A \code{data.frame} of observations,  whose columns contain the components
of \eqn{W} and the value of \eqn{A}.}
  \item{rho}{A non-negative \code{numeric}, with default value equal to one.}
}
\details{This  conditional expectation  of \eqn{Y} given  \eqn{(A,W)} is
always equal to \eqn{\rho}.}
\value{Returns a \code{vector} of conditional variances given \eqn{(A,W)}.}

\author{Antoine Chambaz [aut, cre]}





