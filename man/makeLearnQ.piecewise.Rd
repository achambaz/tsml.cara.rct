\name{makeLearnQ.piecewise}
\alias{makeLearnQ.piecewise}
\title{Builds a Parametric Model Based on Sample Size}
\description{Builds a Parametric Working Model Based on Sample Size}
\usage{makeLearnQ.piecewise()}
\details{This functions builds a sample-size-dependent parametric working
model.   Two fine-tune  parameters  are determined  based  on sample  size
\eqn{n}: \eqn{deg=\min(\lceil  n/100\rceil, 6)}  and \eqn{nlev=\min(\lceil
n/100\rceil, 10)}.  If  \eqn{nlev} equals one, then the model  is given by
the   \code{formula}    \deqn{Y~I(A=0)*V*poly(U)+I(A=1)*V*poly(U)}   where
\eqn{poly(U)} consists  of \eqn{deg} orthogonal  polynoms of degrees  1 to
\eqn{deg}. If \eqn{nlev} is larger than or equal to two, then the model is
given  by   the  \code{formula}   \deqn{Y~I(A=0)*V*(poly(U)+\sum_k  step_k
(U))+I(A=1)*V*(poly(U)+\sum_k step_k (U))} where \eqn{poly(U)} consists of
\eqn{deg} orthogonal  polynoms of degrees  1 to \eqn{deg}  and \eqn{step_k
(U)} consits of \eqn{nlev} ortogonal polynoms of degrees 1 to \eqn{deg} in
\eqn{U} times the  indicator function of the subset of  \eqn{[0,1]} of the
form \eqn{\{x:\frac{k}{nlev}\le x <\frac{k}{nlev}\}}.}
\value{Returns a \code{formula} describing the parametric model.}

\author{Antoine Chambaz [aut, cre]}





