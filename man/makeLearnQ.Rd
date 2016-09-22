\name{makeLearnQ}
\alias{makeLearnQ}
\title{Builds a Parametric Working Model Based on Sample Size}
\description{Builds a Parametric Working Model Based on Sample Size}
\usage{makeLearnQ()}
\details{This functions builds a sample-size-dependent parametric working
model.   Two fine-tune  parameters  are determined  based  on sample  size
\eqn{n}:   \eqn{deg=3  +   \lfloor   n/500\rfloor}  and   \eqn{nlev=\lceil
n/250\rceil}.  If  \eqn{nlev} equals one, then  the model is given  by the
\code{formula}       \deqn{Y~I(A=0)*V*poly(U)+I(A=1)*V*poly(U)}      where
\eqn{poly(U)} consists  of \eqn{deg} orthogonal  polynoms of degrees  1 to
\eqn{deg}. If \eqn{nlev} is larger than or equal to two, then the model is
given                 by                the                 \code{formula}
\deqn{Y~I(A=0)*V*(poly(U)+step(U))+I(A=1)*V*(poly(U)+step(U))}       where
\eqn{poly(U)} consists  of \eqn{deg} orthogonal  polynoms of degrees  1 to
\eqn{deg} and  \eqn{step(U)} consits of \eqn{nlev}  indicator functions of
subsets of \eqn{[0,1]} of the form \eqn{\{x:x<\frac{k}{nlev}\}}.}
\value{Returns a \code{formula} describing the parametric model.}

\author{Antoine Chambaz [aut, cre]}





