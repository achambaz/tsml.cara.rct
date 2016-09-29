\name{getOptVar}
\alias{getOptVar}
\title{Computes the Optimal Variance Given a Parametric Model When Targeting the Average Treatment Effect}
\description{Computes the optimal variance when  targeting the Average Treatment Effect
(ATE). The optimal variance is defined as the variance associated with the
optimal treatment mechanism  within a given parametric  model of treatment
mechanisms. The computation is based on Monte Carlo simulation.}
\usage{getOptVar(tm.model = formula(A ~ 1), Gmin = 0.01, piV = c(1/2, 
    1/3, 1/6), family = c("beta", "gamma"), Qbar = Qbar1, Vbar = Vbar1, 
    n = 1000, tm.control = glm.control(maxit = 500), verbose = FALSE)}
\arguments{
  \item{tm.model}{A parametric  model \eqn{{\cal G}} of treatment  mechanisms. The procedure
targets  the optimal treatment  mechanism within  this model.  Defaults to
\code{formula(A~1)}.}
  \item{Gmin}{A  small positive  \code{numeric}, the  minimum value  of elements  of the
parametric  model \eqn{{\cal  G}}  of treatment  mechanisms (see  argument
\code{tm.model}). The maximum value is \code{1-Gmin}.}
  \item{piV}{Marginal distribution of \eqn{V}. Defaults to \code{c(1/2, 1/3, 1/6)}.}
  \item{family}{A \code{character}, either "beta" (default)  or "gamma", the nature of the
law of outcome.}
  \item{Qbar}{A   \code{function},  the   conditional  expectation   of   \eqn{Y}  given
\eqn{(A,W)}. Defaults to \code{Qbar1}.}
  \item{Vbar}{A   \code{function},   the   conditional   variance   of   \eqn{Y}   given
\eqn{(A,W)}. Defaults to \code{Vbar1}.}
  \item{n}{An \code{integer}, the  common sample size of the  two simulated data sets
used  to compute the  optimal treatment  mechanism and  associated optimal
variance.}
  \item{tm.control}{A \code{list}  of options to pass  to \code{glm} for the  targeting of the
treatment  mechanism  within the  model  defined  by argument  'tm.model'.
Defaults to \code{glm.control(maxit=500)}.}
  \item{verbose}{A \code{logical}  or an \code{integer}  indicating the level  of verbosity
(defaults to 'FALSE').}
}

\value{Returns  the  optimal  variance   associated  with  the  parametric  model
\eqn{{\cal G}} when targeting the Average Treatment Effect (ATE).}

\author{Antoine Chambaz [aut, cre]}



\seealso{getSample, getOptTm}
\examples{##
## testing if the optimal variance is non-negative
##

library("R.utils")
library("tsml.cara.rct")

psi.sd <- getOptVar(n=1e3,
                    tm.model=formula(A~1),
                    piV=c(1/2, 1/3, 1/6),
                    family="gamma",
                    Qbar=Qbar1,
                    Vbar=Vbar1)

if (psi.sd<0) {
  throw("Numeric 'psi.sd' should be non-negative.\n")
}}
