\name{getSample}
\alias{getSample}
\title{Generates  Data}
\description{Generates  a run  of simulated  observations of  the  form \eqn{(W,G,A,Y,Y1,Y0)}
under  the  specified treatment  mechanism.  }
\usage{getSample(n, tm = oneOne, rule = NULL, piV = c(1/2, 1/3, 1/6), 
    family = c("beta", "gamma"), Qbar = Qbar1, Vbar = Vbar1, 
    what, slice.by = n)}
\arguments{
  \item{n}{An \code{integer}, the number of observations to be generated.}
  \item{tm}{A  \code{function}  describing the  treatment  mechanism  to be  employed.
Defaults to  the balanced  (1:1) treatment mechanism,  ie, \code{function}
\code{\link{oneOne}}. The \code{G} column equals the vector \code{tm(W)}.}
  \item{rule}{Either 'NULL' or a \code{function} describing the rule to be employed when
\code{what}  equals "MOR".  In  that case,  it must  be a  function of
\eqn{W} with values in \eqn{\{0,1\}}.}
  \item{piV}{Marginal distribution of \eqn{V}. Defaults to \code{c(1/2, 1/3, 1/6)}.}
  \item{family}{A \code{character}, either "beta" (default)  or "gamma", the nature of the
law of outcome.}
  \item{Qbar}{A   \code{function},  the   conditional  expectation   of   \eqn{Y}  given
\eqn{(A,W)}. Defaults to \code{Qbar1}.}
  \item{Vbar}{A   \code{function},   the   conditional   variance   of   \eqn{Y}   given
\eqn{(A,W)}. Defaults to \code{Vbar1}.}
  \item{what}{A \code{character}. If it  is not missing then it must  be equal to either
"ATE"  (Average  Treatment  Effect)  or  "MOR"  (Mean  under  the  Optimal
treatment Rule).  In  that case, ONLY estimates of the  true parameter and
optimal  standard deviation  under the  specified simulation  scheme (and,
possibly, treatment mechanism) are computed and returned.}
  \item{slice.by}{An \code{integer}. If it is smaller than argument 'n', then the simulation
is  decomposed  into  'n\%/\%slice.by'  smaller  simulations  of  'slice.by'
observations and one  of 'n\%\%slice.by' observations. Defaults  to 'n' (ie,
no decomposition). Mainly for internal use.
}
}
\details{By default,  implements the  same simulation  scheme as  in
Chambaz, van der Laan, Scand.  J. Stat., 41(1):104--140 (2014).}
\value{If   \code{truth}  is   'FALSE'  then   returns  a   \code{data.frame}  of
observations   with  columns   \code{Y1}  and   \code{Y0}  (counterfactual
outcomes),  \code{Y}  (actual  outcome),  \code{A}  (assigned  treatment),
\code{G} (the  conditional probability that  A=1 given W), the  rest being
interpreted  as   covariates.   Otherwise,  returns  the   estimated  true
parameter value  and standard deviation  of the efficient  influence curve
under the specified  simulation scheme and treatment  mechanism, if 'what'
equals "ATE", or the standard  deviationS of the efficient influence curve
under the specified simulation scheme and (i) treatment mechanism and (ii)
the optimal treatment rule, if 'what' equals "MOR".}
\references{Chambaz,  van der Laan, Scand.  J.  Stat., 41(1):104--140
(2014).}
\author{Antoine Chambaz [aut, cre]}




\examples{
## Setting the verbosity parameter
library(R.utils)
log <- Arguments$getVerbose(-8, timestamp=TRUE)

##
obs <- getSample(10)
}
