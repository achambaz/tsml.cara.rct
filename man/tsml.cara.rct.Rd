\name{tsml.cara.rct}
\alias{tsml.cara.rct}
\alias{tsml.cara.rct}
\title{Targeted  Minimum Loss  Covariate-Adjusted  Response-Adaptive RCT Design and Statistical Analysis}
\description{Simulates a targeted minimum loss covariate-adjusted response-adaptive RCT
design and statistical analysis.}
\usage{tsml.cara.rct(what = c("ATE", "MOR"), 
    flavor = c("parametric", "lasso"), ninit = 50, by = 25, nmax = 500,
    tm.init = oneOne, tm.ref = oneOne, learnQ,
    tm.model = formula(A ~ 1), tm.control = glm.control(maxit = 500), 
    Gmin = 0.01, Gexploit = Gmin, Gexplore = 0.01, Qmin = 0.01, 
    conf.level = 0.95, verbose = FALSE, piV = c(1/2, 1/3, 1/6), 
    family = c("beta", "gamma"), Qbar = Qbar1, Vbar = Vbar1, 
    Bn = 1e+05, slice.by = 1e+05)}
\arguments{
  \item{what}{A  \code{character}  indicating the  parameter  of  interest to  estimate.
Either "ATE" for the Average  Treatment Effect, the difference between the
means under  '\eqn{do(A=1)}' and  '\eqn{do(A=0)}', or  "MOR" for  the Mean
under the Optimal treatment Rule '\eqn{do(A=r(W))}'.}
  \item{flavor}{A \code{character}  indicating the  'flavor' of the  procedure. }
  \item{ninit}{An     \code{integer},    number    of     subjects    to     sample    at
initialization. Defaults to 50.}
  \item{by}{An  \code{integer},  number of  subjects  to  sample  at each  step  after
initialization. Defaults to 25.}
  \item{nmax}{An  \code{integer},  maximum  number  of  subjects to  sample  during  the
trial. Must be larger than 'ninit+by'. Defaults to 500.}
  \item{tm.init}{A  \code{function}  describing  the  initial  treatment  mechanism  to  be
employed.   Defaults  to  the  balanced  (1:1)  treatment  mechanism,  ie,
\code{function} \code{\link{oneOne}}.}
  \item{tm.ref}{A  \code{function}  describing the  reference  treatment  mechanism to  be
employed.   Defaults  to  the  balanced  (1:1)  treatment  mechanism,  ie,
\code{function} \code{\link{oneOne}}.}
  \item{learnQ}{A  model  \eqn{{\cal Q}}  of  conditional  expectations of  \eqn{Y}  given
\eqn{(A,W)}  for  both  flavors  'parametric'  and  'lasso',  given  as  a
\code{formula}  or  a  \code{function}  outputing  formulas.  Defaults  to
\code{formula(Y~1)} for flavors  'parametric' and 'lasso'.}
  \item{tm.model}{A parametric model \eqn{{\cal G}}  of treatment mechanisms, used only when
'what'  equals  "ATE".   The   procedure  targets  the  optimal  treatment
mechanism within this model.  Defaults to \code{formula(A~1)}.}
  \item{tm.control}{A  \code{list} of  options for  the targeting  of the  treatment mechanism
within the  model defined by  argument 'tm.model'.  Used only  when 'what'
equals "ATE", it defaults to \code{glm.control(maxit=500)}.}
  \item{Gmin}{A  small positive  \code{numeric}, with  default value  \code{1e-2}.  When
\code{what} equals 'ATE', it is  the minimum value of elements of the
parametric  model \eqn{{\cal  G}}  of treatment  mechanisms (see  argument
\code{tm.model}).  The  maximum value is  \code{1-Gmin}.  When \code{what}
equals 'MOR', it  is the minimum value of  the conditional probability
of \eqn{A=r_n(W)} given \eqn{W}.}
  \item{Gexploit}{A small positive  \code{numeric}, with default value  that of \code{Gmin},
or a  function of sample  size giving such  small numbers, only  used when
\code{what} equals "MOR", in conjunction with \code{Gexplore}.}
  \item{Gexplore}{Either a small positive \code{numeric}, with default value \code{1e-2}, or
a  function of  sample  size giving  such small  numbers,  only used  when
\code{what} equals "MOR", in conjunction with \code{Gexploit}.}
  \item{Qmin}{A  small positive  \code{numeric}, the  minimum value  of  scaled outcomes
\eqn{Y}. The maximum value is \code{1-Qmin}.}
  \item{conf.level}{A  \code{numeric},  the  confidence  level of  the  resulting  confidence
interval.}
  \item{verbose}{A \code{logical}  or an \code{integer}  indicating the level  of verbosity
(defaults to 'FALSE').}
  \item{piV}{Marginal distribution of \eqn{V}. Defaults to \code{c(1/2, 1/3, 1/6)}.}
  \item{family}{A \code{character}, either "beta" (default)  or "gamma", the nature of the
law of outcome.}
  \item{Qbar}{A   \code{function},  the   conditional  expectation   of   \eqn{Y}  given
\eqn{(A,W)}. Defaults to \code{Qbar1}.}
  \item{Vbar}{A   \code{function},   the   conditional   variance   of   \eqn{Y}   given
\eqn{(A,W)}. Defaults to \code{Vbar1}.}
  \item{Bn}{An \code{integer}, the sample size used  to estimate the true value of the
data-adaptive parameter at  each step of the procedure  when 'what' equals
'MOR'. Defaults to 1e5.}
  \item{slice.by}{An \code{integer}. If it is smaller than argument 'n' of 'getSample', then
the simulation  is decomposed  into 'n\%/\%slice.by' smaller  simulations of
'slice.by' observations and one of 'n\%\%slice.by' observations. Defaults to
1e5 (hence, no decomposition if 'n' smaller than 4e5). Mainly for internal
use.}
}
\details{Defines  a  lower-bound   on  the  conditional  probability  of
\eqn{do(A=1-r_n(W))} given \eqn{W}.}
\value{Returns a  \code{TSMLCARA} object  which summarizes the  TSMLCARA undertaken
procedure.}
\references{Chambaz, van der Laan,  Zheng, Chapter  16, Modern Adaptive Randomized Clinical  Trials: Statistical, Operational, and Regulatory    Aspects,    by    A.    Sverdlov    (CRC    Press,    2015). }
\author{Antoine Chambaz [aut, cre]}



\seealso{update, targetPsi, getSample}
\examples{
##
log <- Arguments$getVerbose(-1, timestamp=TRUE)
set.seed(12345)

## ########################
## AVERAGE TREATMENT EFFECT
## ########################
tm.model <- formula(A~.)
psi.sd <- sqrt(getOptVar(n=1e5,
                         tm.model=tm.model,
                         piV=c(1/2, 1/3, 1/6),
                         family="gamma",
                         Qbar=Qbar1,
                         Vbar=Vbar1))
truth <- c(psi=91/72, psi.sd=psi.sd)

## parametric example
learnQ <- formula(Y~I(as.integer(A)):(U+V)+I(as.integer(1-A)):(U+V))
ATE.param <- tsml.cara.rct(what="ATE",
                           flavor="parametric",
                           ninit=200,
                           by=100,
                           nmax=400,
                           tm.init=oneOne,
                           tm.ref=oneOne,
                           learnQ=learnQ,
                           tm.model=tm.model,
                           conf.level=0.95,
                           piV=c(1/2, 1/3, 1/6),
                           family="gamma",
                           Qbar=Qbar1,
                           Vbar=Vbar1)
ATE.param
## Not run:
plot(ATE.param, truth=truth)
## End(**Not run**)

## See the vignette for more examples...
}
