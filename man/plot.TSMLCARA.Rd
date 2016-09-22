\name{plot.TSMLCARA}
\alias{plot.TSMLCARA}
\alias{plot}
\title{Plots a TSMLCARA Object}
\description{Plots a  TSMLCARA object.}
\usage{\method{plot}{TSMLCARA}(x, truth = NULL, regret = NULL, regret.mirror = FALSE, 
    lower.bound = TRUE, xlog = TRUE, colors = c("tomato", "olivedrab4"), 
    file = NULL, ...)}
\arguments{
  \item{x}{An object of class \code{TSMLCARA}.}
  \item{truth}{Defaults  to \code{NULL}.   When estimating  the Average  Treatment Effect
('what' equal to "ATE"), can be set to a \code{vector}, the true values of
both  the  targeted parameter  and  standard  deviation of  the  efficient
influence  curve under  the optimal  treatment mechanism  given a  working
model,  as output  by 'getOptVar'.   When  estimating the  Mean under  the
Optimal treatment  Rule ('what' equal  to 'MOR'), can  be set to  either a
\code{vector} or  a list of two  vectors.  In the former  case, the vector
gives  the  true  values  of  both the  targeted  parameter  and  standard
deviation of  the efficient  influence curve  under the  optimal treatment
rule. In the latter case, the list  contains the same vector as before and
a  second  vector   giving  the  values  of   the  targeted  data-adaptive
parameters.}
  \item{regret}{Defauls to \code{NULL}. When estimating the optimal treatment rule ('what'
equal to 'MOR'), can  be set to a vector or a list  of two vectors. In the
former case,  the vectors gives  the true values  of the empirical  or the
counterfactual regrets. In the latter case, the list contains both.}
  \item{regret.mirror}{Should the regret(s) be mirrored? If  \code{regret} is given, should it be
multiplied by (-1)? (defaults to \code{FALSE}).}
  \item{lower.bound}{Defaults to \code{TRUE}. Should the  confidence lower bound on the regrets
be printed or a confidence interval?}
  \item{xlog}{Defaults  to \code{TRUE}.  Should x-axis be on a log-scale?}
  \item{colors}{A length 2 vector of  colors. Defauls to \code{c("tomato", "olivedrab4")},
colors used for plotting.}
  \item{file}{A \code{character}.   Defaults to NULL,  and if a \code{character}  then a
PDF plot is produced, with name "'file'.pdf".}
  \item{\dots}{Not used.}
}



\author{Antoine Chambaz [aut, cre]}



\seealso{tsml.cara.rct}
\examples{FALSE}
