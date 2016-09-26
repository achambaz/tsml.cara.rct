library("tsml.cara.rct")

truth <- getSample(1e3,
                   tm=oneOne,
                   rule=NULL,
                   piV=c(1/2, 1/3, 1/6),
                   family="beta",
                   Qbar=Qbar2,
                   Vbar=Vbar2,
                   what="MOR")
truth

Gexploit <- function(nn) {
  if (nn <= 100) {
    exploit <- 0.2
  } else if (nn <= 125) {
    exploit <- 0.15
  } else if (nn <= 150) {
    exploit <- 0.1
  } else {
    exploit <- 0.05
  }
  return(exploit)
}
Gexplore <- Gexploit

MOR.lasso <- tsml.cara.rct(what="MOR",
                           flavor="parametric",
                           ninit=100,
                           by=50,
                           nmax=200,
                           learnQ=makeLearnQ,
                           tm.init=oneOne,
                           tm.ref=oneOne,
                           Gmin=0.10,
                           Gexploit=Gexploit,
                           Gexplore=Gexplore,
                           conf.level=0.95,
                           piV=c(1/2, 1/3, 1/6),
                           family="beta",
                           Qbar=Qbar2,
                           Vbar=Vbar2)                           
MOR.lasso

