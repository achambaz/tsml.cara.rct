library("tsml.cara.rct")

truth <- getSample(1e5,
                   tm=oneOne,
                   rule=NULL,
                   piV=c(1/2, 1/3, 1/6),
                   family="beta",
                   Qbar=Qbar2,
                   Vbar=Vbar2,
                   what="MOR")
truth

Gexploit <- function(nn) {
  if (nn <= 200) {
    exploit <- 0.2
  } else if (nn <= 500) {
    exploit <- 0.15
  } else if (nn <= 1000) {
    exploit <- 0.1
  } else {
    exploit <- 0.05
  }
  return(exploit)
}

Gexplore <- function(nn) {
  if (nn <= 200) {
    explore <- 0.2
  } else if (nn <= 500) {
    explore <- 0.15
  } else if (nn <= 1000) {
    explore <- 0.1
  } else {
    explore <- 0.05
  }
  return(explore)
}


MOR.lasso <- tsml.cara.rct(what="MOR",
                           flavor="lasso",
                           ninit=200,
                           by=200,
                           nmax=2000,
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

