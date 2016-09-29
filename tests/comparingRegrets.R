##
## comparing  the empirical  cumulative  pseudo-regrets under  optimal TR  and
## balanced TR
##

library("R.utils")
library("tsml.cara.rct")

piV <- c(1/2, 1/3, 1/6)
family <- "beta"

Qbar2.adapt <- function(A,W) {
  AW <- cbind(A=A, W)
  Qbar2(AW)
}
optRule <- tsml.cara.rct:::ruleQbar(Qbar2.adapt)

obs.bTR <- getSample(1000,
                     tm=oneOne,
                     piV=piV, family=family,
                     Qbar=Qbar2, Vbar=Vbar2)

obs.oTR <- getSample(1000,
                     tm=optRule,
                     piV=piV, family=family,
                     Qbar=Qbar2, Vbar=Vbar2)

## balanced TR
W <- tsml.cara.rct:::extractAW(obs.bTR)
rA <- as.numeric(optRule(W))
eregret.bTR <- mean(obs.bTR[, "Y"] - Qbar2(cbind(A=rA, W)))
## ## what we lose by not adapting

## optimal TR
W <- tsml.cara.rct:::extractW(obs.oTR)
rA <- as.numeric(optRule(W))
eregret.oTR <- mean(obs.oTR[, "Y"] - Qbar2(cbind(A=rA, W)))
## ## what we lose despite being optimal (goes to zero quickly!)
