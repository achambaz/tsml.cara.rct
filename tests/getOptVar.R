##
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
}
