library("tsml.cara.rct")

tm.model <- formula(A~.)
learnQ <- formula(Y~I(as.integer(A)):(U+V)+I(as.integer(1-A)):(U+V))

ATE.param <- tsml.cara.rct(what="ATE",
                           flavor="parametric",
                           ninit=100,
                           by=50,
                           nmax=200,
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
