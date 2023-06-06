install.packages(
    "~/Desktop/PANPRS",
    type = "source",
    repos = NULL,
    force = TRUE
)

library(devtools)
load_all("~/Desktop/PANPRS/")

library("PANPRS")
data("summaryZ")
data("Nvec")
data("plinkLD")
data("funcIndex")

output <- gsfPEN(
    summaryZ = summaryZ,
    Nvec = Nvec,
    plinkLD = plinkLD,
    funcIndex = funcIndex,
    numfunc = ncol(funcIndex)
)
