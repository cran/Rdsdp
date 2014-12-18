# Test Rdsdp for sdpa data sets from $DSDPROOT/bin
library(Rdsdp)

datapath="../Rdsdp/extdata/"
tol=0.000001

target=c(8.999996,-17.784626,-40.819012,-226.1573)
#names(target) = c("truss1", "control1", "vibra1", "mcp100")

truss1 = dsdp.readsdpa(paste(datapath, "truss1.dat-s", sep="/"))
all.equal.numeric(target[1], truss1[[3]][[2]], tolerance=tol)

control1 = dsdp.readsdpa(paste(datapath, "control1.dat-s", sep="/"))
all.equal.numeric(target[2], control1[[3]][[2]], tolerance=tol)

vibra1 = dsdp.readsdpa(paste(datapath, "vibra1.dat-s", sep="/"))
all.equal.numeric(target[3], vibra1[[3]][[2]], tolerance=tol)

mcp100 = dsdp.readsdpa(paste(datapath, "mcp100.dat-s", sep="/"))
all.equal.numeric(target[4], mcp100[[3]][[2]], tolerance=tol*100)
