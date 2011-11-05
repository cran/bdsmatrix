attach(".RData")
dyn.load("bdsmatrix.so")
options(na.action='na.exclude', contrasts=c('contr.treatment', 'contr.poly'))
date()

aeq <- function(x,y, ...) all.equal(as.vector(x), as.vector(y),...)
