pkgname <- "PMD"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('PMD')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("dpmd")
### * dpmd

flush(stderr()); flush(stdout())

### Name: dpmd
### Title: Probability Mass of Poisson-Multinomial Distributions
### Aliases: dpmd

### ** Examples


pp=matrix(c(.1, .1, .1, .7, .1, .3, .3, .3, .5, .2, .1, .2), nrow=3, byrow=TRUE)
dpmd(pp)
dpmd(pp,"simulation",B=10^3)
dpmd(pp,"NA by demands", vec = c(0,0,1,2))
dpmd(pp,"simulation by demands", vec = c(0,0,1,2), B=10^3)



cleanEx()
nameEx("pmatrix")
### * pmatrix

flush(stderr()); flush(stdout())

### Name: pmatrix
### Title: pmatrix
### Aliases: pmatrix

### ** Examples

pp = pmatrix(2,2)
pp



cleanEx()
nameEx("ppmd")
### * ppmd

flush(stderr()); flush(stdout())

### Name: ppmd
### Title: cumulative mass function of PMN
### Aliases: ppmd

### ** Examples

pp=matrix(c(.1, .1, .1, .7, .1, .3, .3, .3, .5, .2, .1, .2), nrow=3, byrow=TRUE)
ppmd(pp,c(3,2,1,3))



cleanEx()
nameEx("rpmd")
### * rpmd

flush(stderr()); flush(stdout())

### Name: rpmd
### Title: generate random number from PMD
### Aliases: rpmd

### ** Examples

pp=matrix(c(.1, .1, .1, .7, .1, .3, .3, .3, .5, .2, .1, .2), nrow=3, byrow=TRUE)
rpmd(pp)



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
