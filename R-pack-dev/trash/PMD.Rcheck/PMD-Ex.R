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
### Title: Probability Mass Function of Poisson-Multinomial Distributions
### Aliases: dpmd

### ** Examples

pp=matrix(c(.1, .1, .1, .7, .1, .3, .3, .3, .5, .2, .1, .2), nrow = 3, byrow = TRUE)

dpmd(pmat = pp)
dpmd(pmat = pp, method = "SIM-ALL", B = 1e3)
dpmd(pmat = pp, x = c(0,0,1,2), method = "NA" )
dpmd(pmat = pp, x = c(0,0,1,2), method = "SIM", B = 1e3)




cleanEx()
nameEx("ppmd")
### * ppmd

flush(stderr()); flush(stdout())

### Name: ppmd
### Title: Cumulative Distribution Function of Poisson-Multinomial
###   Distribution
### Aliases: ppmd

### ** Examples

pp=matrix(c(.1, .1, .1, .7, .1, .3, .3, .3, .5, .2, .1, .2), nrow = 3, byrow = TRUE)

ppmd(pmat = pp, x = c(3,2,1,3))
ppmd(pmat = pp, x = c(3,2,1,3), method = "NA")
ppmd(pmat = pp, x = c(3,2,1,3), method = "SIM-ALL", B = 1e3)



cleanEx()
nameEx("rpmd")
### * rpmd

flush(stderr()); flush(stdout())

### Name: rpmd
### Title: Poisson-Multinomial Distribution Random Number Generator
### Aliases: rpmd

### ** Examples

pp=matrix(c(.1, .1, .1, .7, .1, .3, .3, .3, .5, .2, .1, .2), nrow = 3, byrow = TRUE)
 
rpmd(pmat = pp, s = 5)




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
