pkgname <- "PoissonMultinomial"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
base::assign(".ExTimings", "PoissonMultinomial-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('PoissonMultinomial')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("dpmd")
### * dpmd

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: dpmd
### Title: Probability Mass Function of Poisson-Multinomial Distributions
### Aliases: dpmd

### ** Examples

pp <- matrix(c(.1, .1, .1, .7, .1, .3, .3, .3, .5, .2, .1, .2), nrow = 3, byrow = TRUE)
x <- matrix(c(0,0,1,2), nrow=1) 
x1 <- matrix(c(0,0,1,2,2,1,0,0),nrow=2,byrow=TRUE)

dpmd(pmat = pp)
dpmd(pmat = pp, xmat = x1)
dpmd(pmat = pp, xmat = x)

dpmd(pmat = pp, xmat = x, method = "NA" )
dpmd(pmat = pp, xmat = x1, method = "NA" )

dpmd(pmat = pp, method = "SIM", B = 1e3)
dpmd(pmat = pp, xmat = x, method = "SIM", B = 1e3)
dpmd(pmat = pp, xmat = x1, method = "SIM", B = 1e3)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("dpmd", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("ppmd")
### * ppmd

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: ppmd
### Title: Cumulative Distribution Function of Poisson-Multinomial
###   Distribution
### Aliases: ppmd

### ** Examples

pp <- matrix(c(.1, .1, .1, .7, .1, .3, .3, .3, .5, .2, .1, .2), nrow = 3, byrow = TRUE)
x <- matrix(c(3,2,1,3),nrow=1)
x1 <- matrix(c(0,0,1,2,2,1,0,0),nrow=2,byrow=TRUE)

ppmd(pmat = pp, xmat = x)
ppmd(pmat = pp, xmat = x1)

ppmd(pmat = pp, xmat = x, method = "NA")
ppmd(pmat = pp, xmat = x1, method = "NA")

ppmd(pmat = pp, xmat = x, method = "SIM", B = 1e3)
ppmd(pmat = pp, xmat = x1, method = "SIM", B = 1e3)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("ppmd", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("rpmd")
### * rpmd

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: rpmd
### Title: Poisson-Multinomial Distribution Random Number Generator
### Aliases: rpmd

### ** Examples

pp <- matrix(c(.1, .1, .1, .7, .1, .3, .3, .3, .5, .2, .1, .2), nrow = 3, byrow = TRUE)
 
rpmd(pmat = pp, s = 5)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("rpmd", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
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
