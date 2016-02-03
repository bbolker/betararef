library(betararef)
packageVersion("betararef")
require(plyr)
require(vegan)
require(reshape2)
fn <- "rareftest-batch17.RData"
descr  <- c(short="raref_power",
            long="POWER 4: pmix*sizediff*popsize*popsizediff; rarefy within permutations, median P (corrected 2)")
nsim <- 250

## 10 indiv/site (default) mean 20 patches/ttt, 2 treatments
## 50 raref, 200 permutations each
pmixdiffvec <- seq(0,0.5,by=0.1)
sizediffvec <- c(0,20,30)  ## pick sizes that allow 
                           ## even multiples in each treatment
                           ## (20/20, 10/30, 5/35)
popsizevec <- c(10,20)  ## only popsize 10/20
popsizediffvec <- c(0,0.05,0.1,0.5)
tmpf <- function(...) {
    sx <- betararef:::simComm(...)
    tt <- try(raref_test2(sx[,-1],sx[,1],
                          ptype="raref_perm_medianp")$p.value)
    if (is(tt,"try-error")) NA else tt
}
set.seed(102)
power2arr <- array(NA,dim=c(nsim,
                      length(pmixdiffvec),
                      length(sizediffvec),
                      length(popsizevec),
                      length(popsizediffvec)),
                   dimnames=list(rep=seq(nsim),
                   pmixdiff=pmixdiffvec,
                   sizediff=sizediffvec,
                   popsize=popsizevec,
                   popsizediff=popsizediffvec))

for (i in seq_along(pmixdiffvec)) {
    cat("pmix",i,"\n")
    p.mixes <- 0.5+c(-1,1)*pmixdiffvec[i]/2
    for (j in seq_along(sizediffvec)) {
        cat("sizediff",j,"\n")
        n.sites <- 20+c(-1,1)*sizediffvec[j]/2
        for (k in seq_along(popsizevec)) {
            for (l in seq_along(popsizediffvec)) {
                cat(i,j,k,l,"\n")
                popsize <- popsizevec[k]
                pdiff <- popsizediffvec[l]
                p.sizes <- round(c(popsize*(1-pdiff/2),popsize*(1+pdiff/2)))
                power2arr[,i,j,k,l] <-
                    raply(nsim,tmpf(n.site=n.sites,
                                    p.mix=p.mixes,
                                    n.indiv.site=p.sizes,
                                    totsp=c(20,20)))
                if (FALSE) {
                    ## testing; this fails on seed=102, but haven't
                    ## figured out why right now -- just patched over
                    ## with a try() in tmpf() above.
                    for (i in 1:50) {
                        cat(i,"\n")
                        set.seed(100+i)
                        tmpf(n.site=n.sites,
                             p.mix=p.mixes,
                             n.indiv.site=p.sizes)
                    }
                }
            }
            save("power2arr",file=fn)                        
        }
        
    }
}
    
save("power2arr",file=fn)

sessionInfo()
