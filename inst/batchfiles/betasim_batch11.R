## see betasim_batch*gen
BATCHNUM <- 11
setwd("/work/bolker/projects/beta")  ## for SHARCnet
batchfile <- paste0("betasim_batch",BATCHNUM,"_ABUND.RData")

library(betararef)
library(reshape2) ## for melt
library(vegan)

## p
abundvec <- 10^(1+(ABUND-1)*0.1)
## batch runs ABUND from 1 to 41 (same as below)
## abundvec <- round(10^seq(from=1,to=5,by=0.1))  ## how many individuals?
stvec <- c(5,10,20)        ## number of sites
pmixcommonvec <- pmixrarevec <- c(1,0.5,0.1,0) ## redistribution/mixing across reefs
spcurvevec <- c(0.5,0.2,0.05,0.02)               ## rank-abundance curve
rep <- 400
gammadiv <- 20
verbose <- TRUE
distvec1 <- c("bray","morisita","horn")
distvec2 <- c("pairwise","centroid")
beta_mix <- array(NA,
                  dim=c(
                  length(distvec1),
                  length(distvec2),
                  length(abundvec),
                  length(stvec),
                  length(pmixcommonvec),
                  length(pmixrarevec),
                  length(spcurvevec),
                  rep),
                  dimnames=list(
                  dist1=distvec1,
                  dist2=distvec2,
                  abund=abundvec,
                  st=stvec,
                  pmixcommon=pmixcommonvec,
                  pmixrare=pmixrarevec,
                  spcurve=spcurvevec,
                  rep=seq(rep)))
set.seed(101)
for (d1 in seq_along(distvec1)) {
    for (d2 in seq_along(distvec2)) {
        for(i in seq_along(abundvec)){
            for (j in seq_along(stvec)) {
                for(k in seq_along(pmixcommonvec)) {
                    for(l in seq_along(pmixrarevec)) {
                        for(m in seq_along(spcurvevec)) {
                            ## verbose: don't print *every* rep ...
                            if (verbose) cat(d1,"/",length(distvec1)," ",
                                             d2,"/",length(distvec2)," ",
                                             i,"/",length(abundvec)," ",
                                             j,"/",length(stvec)," ",
                                             k,"/",length(pmixcommonvec)," ",
                                             l,"/",length(pmixrarevec)," ",
                                             m,"/",length(spcurvevec)," ",
                                             "\n",
                                             sep="")
                            for (n in seq(rep)) {
                                bb <- betasim(n.abund=2,
                                              p.abund=spcurvevec[m],
                                              spcat=gammadiv/stvec[j],
                                              n.site=stvec[j],
                                              n.indiv.site=abundvec[i],
                                              rand="poisson",
                                              rarefy=1,
                                              p.mix=c(pmixcommonvec[k],pmixrarevec[l]))
                                beta_mix[d1,d2,i,j,k,l,m,n] <- mean(calcbeta(bb,method=distvec1[d1],distances=distvec2[d2]))

                            } ## reps
                        }  ## spcurve
                    } ## mix1 
                } ## mix2
                save("beta_mix",file=batchfile) ## checkpoint
            } ## nsites
        }  ## abundance
    } ## distvec2
} ## distvec1

sessinf <- sessionInfo()
save("beta_mix",file=batchfile)
repdim <- which(names(dimnames(beta_mix))!="rep")
beta_mix_mean <- apply(beta_mix,repdim,mean,na.rm=TRUE)  ## collapse over reps
beta_mix_tab <- melt(beta_mix_mean)
save("beta_mix_tab","beta_mix_mean","beta_mix","sessinf",file=batchfile)
