## code for assembling batch11 runs
## (as general as possible)
## module unload intel; module load r/3.1.1
## *could* assemble data into a giant array, but beta_mix_tab
##  already has it condensed the way we want -- in long format,
##  with replicates averaged to a mean value -- so we just have
##  to read everything and rbind() it together -- or equivalently
##  use plyr::ldply
library("betararef") ## not actually needed, I think
library("gtools") ## mixedsort
library("plyr")  ## ldply
run <- "betasim_batch11_"
outfn <- paste0(run,"allvals.RData")
lf <- mixedsort(list.files(pattern=paste0(run,".*\\.RData")))
lf2 <- gsub(paste0("(",run,"|\\.RData)"),"",lf)
dvals <- as.data.frame(do.call(rbind,strsplit(lf2,"_")))
getTab <- function(x) { load(x); return(beta_mix_tab) }
allvals <- ldply(lf,getTab,.progress="text")
save("allvals",file=outfn)
file.info(outfn)$size