##' Find target sizes for rarefaction of samples for multiple treatments
##' down to a common size.
##'
##' @importFrom ggplot2 cut_number
##' @param nList a list of numeric vectors of patch sizes (total abundances) for each treatment
##' @param method (character): 'smallest' to use the smallest abundance found in any patch in any treatment; 'rankMatch' to match patches by rank within treatment and use the smallest abundance within a rank category; 'randomMatch' to match patches into groups at random. The idea of the default 'rankMatch'
##' is that by matching patches by rank we rarefy as
##' little as necessary to get equal patch sizes in the two
##' groups -- i.e. the minimum in both groups is the same, etc.
##' (I can't prove that this minimizes the required
##' amount of rarefaction, and I'm not even sure it's always
##' true, but in any case it seems to be a reasonable approach.)
##' @param maxit number of iterations to try for randomized algorithms
##' @param debug debugging output?
##' @return a list of vectors of target sizes corresponding to the patches described in \code{nList}
##' @examples
##' ## make sure that the order of the patches (i.e. order of number
##' ## of fish in each patch) doesn't matter
##' nList <- list(1:10,10:1)
##' findTargets(nList)  ## unchanged
##' ## three treatments
##' nList <- list(5:10,9:4,c(2,7,8,1,20,8))
##' f1 <- findTargets(nList)
##' ## unequal numbers of patches
##' nList <- list(5:10,1:10)
##' set.seed(101)
##' (f2 <- findTargets(nList))
##' @importFrom stats anova lm median pbinom
##' @importFrom vegan permutest
##' @export
findTargets <- function(nList,method=c("rankMatch","randomMatch","smallest"),
                        maxit=1000, debug=FALSE) {
    method <- match.arg(method)
    if (method=="smallest") {
        ## target for ALL patches is smallest abundance in entire group
        smallVal <- min(unlist(nList))
        return(lapply(nList,
                      function(x) { x[] <- smallVal; x}))
    }
    if (method=="rankMatch") {
        nTreat <- length(nList)
        if (length(unique(nPatches <- sapply(nList,length)))>1) {
            ## unequal patch size case
            minPatch <- min(nPatches)
            ## divide patches in each treatment evenly into categories
            ## corresponding to the number of patches in the min-patch-ttt
            patchCat <- lapply(nList,
               function(x) {
                   n <- minPatch
                   if (length(x)==n) return(seq_along(x))
                   ## FIX: cut seq_along(x), not x
                   as.numeric(cut_number(seq_along(x),n=minPatch))
               })
        } else {
            ## if all the same, the categories are just 1:npatch
            minPatch <- nPatches[1]
            patchCat <- replicate(nTreat,seq(minPatch),simplify=FALSE)
        }
        sortList <- lapply(nList,sort)
        ## break ties randomly
        rankList <- lapply(nList,rank,ties.method="random")
        ## split each sorted list into patch categories
        splitSortList <- mapply(split,sortList,patchCat,SIMPLIFY=FALSE)
        ## set up structure (same shape) for results
        splitTargetList <- splitSortList
        for (i in seq(minPatch)) {
            list1 <- lapply(splitSortList,"[[",i)
            ## find which treatment has min AVERAGE patch pop size
            ##   in this category
            catMean <- sapply(list1,mean)
            whichMinTrt <- which.min(catMean)
            ## min-mean-treatment stays the same
            ##  (target=original)
            splitTargetList[[whichMinTrt]][[i]] <-
                splitSortList[[whichMinTrt]][[i]]
            ## other treatments get a sample with replacement
            ## from min-mean-treatment pop sizes
            for (j in seq(nTreat)[-whichMinTrt]) {
                t <- splitSortList[[whichMinTrt]][[i]]
                if (debug) cat(i,j,t,"\n")
                tlen <- length(splitTargetList[[j]][[i]])
                if (length(t)==1) {
                    ## sometimes there is one sample > the
                    ## corresponding min-category size ...
                    s <- pmin(splitSortList[[j]][[i]],rep(t,tlen))
                } else {
                    repeat {
                        s <- sample(t,
                               size=tlen,
                               replace=TRUE)
                        ## keep trying until all samples are <=
                        ##  the corresponding original pop size
                        if (all(s<=splitSortList[[j]][[i]])) break
                    }
                }
                if (debug) cat("...",s,"\n")
                splitTargetList[[j]][[i]] <- s
                ## cat(i,j,"\n")
                stopifnot(length(unlist(splitTargetList))==length(unlist(splitSortList)))
            }
        }
        ## collapse lists
        tVals <- lapply(splitTargetList,unlist)
        ## PREVIOUSLY: find min val for each set
        ## tVals <- do.call(mapply,c(list(FUN=min),sortList))
        ## assign them to the corresponding element of each set
        tList <- mapply(function(x,y) y[x],rankList,tVals,SIMPLIFY=FALSE)
        return(tList)
    }
    if (method=="randomMatch") {
        ## match patches by rank
        if (length(unique(sapply(nList,length)))>1)
            stop("randomMatch not yet implemented for unequal numbers of patches")
        ## what is inverse permutation?
        permList <- lapply(nList,function(x) sample(length(x)))
        orderList <- lapply(permList,order)
        ## find min val for each set
        tVals <- do.call(mapply,c(list(FUN=min),permList))
        ## assign them to the corresponding element of each set
        tList <- lapply(orderList,function(x) tVals[x])
        return(tList)
    }
    stop("unimplemented method")
}

    
## rarefies species abundance vector single treatment down to 'n.target' individuals
sampleRow <- function(xmat,n.target,spnames=names(xmat)) {
  spvec <- rep(spnames,xmat)  ## expand species counts to a vector of individuals for each row
  if (n.target==length(xmat)) return(xmat)  ## don't really rarefy
  ## subsample and collapse back down to a species count
  table(factor(sample(spvec,replace=FALSE,size=n.target),
               levels=spnames))
}

##' Individual-based resampling/rarefaction of a community matrix
##' 
##' @param x a species abundance matrix to rarefy (rows=sites, columns=species)
##' @param n.target rarefied size for each sample
##' @param debug (logical) debugging flag
##' @return a matrix of the same structure as the original
##' @export
sampleMat <- function(x,n.target,debug=FALSE) {
    which(rowSums(x)<n.target)
    for (i in 1:nrow(x))
        x[i,] <- sampleRow(x[i,],n.target[i],colnames(x))
    x
}

sampleList <- function(matList,targetList) {
    mapply(sampleMat,matList,targetList,SIMPLIFY=FALSE)
}

## single rarefaction: take list of matrices and target sizes,
##  return a single combined community matrix
rarefy1 <- function(matList,tList) {
  as.matrix(do.call(rbind,sampleList(matList,tList)))
}

##' Calculate beta-diversity statistic for a single realization
##' 
##' @param m community matrix
##' @param ttt treatment/grouping factor
##' @param method beta diversity metric
##' @param binary compute beta diversity based on presence/absence?
##' @param stat beta diversity statistic to store: \code{Fstat}=F-statistic from \code{\link{anova.betadisper}}; \code{Fstat_perm}=F-statistic from \code{\link{permutest.betadisper}}; \code{pval_perm}=p-value from \code{\link{permutest.betadisper}}; \code{pval_Fstat_perm}=F- and p-values as above
##' @export
calcbeta_stat <- function(m,ttt,method="jaccard",binary=TRUE,
                     stat=c("Fstat","Fstat_perm","pval_perm","pval_Fstat_perm")) {
    stat <- match.arg(stat)
    b <- betadisper(vegdist(m,method=method,binary=binary,
                            bias.adjust=TRUE),ttt)
    if (stat=="sum") return(sum(b$distances^2))
    else if (stat=="Fstat") {
        ## return observed F statistic
        return(anova(b)["Groups","F value"])
    } else if (stat=="Fstat_perm") {
        ## return vector: first element is observed,
        ##  remainder are permuted F-statistics
        pp <- permutest(b)
        vals <- c(pp$statistic[1],pp$permutations[,1])
        return(vals)
    } else if (stat=="pval_perm") {
        p <- permutest(b)$tab[["Pr(>F)"]][1]
        return(p)
    } else if (stat=="pval_Fstat_perm") {
        pp <- permutest(b)
        return(c(pval=pp$tab[["Pr(>F)"]][1],
                 F=pp$tab[["F"]][1]))
    }
}

## permute species matrix (not currently used)
perm_spMat <- function(spMat) {
  spMat[sample(nrow(spMat)),]
}


##' rarefy once (maybe) and spit out F statistic(s)
##'
##' @rdname simfun2
simfun0 <- function(matList,tList,ttt, rarefy=TRUE,bstat="Fstat") {
    r1 <- if (rarefy) rarefy1(matList,tList) else do.call(rbind,matList)
    cc <- calcbeta_stat(r1,ttt,stat=bstat)
}

##'  return vector (or matrix) of F statistics
##'
##' @rdname simfun2
##' @param matList list of community matrices
##' @param tList list of target sizes for each community
##' @param ttt treatment vector
simfun1 <- function(matList, tList, ttt, bstat, rarefy=TRUE,
                    .progress="text",
                    nsim=100) {
    simvec <- raply(nsim,simfun0(matList,tList, ttt,
                                 bstat=bstat,
                                 rarefy=rarefy),.progress=.progress)
    simvec
}

##' simulate multiple communities with different beta/abundance/etc.
##' characteristics
##'
##' @inheritParams betasim
##' @param ncomm number of communities
##' @param commlabs labels for communities/treatments
##' @param retval return communities as a single metacommunity matrix or as a list of community matrices?
##' @param seed random-number seed
##' @importFrom plyr rbind.fill
##' @examples
##' plot(simComm(totsp=c(30,30),p.mix=c(0.2,0.7)))
##' @export
simComm <- function(n.indiv.site=rep(10,ncomm),
                    n.site=rep(3,ncomm),
                    spcat=rep(1,ncomm),
                    ncomm=2,
                    p.mix=rep(0.5,ncomm),
                    n.abund=rep(1,ncomm),
                    totsp=rep(3,ncomm),
                    rand="poisson",
                    commlabs=letters[1:ncomm],
                    retval=c("combined","list"),
                    seed=NULL) {
    if (!is.null(seed)) set.seed(seed)
    retval <- match.arg(retval)
    ## (allow combination of data frames w/ different numbers of
    ##  columns, filling in NA as appropriate)
    ## generate all communities
    ## FIXME: are there other parameters we want to vary between
    ## communities/treatments?
    L <- mapply(betasim,
                n.site      =as.list(n.site),
                n.indiv.site=as.list(n.indiv.site),
                ## spcat       =as.list(spcat),
                p.mix       =as.list(p.mix),
                n.abund     =as.list(n.abund),
                totsp       =as.list(totsp),
                MoreArgs=list(rand=rand),SIMPLIFY=FALSE)
    ## treatment ID vector
    names(L) <- commlabs
    ttt <- rep(commlabs,n.site)
    ## put it all together
    if (retval=="combined") {
        ## return as data frame
        dd <- data.frame(ttt,
                   do.call(rbind.fill,lapply(L,as.data.frame)),
                   row.names=seq(length(ttt)))
        ## rbind.fill fills mismatches with NA
        ## ?? do I want to fill them with 0?  Does it matter?
        dd[is.na(dd)] <- 0
        class(dd) <- c("commframe","data.frame")
        return(dd)
    } else return(list(matList=L,ttt=ttt))
}

## simComm(n.site=c(5,10),spcat=c(2,1))
## simComm(n.site=c(5,10))

##' simulate one community, then
##' rarefy multiple times and return the vector (or matrix) of F statistics
##'
##' @param bstat beta diversity statistic (see \code{\link{calcbeta_stat}})
##' @param rarefy (logical) perform rarefaction?
##' @param .progress passed to \code{simfun1}
##' @param nsim number of simulations
##' @param \dots arguments passed to \code{simfun1}
##' @export
simfun2 <- function(..., bstat, rarefy=TRUE,
                    .progress="text",nsim=100) {
   L <- simComm(...,retval="list")
   nList <- lapply(L$matList,rowSums)  ## numbers per patch
   tList <- findTargets(nList)
   simfun1(L$matList,tList,L$ttt, bstat, rarefy, .progress=.progress,
                    nsim=nsim)
}


##' Rarefaction testing with rarefactions within permutations
##' @importFrom plyr raply
##' @importFrom permute shuffle
##' @param comm community matrix (patches=rows, species=columns)
##' @param ttt treatment vector (treatment id of each patch)
## @param .progress progress bar type -- passed through to \code{\link{raply}}
##' @param n_raref number of rarefactions
##' @param nperm number of permutations
##' @param method beta diversity metric (passed to vegan::betadisper)
##' @param binary reduce data to presence/absence? (passed to vegan::betadisper)
##' @param ptype p-value calculations: "raref_perm_meanF", calculate mean F-statistics within permutations; "raref_perm_meanp", calculate mean p-value across permutations; "raref_perm_medianp" (default), calculate median p-value across permutations.
##' @param return.all (logical) return all information rather than summarized test statistics?
##' @param dummy_raref (logical) [TESTING ONLY] fake rarefaction? (i.e., don't actually do any rarefaction)
##' @param seed random-number seed
##' @return A list with class \code{"htest"} containing the following components:
##' \item{method}{a string giving the method and number of rarefactions used}
##' \item{estimate}{the mean beta diversity (measured as median distance from centroid) across all sites and rarefactions}
##' \item{conf.int}{the confidence interval of the mean beta diversities across rarefactions}
##' \item{statistic}{the median F-statistic across permutations}
##' \item{p.value}{the p-value corresponding to the statistic}
##' @examples
##' set.seed(101)
##' sx <- simComm(n.site=20,p.mix=0.5,n.indiv.site=c(10,20),totsp=c(20,40))
##' plot(sx)
##' raref_test(sx[,-1],sx[,1])
##' tmpf <- function(...) {
##'    sx <- simComm(...)
##'    tt <- try(raref_test(sx[,-1],sx[,1])$p.value)
##'    if (is(tt,"try-error")) return(NA) else return(tt)
##' }
##' \dontrun{
##' ## compute power for this case
##' set.seed(101)
##' system.time(pvec <- replicate(100,tmpf(n.site=20,p.mix=0.5,
##'                        n.indiv.site=c(10,20),totsp=c(20,40))))
##' ## ~ 10 minutes
##' hist(pvec,breaks=20,xlim=c(0,1),col="gray")
##' mean(pvec<0.05)  ## 0; conservative (but very small sample!)
##' }
##' @export
raref_test <- function(comm,ttt,n_raref=50,nperm=200,
                       method="jaccard",binary=TRUE,
                       ptype=c( "raref_perm_medianp",
                               "raref_perm_meanF","raref_perm_meanp"),
                       return.all=FALSE,
                       dummy_raref=FALSE,
                       seed=NULL) {
    ## FIXME: add progress bar?
    if (!is.null(seed)) set.seed(seed)
    ptype <- match.arg(ptype)
    matList <- split.data.frame(comm,ttt)
    nList <- lapply(matList,rowSums)  ## numbers per patch
    ff <- findTargets(nList)
    ## need to match targets list back with community matrix
    ## *in the correct order*!  simple unlist() won't work
    tvec <- unsplit(ff,ttt)
    if (any(tvec==0)) stop("uh-oh")
    dmat <- raply(n_raref,
              {
                  rr <- if (dummy_raref) comm else sampleMat(comm,tvec)
                  betadisper(vegdist(rr,method=method,binary=binary),
                             ttt,bias.adjust=TRUE)$distances
              })
    if (n_raref==1) dmat <- t(matrix(dmat)) ## hack
    ## now set up permutations ... (from permutest.betadisper)
    b0 <- betadisper(vegdist(comm,method=method,binary=binary,
                             bias.adjust=TRUE),ttt)
    mod.aov <- anova(b0)
    nobs <- length(b0$distances)
    res <- matrix(nrow=n_raref,ncol=nperm+1)
    ## calculate the order of all the permutations
    permMat <- t(raply(nperm,shuffle(nobs)))
    for (i in seq(n_raref)) {
        ## for each rarefaction, compute the F statistic for each permutation
        mod <- lm(dmat[i,] ~ ttt)
        mod.Q <- mod$qr
        p <- mod.Q$rank
        rdf <- nobs - p
        resids <- qr.resid(mod.Q, dmat[i,])
        res[i,1] <- summary(mod)$fstatistic[1]
        for (j in seq(nperm)) {
            perm <- permMat[,j]
            perm.resid <- resids[perm]
            f <- qr.fitted(mod.Q, perm.resid)
            mss <- sum((f - mean(f))^2)
            r <- qr.resid(mod.Q, perm.resid)
            rss <- sum(r^2)
            resvar <- rss / rdf
            res[i,j+1] <- (mss / (p - 1)) / resvar
        }
    }
    if (return.all) {
        dimnames(res) <- list(raref=seq(n_raref),
                              perm=c("obs",paste0("p",seq(nperm))))
        return(res)
    }
    ## now calculate mean Fstatistics **within permutations**
    ## still need to sort before taking the mean!
    if (ptype=="raref_perm_meanF")  {
        res[,-1] <- t(apply(res[,-1],1,sort))
        resmeans <- colMeans(res)
        pval <- mean(resmeans[1]<=resmeans[-1])
    }
    else if (ptype=="raref_perm_meanp") {
        pval <- mean(apply(res,1,function(x) mean(x[1]<=x[-1])))
    } else if (ptype=="raref_perm_medianp") {
        pval <- median(apply(res,1,function(x) mean(x[1]<=x[-1])))
    } else stop("unknown ptype")
    obs <- mean(res[,1])
    col_means <- colMeans(dmat)
    grand_mean <- mean(col_means)
    ci <- quantile(col_means,c(0.025,0.975))
    attr(ci,"conf.level")=0.95
    structure(list(
        method=paste0("hierarchical rarefaction with permDISP (n_raref=",
        n_raref,",ptype=",ptype,")"),
        estimate = grand_mean,
        conf.int = ci,
        data.name=deparse(substitute(comm)),
        statistic = c("F"=obs),
        p.value = pval),
              class = "htest")
}


