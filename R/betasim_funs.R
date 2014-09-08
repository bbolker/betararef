
##' Simulate communities with specified diversity characteristics
##' @param n.abund number of abundance categories
##' @param p.abund relative ranking of subsequent abundance categories
##' @param diff.abund alternative parameterization: difference between most and least abundant category
##' @param spcat number of species per abundance category PER SITE (overdetermined)
##' @param totsp total species pool size (FIXME)
##' @param n.site number of sites
##' @param n.indiv.site individuals per site
##' @param n.indiv.tot total number of individuals
##' @param p.mix mixing parameter: probability of among-site mixing for each abundance category: 0=endemic, 1=homogeneous
##' @param seed random number seed
##' @param rand randomization type: 'none' to keep numbers==expected numbers; 'multinom' for fixed total per site; 'poisson' for Poisson sampling, 'trpoisson' for truncated Poisson sampling
##' @param rarefy degree of rarefaction (1=none)
##' @param pred.spec predator specification (stub?)
##' @param pred.rand predator randomness
##' @param pred.gamma predation intensity (-Inf==none)
##' @param pred.alpha predator preference parameter; if \code{pred.spec=="abund"} 0=neutral, -1=pref. rare, +1=pref. common; if pred.spec=="focal" 0=neutral, -1=pref non-endemic, +1=pref endemic
##' @param pred.gamma.sd patch-to-patch variability in predation intensity (0=none)
##' @param no.zero.patch force all patches to have positive occupancy?
##' @importFrom abind abind
##' @importFrom reshape2 dcast
##' @importFrom vegan betadisper raupcrick vegdist
##' @export
betasim <- function(n.abund=5,        
                    p.abund=0.5,      
                    diff.abund=NULL,  
                    spcat=2,          
                    totsp=30,         
                    n.site=3,         
                    n.indiv.site=100,    
                    n.indiv.tot=NULL,    
                    p.mix=rep(0,n.abund),
                    seed=NULL,
                    rand=c("none","multinom","poisson","trpoisson"),
                    rarefy=1,
                    pred.spec=c("none","focal","abund"),
                    pred.rand=c("binom","none"),
                    pred.gamma=-Inf,      
                    pred.alpha=0,         
                    pred.gamma.sd=0,      
                    no.zero.patch=TRUE   
                    ) {
    ## species category/site combinations
    ## was: LETTER.spcat where LETTER was associated with site,
    ## e.g.
    if ((n.mix <- length(p.mix))>1 && n.mix!=n.abund) {
        stop("mismatch between mixing proportions and number of abundance classes")
    }
    if (spcat<1 && !isTRUE(all.equal(round(n.site*spcat),n.site*spcat)))
        stop("if spcat<1, n.site*spcat must be an integer")
    if (!is.null(n.indiv.tot)) {
        if (!missing(n.indiv.site)) stop("must specify at most one of n.indiv.site and n.indiv.tot")
        n.indiv.site <- n.indiv.tot/n.site
        if (floor(n.indiv.site) != n.indiv.site) {
            warning("rounding number of individuals per site")
            n.indiv.site <- round(n.indiv.site)
        }
    }
    p.mix <- rep(p.mix,length.out=n.abund)  ## replicate p.mix (if necessary)
    rand <- match.arg(rand)  ## check/expand 'rand' argument
    pred.spec <- match.arg(pred.spec)
    pred.rand <- match.arg(pred.rand)
    if (!is.null(seed)) set.seed(seed)
    ## Generate proportions in each category (on each reef):
    if (missing(p.abund) && !is.null(diff.abund))
        p.abund <- exp(log(diff.abund)/n)  ## n^th root of diff.abund
    avec <- p.abund^(0:(n.abund-1))
    avec <- avec/sum(avec)
    ## FIXME: adjust warning
    ## if (totsp %% (n.site * n.abund) != 0) 
    ##        warning("totsp is not an even multiple of n.site*n.abund")
    ## FIXME: we could round() if necessary ...
    if (missing(spcat)) {
        spcat <- totsp/(n.site*n.abund)  ## number of endemic species per abundance class per site
    } else if (!missing(n.abund) && !missing(n.site) && !missing(totsp))
        stop("must specify at most three of spcat, n.abund, n.site, totsp")
    ## FIXME: allow filling-in accordingly
    
    ## Set up three-way array:
    ##  dim 1: species class and number within class
    ##         [lower-case letters + lower case roman numerals]: see spcat
    ##         if spcat >= 1, then lower-case letters correspond to the original endemic site
    ##         (a=1, b=2, etc.); otherwise, lower-case letters may be found at a range of sites
    ##  dim 2: abundance class [0 to $n-1$]
    ##  dim 3: site [numbers -- WAS upper case letters]
    ##
    totspecies <- round(n.site*spcat)
    ## n.site <- 20; spcat <- 0.25; totspecies <- 5  ## 5 LETTERS, 1 roman code
    ## n.site <- 5; spcat <- 1; totspecies <- 5   ## 5 LETTERS, 1 roman code
    ## n.site <- 5; spcat <- 2; totspecies <- 10  ## 5 LETTERS, 2 roman codes
    ## ugh.  there should be a better way to get this pattern ...
    betw <- function(x,a,b) {
        y <- x > a & x <= b
        storage.mode(y) <- "numeric"
        y
    }
    ## http://stackoverflow.com/questions/21681785/repeating-vector-of-letters/21682003#21682003
    myletters <- function(n) 
        unlist(Reduce(paste0, 
                      replicate(n %/% length(letters), letters, simplify=FALSE),
                      init=letters,
                      accumulate=TRUE))[1:n]
    xletters <- myletters(500)
    ## generate a matrix with 1:n.site in each row
    m <- matrix(rep(seq(n.site),totspecies),ncol=totspecies)
    ## ?? generating indices for which species category each species is in
    ## at each site ??
    m2 <- t(betw(ceiling(m-(col(m)-1)/spcat),0,ceiling(1/spcat)))
    sp1 <- rep(xletters[seq(totspecies/ceiling(spcat))],each=ceiling(spcat))
    sp2 <- rep(tolower(as.roman(seq(ceiling(spcat)))),length.out=totspecies)
    dimnames(m2) <- list(species=paste(sp1,sp2,sep="."),site=1:n.site)
    ## create array: n.abund copies of 'm2'
    a1 <- abind(replicate(n.abund, m2, simplify=FALSE), along=3)
    a2 <- aperm(a1,c(1,3,2))
    dimnames(a2)[[2]] <- 0:(n.abund-1)
    names(dimnames(a2)) <- c("species","abund","site")
    ## array of which species are endemic in which sites
    sp_array <- d <- a2
    ## Assign proportions to endemic sites:
    d <- sweep(d,2,avec,"*")
    ## Normalize:
    d <- sweep(d,3,apply(d,3,sum),"/")
    if (FALSE) {
        ## examine:
        as.data.frame.table(d["a.i",,"1",drop=FALSE])  ## endemic species
        as.data.frame.table(d["b.i",,"1",drop=FALSE])  ## non-endemic species
    }
    ## mix among sites
    for (j in 1:n.abund) {
        pvec0 <- rep(p.mix[j]/n.site,n.site) ## redistributed individuals
        for (i in 1:totspecies) {
            n.exp <- sum(d[i,j,])
            endem <- (d[i,j,]>0)
            pvec <- pvec0+endem*((1-p.mix[j])/sum(endem))
            d[i,j,] <- n.exp*pvec
        }
    }
    
    noZeros <- FALSE  ## initialize so we always go through once
    ## argh, my brain's not working.
    ## want to run this loop
    ##   once no matter what
    ##   repeat *unless* noZeros=TRUE OR no.zero.patch=FALSE
    ##   stop if noZeros OR !no.zero.patch
    ##   i.e., keep going if !(noZeros && no.zero.patch)
    d2  <- d
    while (!(noZeros || !no.zero.patch)) {
    
        if (rand=="multinom") {  ## fixed number of individuals per site
            for (i in 1:n.site) {
                d2[,,i] <- rmultinom(1,size=n.indiv.site,prob=d[,,i])
            }
        } else if (rand=="poisson") {  ## allow variation
            d2[] <- rpois(length(d),lambda=n.indiv.site*d)
        } else if (rand=="trpoisson") {  ## truncated Poisson
            d2[] <- rtrpois(length(d),lambda=n.indiv.site*d)
        }  else if (rand=="none") {
            d2[] <- n.indiv.site*d
        }

        ## rarefy (should be unnecessary?? can just reduce n, use poisson ...)
        if (rarefy<1) d2[] <- rbinom(length(d),size=d2,prob=rarefy)

        ## predator effects!
        if (pred.gamma>(-Inf)) {
            if (any(floor(d2)!=d2)) {
                warning("rounding to do predator selection")
                d2 <- round(d2)
        }
        ## rules for predator impact:
        ## if non-specific:
        ##   if pred.gamma.sd==0:
        ##      equivalent to rarefying with prob=plogis(pred.gamma)
        ##   if pred.gamma.sd>0:
        ##      rarefy with prob(site)=plogis(rnorm(nsite,pred.gamma,pred.gamma.sd))
        ## if focal:
        ##      as above, but add +pred.alpha to endemic species, -pred.alpha to other species
        ## if abund:
        ##      as above, but add +pred.alpha*(scaled p.abund)
        pred.logis <- d2             ## copy shape of species table
        pred.logis[] <- pred.gamma  ## replace values with baseline predation value
        if (pred.gamma.sd>0) {
            sitepred <- rnorm(n.site,sd=pred.gamma.sd)
            ## dim 3 == site
            pred.logis <- pred.logis + sitepred[slice.index(pred.logis,3)]
        }
        if (pred.spec=="focal") {
            stop("fixme for new organization: dim 1 != origin?")
            ## dim 3 == site, dim 1 == origin
            origsite <- (slice.index(pred.logis,1)==slice.index(pred.logis,3))
            pred.logis[] <- pred.logis[] + ifelse(origsite,pred.alpha,-pred.alpha)
        } else if (pred.spec=="abund") {
            ## dim 3 == abund
            stop("fixme: assign predation on basis of actual local abundance")
            ## FIXME: predation is assigned on the basis of abundance category
            ##   not whether actually abundant at that site
            abundscore <- seq(1,-1,length=n.abund)[slice.index(pred.logis,2)]
            pred.logis[] <- pred.logis[] + pred.alpha*abundscore
        }
        ## survival probability is COMPLEMENTARY plogis() ...
        survprob <- plogis(pred.logis,lower.tail=FALSE)
        d2[] <- if (pred.rand=="binom") {
            rbinom(length(d),size=d2,prob=survprob)
        } else d2*survprob
        }
        
        noZeros <- all(apply(d2,3,sum)>0)
        ## if (!noZeros) browser()
    }
    
    ## rearrange to data frame
    d3 <- as.data.frame.table(d2)

    ## rearrange data frame to species matrix
    species <- n <- abund <- NULL ## shut off false-positive warning about global variables
    d4 <- transform(d3,
                    sp=paste(species,abund,sep="")) ## species names
    ## reorder properly (i.e. in order of occurrence)
    d4$sp <- factor(d4$sp,levels=unique(d4$sp))
    if (any(with(d4,table(site,sp))>1)) stop("uh-oh: repeated species/site combinations")
    d5 <- dcast(d4,site~sp,value.var="Freq")[,-1] ## recast & drop site column
    dnames <- list(dimnames(sp_array)[["site"]],colnames(d5))
    d6 <- as.matrix(d5)
    dimnames(d6) <- dnames  ## restore row/column names
    d7 <- d6[,colSums(d6)>0]
    class(d7) <- c("spmat","matrix")
    d7
}

##' compute beta diversity summary statistic
##'
##' @param m community matrix
##' @param method diversity measure
##' @param distances calculate average distance to centroid or average pairwise distance?
##' @param trap.errors replace error returns with NA values?
##' @param transFun transform (e.g. for taking square-root of abundances)
##' @param binaryFlag compute presence.absence statistics?
##' @param \dots additional parameters for \code{raupcrick}
##' @export
calcbeta <- function(m,method="jaccard",
                     distances=c("centroid","pairwise"),
                     trap.errors=TRUE,
                     transFun=identity,
                     binaryFlag=FALSE,
                     ...) {
    ## calculate beta distances
    if (method=="raupcrick2") {
        return(c(unclass(raupcrick(m, ...))))
    }
    distances <- match.arg(distances)
    nsites <- nrow(m)
    if (missing(binaryFlag) && method %in% c("jaccard","raup"))
        binaryFlag <- TRUE
    if (missing(transFun)) {
        transFun <- switch(method, manhattan=sqrt, gower=, altGower=log,
                           identity)
    }
    vv <- vegdist(transFun(m),method=method,binary=binaryFlag)
    nvals <- switch(distances,centroid=nsites,pairwise=nsites*(nsites-1)/2)
    if (all(vv==0)) return(rep(0,nvals)) ## all identical sites
    ##  if (any(vv==0)) NA else
    if (distances=="centroid") {
        if (!trap.errors) {
            retval <- betadisper(vv,rep("blank",nsites))$distances
        } else {
            tt <- try(betadisper(vv,rep("blank",nsites)))
            if (inherits(tt,"try-error")) {
                retval <- rep(NA,nsites)
            } else retval <- tt$distances
        } 
    } else {
        retval <- c(unclass(vv)) 
    }
    retval
}

## https://stat.ethz.ch/pipermail/r-help/2005-May/070678.html
## m1 = m2/(1-exp(-m2))  ## pre-truncation mean to truncated mean
## m1-m1*exp(-m2)=m2
## m1 = m2*(1+m1*exp(-m2))
## according to Wolfram Alpha,
## m2 = W(-m1*exp(-m1))+m1
##
## m2 <- 5
## m1 <- m2/(1-exp(-m2))
## emdbook::lambertW(-m1*exp(-m1))+m1

rtrpois <- function(n,lambda) {
    ## recover pre-truncation mean from truncated mean
    lambda0 <- lambertW(-lambda*exp(-lambda))+lambda
    U <- runif(n)   # the uniform sample
    t0 = -log(1 - U*(1 - exp(-lambda0))) # the "first" event-times
    T1<-(lambda0 - t0)   # the set of (T-t)
    rpois(n,T1)+1 # the final truncated Poisson sample
}
## test
## mean(rtrpois(100000,lambda=2.5))
## mean(rtrpois(100000,lambda=17))
