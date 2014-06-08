##' Plot species matrices and community frames
##'
##' @param x a species matrix (class \code{spmat}) or community data frame (class \code{commframe})
##' @param y ignored, for generic compatibility
##' @param \dots ignored, for generic compatibility
##' @param aspect aspect ratio (\code{"fill"} or \code{"iso"})
##' @importFrom Matrix Matrix
##' @importClassesFrom Matrix Matrix
##' @importMethodsFrom Matrix image
##' @export plot.spmat
plot.spmat <- function(x,y, aspect="fill", ...) {
    class(x) <- "matrix" ## strip spmat
    image(Matrix(x),xlab="Species",ylab="Patch",sub="",aspect=aspect,
          scales=list(y=list(at=1:nrow(x)),
                      x=list(at=1:ncol(x),labels=colnames(x))))
}

##' @importFrom lattice levelplot
##' @importFrom reshape2 melt
##' @export plot.commframe
##' @rdname plot.spmat
plot.commframe <- function(x,y,aspect="fill", ...) {
  m <- x[,-1]
  ttt <- x[,1]
  dd <-   data.frame(ttt=rep(ttt,ncol(m)),
                    melt(as.matrix(m)))
  gg <- grey(seq(from = 0.7, to = 0, length = 100))
  levelplot(value~Var2*Var1|ttt,data=dd,col.regions=gg,
      xlab="Species",ylab="Patch",sub="",aspect=aspect,
           scales=list(x=list(at=1:ncol(m),labels=colnames(m)),
           relation="free"))
}

