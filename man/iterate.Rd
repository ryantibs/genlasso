\name{iterate}
\alias{iterate}
\title{
    Iterate a genlasso object
}
\description{
    Given an incomplete \code{genlasso} path object, this function continues 
    the path computation from the last computed knot, either until 
    the complete path has been computed or the step limit specified by
    \code{moresteps} has been reached. All options are assumed
    to be the same as those in the initial call to a genlasso function (as in
    \code{\link{genlasso}}, \code{\link{fusedlasso}}, or \code{\link{trendfilter}}), 
    with the exception of \code{minlam} and \code{verbose}, which can be changed 
    with a call to iterate.
}
\usage{
iterate(object, moresteps=200, minlam=0, verbose=FALSE)
}
\arguments{
  \item{object}{
    a genlasso object with an incomplete path.
  }
  \item{moresteps}{
    an integer specifying the number of additional steps to take, starting
    from termination point of the passed (incomplete) path object.
  }
  \item{minlam}{
    a numeric variable indicating the value of lambda at which the path
    should terminate. Default is 0.
  }
  \item{verbose}{
    a logical variable indicating if progress should be reported after
    each knot in the path.
  }
}
\value{
  Returns an list of the same class typing and same structure as the passed
  \code{object}.
}
\seealso{
  \link{genlasso},
  \link{trendfilter}, \link{fusedlasso}
}
\examples{
  # Sparse 2d fused lasso
  library(genlasso)
  set.seed(1)
  dim1 = dim2 = 10
  n = 100
  y = as.numeric(row(diag(dim1)) > 5 & col(diag(dim2)) > 5) * 3 + rnorm(n)

  a10 = fusedlasso2d(y, dim1=dim1, dim2=dim2, gamma=0.5, maxsteps=10)
  a20 = fusedlasso2d(y, dim1=dim1, dim2=dim2, gamma=0.5, maxsteps=20)
  a30 = fusedlasso2d(y, dim1=dim1, dim2=dim2, gamma=0.5, maxsteps=30)
  b20 = iterate(a10, moresteps=10)
  b30 = iterate(b20, moresteps=10)

  # Check for equality; should match on all but 'call'
  b20$call = a20$call
  b30$call = a30$call
  all.equal(target=a20, current=b20)
  all.equal(target=a30, current=b30)
}
\keyword{models}
