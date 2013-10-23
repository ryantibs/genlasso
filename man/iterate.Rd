\name{iterate}
\alias{iterate}
\title{
    Iterate genlasso object
}
\description{
    Given an incomplete genlasso path object, this function computes
    additional values of the path. All options are assumed to be the
    same as those in the initial call to a genlasso function, with
    the exception of maxsteps, minlam, and verbose, all of which can
    be changed with a call to iterate.
}
\usage{
iterate(object, maxsteps=2000, minlam=0, verbose=FALSE)
}
\arguments{
  \item{object}{
    A genlasso object with an incomplete path.
  }
  \item{maxsteps}{
    an integer specifying the maximum number of steps for the algorithm
    to take before termination. This is cumulative, so that steps in the
    initial path are counted towards the maxsteps number. Default is 2000. 
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
  Returns an object of class "genlasso", a list with at least following
  components:
  \item{lambda}{
    values of lambda at which the solution path changes slope,
    i.e., kinks or knots.
  }
  \item{beta}{
    a matrix of primal coefficients, each column corresponding to a knot
    in the solution path.
  }
  \item{fit}{
    a matrix of fitted values, each column corresponding to a knot in
    the solution path.
  }
  \item{u}{
    a matrix of dual coefficients, each column corresponding to a knot
    in the solution path.
  }
  \item{hit}{
    a vector of logical values indicating if a new variable in the dual
    solution hit the box contraint boundary. A value of \code{FALSE}
    indicates a variable leaving the boundary. 
  }
  \item{df}{
    a vector giving an unbiased estimate of the degrees of freedom of
    the fit at each knot in the solution path.
  }
  \item{y}{
    the observed response vector. Useful for plotting and other
    methods.
  }
  \item{completepath}{
    a logical variable indicating whether the complete path was
    computed (terminating the path early with the \code{maxsteps} or
    \code{minlam} options results in a value of \code{FALSE}).
  }
  \item{bls}{
    the least squares solution, i.e., the solution at lambda = 0. 
  }
  \item{call}{
    the matched call.
  }
}
\author{
  Ryan J. Tibshirani and Taylor B. Arnold
}
\references{
  Tibshirani, R. J. and Taylor, J. (2011), "The solution path of the
  generalized lasso", Annals of Statistics 39 (3) 1335--1371.
}
\seealso{
  \code{\link{trendfilter}}, \code{\link{fusedlasso}},
  \code{\link{genlasso}}
}
\examples{
   # Wide with an X matrix:
     library(genlasso)
     set.seed(1)
     n = 100
     p = 35
     X = matrix(rnorm(n*p),nrow=n)
     y = 3*X[,1] + rnorm(n)
     D = diag(1,p)

     steps10 = genlasso(y, X=X, D=D, maxsteps=10)
     steps20 = genlasso(y, X=X, D=D, maxsteps=20)
     steps30 = genlasso(y, X=X, D=D, maxsteps=30)
     steps20_iterated = iterate(steps10, maxsteps=20)
     steps30_iterated = iterate(steps20_iterated, maxsteps=30)

     # Check for equality; should match on all but 'call'
     steps20$call = steps20_iterated$call
     steps30$call = steps30_iterated$call
     all.equal(target=steps20, current=steps20_iterated)
     all.equal(target=steps30, current=steps30_iterated)

   # Fused w/ L1 penalty:
     library(genlasso)
     set.seed(1)
     n = 100
     y = as.numeric(row(diag(10)) > 5 & col(diag(10)) > 5) * 3 + rnorm(100)
     D = getD2dSparse(10, 10)

     steps10 = fusedlasso(y, D=D, maxsteps=10, gamma=0.5)
     steps20 = fusedlasso(y, D=D, maxsteps=20, gamma=0.5)
     steps30 = fusedlasso(y, D=D, maxsteps=30, gamma=0.5)
     steps20_iterated = iterate(steps10, maxsteps=20)
     steps30_iterated = iterate(steps20_iterated, maxsteps=30)

     # Check for equality; should match on all but 'call'
     steps20$call = steps20_iterated$call
     steps30$call = steps30_iterated$call
     all.equal(target=steps20, current=steps20_iterated)
     all.equal(target=steps30, current=steps30_iterated)
}
\keyword{models}
