\name{genlasso}
\alias{genlasso}
\alias{print.genlasso}
\alias{print.summary.genlasso}
\alias{summary.genlasso}
\title{
  Compute the generalized lasso solution path for arbitrary penalty matrix
}
\description{
  This function computes the solution path of the generalized lasso
  problem for an arbitrary penalty matrix. Speciality functions exist
  for the trend filtering and fused lasso problems; see
  \code{\link{trendfilter}} and \code{\link{fusedlasso}}.
}
\usage{
genlasso(y, X, D, approx = FALSE, maxsteps = 2000, minlam = 0,
         rtol = 1e-07, btol = 1e-07, eps = 1e-4, verbose = FALSE,
         svd = FALSE)
}
\arguments{
  \item{y}{
    a numeric response vector.
  }
  \item{X}{
    an optional matrix of predictor variables, with observations along
    the rows, and variables along the columns. If missing, \code{X} is
    assumed to be the identity matrix. If the passed \code{X} does not
    have full column rank, then a warning is given, and a small ridge
    penalty is added to the generalized lasso criterion before the path
    is computed.
  }
  \item{D}{
    a penalty matrix. Its number of columns must be equal to the number
    of columns of \code{X}, or if no \code{X} is given, the length of 
    \code{y}. This can be a sparse matrix from \code{Matrix} package,
    but this will be ignored (converted to a dense matrix) if \code{D}
    is row rank deficient or if \code{X} is specified. See "Details"
    below. 
  }
  \item{approx}{
    a logical variable indicating if the approximate solution path
    should be used (with no dual coordinates leaving the boundary).
    Default is \code{FALSE}. 
  }
  \item{maxsteps}{
    an integer specifying the maximum number of steps for the algorithm
    to take before termination. Default is 2000. 
  }
  \item{minlam}{
    a numeric variable indicating the value of lambda at which the path
    should terminate. Default is 0.
  }
  \item{rtol}{
    a numeric variable giving the tolerance for determining the rank of
    a matrix: if a diagonal value in the R factor of a QR decomposition
    is less than R, in absolute value, then it is considered zero. Hence
    making rtol larger means being less stringent with determination of
    matrix rank. In general, do not change this unless you know what you
    are getting into! Default is 1e-7.
  }
  \item{btol}{
    a numeric variable giving the tolerance for accepting "late" hitting
    and leaving times: future hitting times and leaving times should always 
    be less than the current knot in the path, but sometimes for numerical
    reasons they are larger; any computed hitting or leaving time larger 
    than the current knot + btol is thrown away. Hence making btol larger
    means being less stringent withthe determination of hitting and leaving 
    times. Again, in general, do not change this unless you know what you 
    are getting into! Default is 1e-7.
  }
  \item{eps}{
    a numeric variable indicating the multiplier for the ridge penalty,
    in the case that \code{X} is column rank deficient. Default is
    1e-4. 
  }
  \item{verbose}{
    a logical variable indicating if progress should be reported after
    each knot in the path. 
  }
  \item{svd}{
    a logical variable indicating if the genlasso function should use 
    singular value decomposition to solve least squares problems at each
    path step, which is slower, but should be more stable.
  }
}
\details{
  The generalized lasso estimate minimizes the criterion
  \deqn{
    1/2 \|y - X \beta\|_2^2 + \lambda \|D \beta\|_1.
  }
  The solution \eqn{\hat{\beta}} is computed as a function of
  the regularization parameter \eqn{\lambda}. The advantage of the
  \code{genlasso} function lies in its flexibility, i.e., the user can
  specify any penalty matrix \code{D} of their choosing. However, for a
  trend filtering problem or a fused lasso problem, it is strongly
  recommended to use one of the speciality functions,
  \code{\link{trendfilter}} or \code{\link{fusedlasso}}. When compared
  to these functions, \code{genlasso} is not as numerically stable and
  much less efficient.

  Note that, when \code{D} is passed as a sparse matrix, the linear 
  systems that arise at each step of the path algorithm are solved
  separately via a sparse solver. The usual strategy (when \code{D} is 
  simply a matrix) is to maintain a matrix factorization of \code{D},
  and solve these systems by (or downdating) this factorization, as
  these linear systems are highly related. Therefore,
  when \code{D} is sufficiently sparse and structured, it can be
  advantageous to pass it as a sparse matrix; but if \code{D} is truly
  dense, passing it as a sparse matrix will be highly inefficient. 
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
  Taylor B. Arnold and Ryan J. Tibshirani
}
\references{
  Tibshirani, R. J. and Taylor, J. (2011), "The solution path of the
  generalized lasso", Annals of Statistics 39 (3) 1335--1371.

  Arnold, T. B. and Tibshirani, R. J. (2014), "Efficient implementations
  of the generalized lasso dual path algorithm", arXiv: 1405.3222.
}
\seealso{
  \code{\link{trendfilter}}, \code{\link{fusedlasso}},
  \code{\link{coef.genlasso}}, \code{\link{predict.genlasso}},
  \code{\link{plot.genlasso}} 
}
\examples{
# Using the generalized lasso to run a standard lasso regression
# (for example purposes only! for pure lasso problems, use LARS
# instead)
set.seed(1)
n = 100
p = 10
X = matrix(rnorm(n*p),nrow=n)
y = 3*X[,1] + rnorm(n)
D = diag(1,p)
out = genlasso(y,X,D)
coef(out, lambda=sqrt(n*log(p)))
}
\keyword{models}
