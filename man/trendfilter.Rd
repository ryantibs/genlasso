\name{trendfilter}
\alias{trendfilter}
\title{
  Compute the trend filtering solution path for any polynomial order
}
\description{
  This function computes the solution path for the trend filtering
  problem of an arbitrary polynomial order. When the order is set to
  zero, trend filtering is equivalent to the 1d fused lasso, see
  \code{\link{fusedlasso1d}}.
}
\usage{
trendfilter(y, pos, X, ord = 1, approx = FALSE, maxsteps = 2000,
            minlam = 0, rtol = 1e-07, btol = 1e-07, eps = 1e-04,
            verbose = FALSE)
}
\arguments{
  \item{y}{
    a numeric response vector.
  }
  \item{pos}{
    an optional numeric vector specifying the positions of the
    observations, and missing \code{pos} is assumed to mean unit spacing.
  }
  \item{X}{
    an optional matrix of predictor variables, with observations along
    the rows, and variables along the columns. If the passed \code{X} 
    has more columns than rows, then a warning is given, and a small ridge
    penalty is added to the generalized lasso criterion before the path
    is computed. If \code{X} has less columns than rows, then its rank is
    not checked for efficiency, and (unlike the \code{genasso} function) a 
    ridge penalty is not automatically added if it is rank deficient. 
    Therefore, a tall, rank deficient \code{X} may cause errors. 
  }
  \item{ord}{
    an integer specifying the desired order of the piecewise polyomial
    produced by the solution of the trend filtering problem. Must be
    non-negative, and the default to 1 (linear trend filtering).
  }
  \item{approx}{
    a logical variable indicating if the approximate solution path
    should be used (with no dual coordinates leaving the boundary).
    Default is \code{FALSE}. Note 
    that for the 1d fused lasso (zeroth order trend filtering), with
    identity predictor matrix, this approximate path is the same as 
    the exact solution path. 
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
    in the case that \code{X} is wide (more columns than rows). If numeric
    problems occur, make \code{eps} larger. Default is 1e-4.
  }
  \item{verbose}{
    a logical variable indicating if progress should be reported after
    each knot in the path.
  }
}
\details{
  When the predictor matrix is the identity, trend filtering fits a
  piecewise polynomial to linearly ordered observations. The result is 
  similar to that of a polynomial regression spline or a smoothing
  spline, except the knots in the piecewise polynomial (changes in the 
  (k+1)st derivative, if the polynomial order is k) are chosen
  adaptively based on the observations. This is in contrast to
  regression splines, where the knots are prespecified, and smoothing
  splines, which place a knot at every data point. 

  With a nonidentity predictor matrix, the trend filtering problem
  enforces piecewise polynomial smoothness along successive components
  of the coefficient vector. This can be used to fit a kind of varying
  coefficient model.

  We note that, in the signal approximator (identity predictor matrix) 
  case, fitting trend filtering estimate with arbitrary positions \code{pos} 
  is theoretically no harder than doing so on an evenly spaced grid. However
  in practice, with differing gaps between points, the algorithm can 
  become numerically unstable even for large (or moderately large) problems. 
  This is especially true as the polynomial order increases. Hence, use the
  positions argument \code{pos} with caution. 
}
\value{
  Returns an object of class "trendfilter", a subclass of
  "genlasso". This is a list with at least following components:
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
    the least squares solution, i.e., the solution at lambda = 0. This
    can be \code{NULL} when \code{completepath} is \code{FALSE}. 
  }
  \item{ord}{
    the order of the piecewise polyomial that has been fit.
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

  Tibshirani, R. J. (2014), "Adaptive piecewise polynomial estimation
  via trend filtering", Annals of Statistics 42 (1): 285--323.

  Arnold, T. B. and Tibshirani, R. J. (2014), "Efficient implementations
  of the generalized lasso dual path algorithm", arXiv: 1405.3222.

  Kim, S.-J., Koh, K., Boyd, S. and Gorinevsky, D. (2009), "l1 trend
  filtering", SIAM Review 51 (2), 339--360.
}
\seealso{
  \code{\link{fusedlasso1d}}, \code{\link{genlasso}},
  \code{\link{cv.trendfilter}}, \code{\link{plot.trendfilter}}
}
\examples{
# Constant trend filtering (the 1d fused lasso)
set.seed(0)
n = 100
beta0 = rep(sample(1:10,5),each=n/5)
y = beta0 + rnorm(n,sd=0.8)
a = fusedlasso1d(y)
plot(a)

# Linear trend filtering
set.seed(0)
n = 100
beta0 = numeric(n)
beta0[1:20] = (0:19)*4/19+2
beta0[20:45] = (25:0)*3/25+3
beta0[45:80] = (0:35)*9/35+3
beta0[80:100] = (20:0)*4/20+8
y = beta0 + rnorm(n)
a = trendfilter(y,ord=1)
plot(a,df=c(2,3,4,10))

# Cubic trend filtering
set.seed(0)
n = 100
beta0 = numeric(100)
beta0[1:40] = (1:40-20)^3
beta0[40:50] = -60*(40:50-50)^2 + 60*100+20^3
beta0[50:70] = -20*(50:70-50)^2 + 60*100+20^3
beta0[70:100] = -1/6*(70:100-110)^3 + -1/6*40^3 + 6000
beta0 = -beta0
beta0 = (beta0-min(beta0))*10/diff(range(beta0))
y = beta0 + rnorm(n)
a = trendfilter(y,ord=3)
plot(a,nlam=5)
}
\keyword{models}
