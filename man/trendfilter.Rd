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
trendfilter(y, X, ord = 1, approx = FALSE, maxsteps = 2000,
            minlam = 0, tol = 1e-11, verbose = FALSE)
}
\arguments{
  \item{y}{
    a numeric response vector.
  }
  \item{X}{
    an optional matrix of predictor variables, with observations along
    the rows, and variables along the columns. It is assumed to have
    full column rank (linear independent predictors), but this is not
    checked, for efficiency. A rank deficient \code{X} may or may not
    cause errors, depending on the type of rank deficiency. If missing,
    \code{X} is assumed to be the identity matrix. 
  }
%  \item{z}{
%    an optional numeric vector specifying the positions of the
%    observations, and missing \code{z} is assumed to mean unit spacing.
%    This can only be used when the predictor matrix \code{X} is missing
%    (i.e., identity predictor matrix).
%  }
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
    identity predictor matrix,
%    and evenly spaced grid positions, 
    this approximate path is the same as the exact solution path. 
  }
  \item{maxsteps}{
    an integer specifying the maximum number of steps for the algorithm
    to take before termination. Default is 2000. 
  }
  \item{minlam}{
    a numeric variable indicating the value of lambda at which the path
    should terminate. Default is 0.
  }
  \item{tol}{
    a numeric variable giving the tolerance used in the calculation of
    the hitting and leaving times. A larger value is more conservative,
    and may cause the algorithm to miss some hitting or leaving events
    (do not change unless you know what you're getting into!). Default
    is 1e-11. 
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
  regression splines, where the knots are pre-specified, and smoothing
  splines, which place a knot at every data point.

  With a non-identity predictor matrix, the trend filtering problem
  enforces piecewise polynomial smoothness along successive components
  of the coefficient vector. This can be used to fit a kind of varying
  coefficient model. 
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
  \item{trendorder}{
    the order of the piecewise polyomial that has been fit.
  }
  \item{call}{
    the matched call.
  }
}
\author{
  Ryan J. Tibshirani and Taylor B. Arnold
}
\references{
  Tibshirani, R. J. (2013), "Adaptive piecewise polynomial estimation
  via trend filtering", arXiv: 1304.2986.
	
  Tibshirani, R. J. and Taylor, J. (2011), "The solution path of the
  generalized lasso", Annals of Statistics 39 (3) 1335--1371.

  Kim, S.-J., Koh, K., Boyd, S. and Gorinevsky, D. (2009), "l1 trend
  filtering", SIAM Review 51 (2), 339--360.
}
\seealso{
  \code{\link{fusedlasso1d}}, \code{\link{genlasso}},
  \code{\link{cv.trendfilter}}, \code{\link{plot.trendfilter}}
}
\section{Notes}{
  Trend filtering inherently assumes that the input points are evenly
  spaced. Properly extending trend filtering to the case of unevenly 
  spaced inputs is a work in progress.
%  Fitting trend filtering estimate with arbitrary positions \code{z} is
%  theoretically no harder than doing so on an evenly spaced grid, but in
%  practice it can become numerically unstable even for moderately sized
%  problems, especially as the polynomial order increases. Use the
%  positions argument \code{z} with caution. 
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
