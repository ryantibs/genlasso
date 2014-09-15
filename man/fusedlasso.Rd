\name{fusedlasso}
\alias{fusedlasso}
\alias{fusedlasso1d}
\alias{fusedlasso2d}
\title{
  Compute the fused lasso solution path for a general graph, or a 1d
  or 2d grid}
\description{
  These functions produce the solution path for a general fused lasso
  problem. The \code{fusedlasso} function takes either a penalty matrix
  or a graph object from the \code{igraph} package. The
  \code{fusedlasso1d} and \code{fusedlasso2d} functions are convenience
  functions that construct the penalty matrix over a 1d or 2d grid.}
\usage{
fusedlasso(y, X, D, graph, gamma = 0, approx = FALSE, maxsteps = 2000,
           minlam = 0, rtol = 1e-07, btol = 1e-07, eps = 1e-4, 
	   verbose = FALSE)
fusedlasso1d(y, pos, X, gamma = 0, approx = FALSE, maxsteps = 2000,
             minlam = 0, rtol = 1e-07, btol = 1e-07, eps = 1e-4, 
	     verbose = FALSE)
fusedlasso2d(y, X, dim1, dim2, gamma = 0, approx = FALSE, maxsteps = 2000,
	     minlam = 0, rtol = 1e-07, btol = 1e-07, eps = 1e-4, 
	     verbose = FALSE)
}
\arguments{
  \item{y}{
    a numeric response vector. Alternatively, for \code{fusedlasso2d}
    with no matrix \code{X} passed, \code{y} can be a matrix
    (its dimensions corresponding to the underlying 2d grid). Note that
    when \code{y} is given as a vector in \code{fusedlasso2d}, with no
    \code{X} passed, it should be in column major order.
  }
  \item{pos}{
    only for \code{fusedlasso1d}, these are the optional positions of
    the positions in the 1d grid. If missing, the 1d grid is assumed to
    have unit spacing.
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
  \item{D}{
    only for \code{fusedlasso}, this is the penalty matrix, i.e., the
    oriented incidence matrix over the underlying graph (the orientation
    of each edge being arbitrary). Only one of \code{D} or \code{graph}
    needs to be specified.
  }
  \item{graph}{
    only for \code{fusedlasso}, this is the underlying graph as an
    \code{igraph} object from the \code{igraph} package. Only one of
    \code{D} or \code{graph} needs to be specified.
  }
  \item{dim1}{
    only for \code{fusedlasso2d}, this is the number of rows in the
    underlying 2d grid. If missing and \code{y} is given as a matrix, it
    is assumed to be the number of rows of \code{y}.
  }
  \item{dim2}{
    only for \code{fusedlasso2d}, this is the number of columns in the
    underlying 2d grid. If missing and \code{y} is given as a matrix, it
    is assumed to be the number of columns of \code{y}.
  }
  \item{gamma}{
    a numeric variable greater than or equal to 0, indicating the ratio
    of the two tuning parameters, one for the fusion penalty, and the
    other for the pure \eqn{\ell_1} penalty. Default is 0. See
    "Details" for more information.
  }
  \item{approx}{
    a logical variable indicating if the approximate solution path
    should be used (with no dual coordinates leaving the boundary).
    Default is \code{FALSE}. Note
    that for the 1d fused lasso, with identity predicor matrix,
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
  The fused lasso estimate minimizes the criterion
  \deqn{
    1/2 \sum_{i=1}^n (y_i - x_i^T \beta_i)^2 + \lambda \sum_{(i,j) \in E}
    |\beta_i - \beta_j| + \gamma \cdot \lambda \sum_{i=1}^p |\beta_i|,
  }
  where \eqn{x_i} is the ith row of the predictor matrix and \eqn{E} is
  the edge set of the underlying graph. The solution \eqn{\hat{\beta}} is
  computed as a function of the regularization parameter \eqn{\lambda},
  for a fixed value of \eqn{\gamma}. The default is to set
  \eqn{\gamma=0}, which corresponds to pure fusion of the coefficient
  vector \eqn{\beta}. A choice \eqn{\gamma>0} introduces both sparsity
  and fusion in the coefficient vector, with a higher value placing more
  priority on sparsity.

  If the predictor matrix is the identity, and the primal solution path
  \eqn{\beta} is desired at several levels of the ratio parameter
  \eqn{\gamma}, it is much more efficient to compute the solution path
  once with \eqn{\gamma=0}, and then use soft-thresholding via the
  \code{\link{softthresh}} function.

  Finally, for the image denoising problem, i.e., the fused lasso over a 2d
  grid with identity predictor matrix, it is easy to specify a huge graph
  with a seemingly small amount of data. For instance, running the 2d
  fused lasso (with identity predictor matrix) on an image at standard
  1080p HD resolution yields a graph with over 2 million
  edges. Moreover, in image denoising problems---somewhat unlike most
  other applications of the fused lasso (and generalized lasso)---a
  solution is often desired near the dense end of the path
  (\eqn{\lambda=0}) as opposed to the regularized end
  (\eqn{\lambda=\infty}). The dual path algorithm implemented by the
  \code{fusedlasso2d} function begins at the fully regularized end
  and works its way down to the dense end. For a problem with many
  edges (dual variables), if a solution at the dense is desired, then it
  must usually pass through a huge number knots in the path. Hence it is
  not advisable to run \code{fusedlasso2d} on image denoising problems of
  large scale, as the dual solution path is computationally
  infeasible. It should be noted that a faster algorithm for the 2d
  fused lasso solution path (when the predictor matrix is the identity),
  which begins at the dense end of the path, is available in the
  \code{flsa} package.
}
\value{
  The function returns an object of class "fusedlasso", and subclass
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
  \item{gamma}{
    the value of the lambda ratio.
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

  Tibshirani, R., Saunders, M., Rosset, S., Zhu, J. and Knight,
  K. (2005), "Sparsity and smoothness via the fused lasso", Journal of
  the Royal Statistics Society: Series B 67(1), 91--108.
}
\seealso{
  \code{\link{softthresh}}, \code{\link{genlasso}}
}
\examples{
# Fused lasso on a custom graph
set.seed(0)
edges = c(1,2,1,3,1,5,2,4,2,5,3,6,3,7,3,8,6,7,6,8)
gr = graph(edges=edges,directed=FALSE)
plot(gr)
y = c(1,1,0,1,1,0,0,0) + rnorm(8,0.1)

# Can either pass the graph object directly, or
# first construct the penalty matrix, and then
# pass this
a1 = fusedlasso(y,graph=gr)
D = getDgSparse(gr)
a2 = fusedlasso(y,D=D)

plot(a1,numbers=TRUE)

\dontrun{
# The 2d fused lasso with a predictor matrix X
set.seed(0)
dim1 = dim2 = 16
p = dim1*dim2
n = 300
X = matrix(rnorm(n*p),nrow=n)
beta0 = matrix(0,dim1,dim2)
beta0[(row(beta0)-dim1/2)^2 + (col(beta0)-dim2/2)^2 <=
(min(dim1,dim2)/3)^2] = 1
y = X \%*\% as.numeric(beta0) + rnorm(n)

# Takes about 30 seconds for the full solution path
out = fusedlasso2d(y,X,dim1=dim1,dim2=dim2)

# Grab the solution at 8 values of lambda over the path
a = coef(out,nlam=8)

# Plot these against the true coefficients
par(mar=c(1,1,2,1),mfrow=c(3,3))
cols = terrain.colors(30)
zlim = range(c(range(beta0),range(a$beta)))
image(beta0,col=cols,zlim=zlim,axes=FALSE)

for (i in 1:8) {
  image(matrix(a$beta[,i],nrow=dim1),col=cols,zlim=zlim,
  axes=FALSE)
  mtext(bquote(lambda==.(sprintf("\%.3f",a$lambda[i]))))
}
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
