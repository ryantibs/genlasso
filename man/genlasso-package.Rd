\name{genlasso-package}
\alias{genlasso-package}
\docType{package}
\title{
  Package to compute the solution path of generalized lasso problems
}
\description{
  This package is centered around computing the solution path of the
  generalized lasso problem, which minimizes the criterion
  \deqn{
    1/2 \|y - X \beta\|_2^2 + \lambda \|D \beta\|_1.
  }
  The solution path is computed by solving the equivalent Lagrange
  dual problem. The dimension of the dual variable u is the number of
  rows of the penalty matrix D, and the primal (original) and dual
  solutions are related by
  \deqn{
    \hat{\beta} = y - D^T \hat{u}
  }
  when the predictor matrix X is the identity, and
  \deqn{
    \hat{\beta} = (X^T X)^{-1} (X^T y - D^T \hat{u})
  }
  for a full column rank predictor matrix X. For column rank deficient
  matrices X, the solution path is not unique and not computed by this
  package. However, one can add a small ridge penalty to the above
  criterion, which can be re-expressed as a generalized lasso problem
  with full column rank predictor matrix X and hence yields a unique
  solution path.

  Important use cases include the fused lasso, where D is the oriented
  incidence matrix of some underlying graph (the orientations being
  arbitrary), and trend filtering, where D is the discrete difference
  operator of any given order k.

  The general function \code{\link{genlasso}} computes a solution path
  for any penalty matrix D and full column rank predictor matrix X
  (adding a ridge penalty when X is rank deficient). For the fused lasso
  and trend filtering problems, the specialty functions
  \code{\link{fusedlasso}} and \code{\link{trendfilter}} should be used 
  as they deliver a significant increase in speed and numerical
  stability.

  For a walk-through of using the package for statistical modelling see
  the included package vignette; for the appropriate background material
  see the generalized lasso paper referenced below. 
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
  \code{\link{genlasso}}, \code{\link{fusedlasso}},
  \code{\link{trendfilter}} 
}
\keyword{package}

