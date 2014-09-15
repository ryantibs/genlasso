\name{coef.genlasso}
\alias{coef.genlasso}
\title{
  Extract cefficients from a genlasso object
}
\description{
  This function extracts coefficients from a generalized lasso solution
  path object, for any set of tuning parameter values along the path. It
  can return dual coefficients. The requested coefficients can also be
  parametrized by degrees of freedom value instead of tuning parameter
  value.
}
\usage{
\method{coef}{genlasso}(object, lambda, nlam, df, 
     type = c("primal", "dual", "both"), ...)
}
\arguments{
  \item{object}{
    an object of class "genlasso", or an object which inherits this
    class (i.e., "fusedlasso", "trendfilter"). 
  }
  \item{lambda}{
    a numeric vector of tuning parameter values at which coefficients
    should be calculated. The user can choose to specify one of
    \code{lambda}, \code{nlam}, or \code{df}; if none are specified,
    then coefficients are returned at every knot in the solution path. 
  }
  \item{nlam}{
    an integer indicating a number of tuning parameters values at which
    coefficients should be calculated. The tuning parameter values are
    then chosen to be equally spaced on the log scale over the first
    half of the solution path (this is if the full solution path has
    been computed; if only a partial path has been computed, the tuning 
    parameter values are spaced over the entirety of the computed path).  
  }
  \item{df}{
    an integer vector of degrees of freedom values at which coefficients
    should be calculated. In the case that a single degrees of freedom
    value appears multiple times throughout the solution path, the least 
    regularized solution (corresponding to the smallest value
    of lambda) is chosen. If a degrees of freedom value does not appear
    at all in the solution path, this function chooses the solution whose 
    degrees of freedom is largest, subject to being less than or equal to 
    the specified value.
  }
  \item{type}{
    a character string, one of "primal", "dual", or "both", indicating
    whether primal coefficients, dual coefficients, or both, should be
    returned. Default is "primal", which corresponds to the solution of
    the original problem.
  }
  \item{...}{
    additional arguments passed to coef.
  }
}
\value{
  Returns a list with the following components:
  \item{beta}{
    if the type is "primal" or "both", a matrix containing the primal
    coefficients, each column corresponding to a value of lambda.
  }
  \item{u}{
    if the type is "dual" or "both", a matrix containing the dual
    coefficients, each column corresponding to a value of lambda.
  }
  \item{lambda}{
    a numeric vector containing the sequence of tuning parameter values,
    corresponding to the columns of \code{beta} and \code{u}.
  }
  \item{df}{
    an integer vector containing the sequence of degrees of freedom 
    values corresponding to the columns of \code{beta} and \code{u}. 
  }
}
\seealso{
  \code{\link{genlasso}}, \code{\link{predict.genlasso}},
  \code{\link{plot.genlasso}}
}
\examples{
# Constant trend filtering (the 1d fused lasso)
set.seed(0)
n = 20
beta0 = rep(sample(1:10,5),each=n/5)
y = beta0 + rnorm(n,sd=0.5)
a = fusedlasso1d(y)

# Get the coefficients that use 3, 4, and 5 degrees
# of freedom 
coef(a,df=3:5)
}
\keyword{methods}
