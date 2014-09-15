\name{predict.genlasso}
\alias{predict.genlasso}
\title{
  Make predictions given a genlasso object
}
\description{
  This predict method for the genlasso class makes a prediction for the
  fitted values at new predictor measurements. Hence it is really only
  useful when the generalized lasso model has been fit with a
  nonidentity predictor matrix. In the case that the predictor matrix 
  is the identity, it does the same thing as \code{\link{coef.genlasso}}. 
}
\usage{
\method{predict}{genlasso}(object, lambda, nlam, df, Xnew, ...)
}
\arguments{
  \item{object}{
    object of class "genlasso", or an object which inherits this class 
    (i.e., "fusedlasso", "trendfilter").
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
    at all in the solution path, the least regularized solution at which
    this degrees of freedom value is not exceeded is chosen. 
  }
  \item{Xnew}{
    a numeric matrix X, containing new predictor measurements at which
    predictions should be made. If missing, it is assumed to be
    the same as the existing predictor measurements in \code{object}.
  }
  \item{...}{
    additional arguments passed to predict.
  }
}
\value{
  Returns a list with the following components:
  \item{fit}{
    a numeric matrix of predictor values, one column for each value of
    lambda.
  }
  \item{lambda}{
    a numeric vector containing the sequence of tuning parameter values,
    corresponding to the columns of \code{fit}.
  }
  \item{df}{
    if \code{df} was specified, an integer vector containing the
    sequence of degrees of freedom values corresponding to the columns
    of \code{fit}.
  }
}
\seealso{
  \code{\link{coef.genlasso}}
}
\keyword{methods}
