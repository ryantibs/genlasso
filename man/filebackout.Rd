\name{filebackout}
\alias{filebackout}
\title{
  Process file backed fused lasso output
}
\description{
  This function processes output from a call to \code{\link{fusedlasso}}
  or \code{\link{fusedlasso}} when file backing is turned on.
}
\usage{
filebackout(object, y, X, file)
}
\arguments{
  \item{object}{
    a list object returned from either \code{\link{fusedlasso}} or
    \code{\link{fusedlasso2d}} when file backing is turned on. The user
    can specify this object, or \code{y}, \code{X}, and the file name
    \code{file}. 
  }
  \item{y}{
    a numeric response vector.
  }
  \item{X}{
    a numeric matrix of predictors, which if missing, is interpreted as
    the identity matrix.
  }
  \item{file}{
    a character string indicating the name of the file.
  }
}
\details{
  For very large problems, it may be easier to extract the output
  manually via \code{readline}. Another benefit of the file backing is
  that the model, up to the current point, can be extracted even if the
  model fitting process quit for some reason. 
}
\value{
  Returns an object of class "fusedlasso" and subclass "genlasso", which
  should be indentical in all regards to the result of a non-file backed
  output, with the one minor exception of the "call" argument. This is a
  list with at least the following components:
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
    the value of the lambda ratio parameter.
  }
  \item{call}{
    the matched call.
  }
}
\seealso{
  \code{\link{fusedlasso}}, \code{\link{fusedlasso2d}}
}
\keyword{utilities}
