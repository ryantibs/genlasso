\name{softthresh}
\alias{softthresh}
\title{
  Fit a sparse variant of the fused lasso
}
\description{
  This function computes solution path to a fused lasso problem of the
  form
  \deqn{
    1/2 \sum_{i=1}^n (y_i - \beta_i)^2 + \lambda \sum_{(i,j) \in E}
    |\beta_i - \beta_j| + \gamma \cdot \lambda \sum_{i=1}^p |\beta_i|, 
  }
  given the solution path corresponding to \eqn{\gamma=0}. Note that the 
  predictor matrix here is the identity, and in this case the new
  solution path is given by a simple soft-thresholding operation
  (Friedman et al. 2007).
}
\usage{
  softthresh(object, lambda, gamma)
}
\arguments{
  \item{object}{
    an object of class "fusedlasso", fit with no predictor matrix
    \code{X} (taken to mean that the predictor matrix is the identity)
    and with \code{gamma} set to 0. Other objects will issue a warning
    that soft-thresholding does not give the exact primal solution path
    to a sparsified generalized lasso problem.
  }
  \item{lambda}{
    a numeric vector giving the values of lambda at which the solution
    should be computed and returned; if missing, defaults to the knots
    in the solution path stored in \code{object}.
  }
  \item{gamma}{
    a numeric variable giving the ratio of the fusion and sparsity
    tuning parameters, must be greater than or equal to 0. 
  }
}
\value{
  Returns a numeric matrix of primal solutions, one column for each
  value of lambda.
}
\references{
  Friedman J., Hastie T., Hoefling H. and Tibshirani, R. (2007),
  "Pathwise coordinate optimization", Annals of Applied Statistics 
  1 (2) 302--332.
}
\seealso{
  \code{\link{fusedlasso}}
}
\examples{
# The 1d fused lasso
set.seed(0)
n = 100
beta0 = rep(sample(1:10,5),each=n/5)
beta0 = beta0-mean(beta0)
y = beta0 + rnorm(n,sd=0.8)
a = fusedlasso1d(y)

lambda = 4
b1 = coef(a,lambda=lambda)$beta

gamma = 0.5
b2 = softthresh(a,lambda=lambda,gamma=gamma)

plot(1:n,y)
lines(1:n,b1)
lines(1:n,b2,col="red")
legend("topright",lty=1,col=c("black","red"),
       legend=c(expression(gamma==0),expression(gamma==0.5)))
}
\keyword{utilities}
