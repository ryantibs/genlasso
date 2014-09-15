\name{cv.trendfilter}
\alias{cv.trendfilter}
\title{
  Perform k-fold cross-validation to choose a trend filtering model
}
\description{
  This function performs k-fold cross-validation to choose the value of
  the regularization parameter lambda for a trend filtering problem,
  given the computed solution path. This function only applies to trend
  filtering objects with identity predictor matrix (no \code{X} passed).
}
\usage{
cv.trendfilter(object, k = 5, mode = c("lambda", "df"),
               approx = FALSE, rtol = 1e-07, btol = 1e-07,
               verbose = FALSE)
}
\arguments{
  \item{object}{
    the solution path object, of class "trendfilter", as returned by the
    \code{\link{trendfilter}} function.
  }
  \item{k}{
    an integer indicating the number of folds to split the data
    into. Must be between 2 and n-2 (n being the number of
    observations), default is 5. It is generally not a good idea to
    pass a value of k much larger than 10 (say, on the scale of n); see
    "Details" below. 
  }
  \item{mode}{
    a character string, either "lambda" or "df". Specifying "lambda"
    means that the cross-validation error will be computed and reported
    at each value of lambda that appears as a knot in the solution
    path. Specifying "df" means that the cross-validation error will be 
    computed and reported for every of degrees of freedom value
    (actually, estimate) incurred along the solution path. In the case
    that the same degrees of freedom value is visited multiple times,
    the model with the most regularization (smallest value of lambda) is
    considered. Default is "lambda". 
  }
  \item{approx}{
    a logical variable indicating if the approximate solution path
    should be used (with no dual coordinates leaving the boundary).
    Default is \code{FALSE}.
  }
  \item{rtol}{
    a numeric variable giving the relative tolerance used in the
    calculation of the hitting and leaving times. A larger value
    is more conservative, and may cause the algorithm to miss some
    hitting or leaving events (do not change unless you know what
    you're getting into!). Defaultis 1e-7. 
  }
  \item{btol}{
    similar to \code{rtol} but in absolute terms. If numerical
    instability is detected, first change rtol; then adjust btol
    if problems persist.
  }
  \item{verbose}{
    a logical variable indicating if progress should be reported after
    each knot in the path.
  }
}
\details{
  For trend filtering (with an identity predictor matrix), the folds 
  for k-fold cross-validation are chosen by placing every kth point into
  the same fold. (Here the points are implicitly ordered according to their 
  underlying positions---either assumed to be evenly spaced, or explicitly 
  passed through the \code{pos} argument.) 
  The first and last points are not included in any fold and are always
  included in building the predictive model. As an example, 
  with n=15 data points and k=4 folds, the points are assigned to folds 
  in the following way:
  \deqn{
    x \; 1 \; 2 \; 3 \; 4 \; 1 \; 2 \; 3 \;  4 \; 1 \; 2 \; 3 \; 4 \; 1
    \; x
  }{
    x 1 2 3 4 1 2 3 4 1 2 3 4 1 x
  }
  where \eqn{x} indicates no assignment. Therefore, the folds are not
  random and running \code{cv.trendfilter} twice will give the same
  result. In the calculation of the cross-validated error, the
  predicted value at a point is given by the average of the fits at this 
  point's two neighbors (guaranteed to be in a different fold). 

  Running cross-validation in modes "lambda" and "df" often yields very
  similar results. The mode "df" simply gives an alternative
  parametrization for the sequence of cross-validated models and can be
  more convenient for some applications; if you are confused about its
  function, simply leave the mode equal to "lambda".
}
\value{
  Returns and object of class "cv.trendfilter", a list with the
  following components:
  \item{err}{
    a numeric vector of cross-validated errors.
  }
  \item{se}{
    a numeric vector of standard errors (standard deviations of the
    cross-validation error estimates).
  }
  \item{mode}{
    a character string indicating the mode, either "lambda" or "df".
  }
  \item{lambda}{
    if \code{mode="lambda"}, the values of lambda at which the
    cross-validation errors in \code{err} were computed. 
  }
  \item{lambda.min}{
    if \code{mode="lambda"}, the value of lambda at which the
    cross-validation error is minimized.
  }
  \item{lambda.1se}{
    if \code{mode="lambda"}, the value of lambda chosen by the one
    standard error rule (the largest value of lambda such that the
    cross-validation error is within one standard error of the minimum).  
  }
  \item{df}{
    if \code{mode="df"}, the degrees of freedom values at which the
    cross-validation errors in \code{err} were computed.
  }
  \item{df.min}{
    if \code{mode="df"}, the degrees of freedom value at which the
    cross-validation error is minimized.
  }
  \item{df.1se}{
    if \code{mode="df"}, the degrees of freedom value chosen by the 
    one standard error rule (the smallest degrees of freedom value such 
    that cross-validation error is within one standard error of the
    minimum).   
  }
  \item{i.min}{
    the index of the model minimizing the cross-validation error.
  }
  \item{i.1se}{
    the index of the model chosen by the one standard error rule.
  }
  \item{call}{
    the matched call.
  }
}
\seealso{
  \code{\link{trendfilter}}, \code{\link{plot.cv.trendfilter}},
  \code{\link{plot.trendfilter}}
}
\examples{
# Constant trend filtering (the 1d fused lasso)
set.seed(0)
n = 50
beta0 = rep(sample(1:10,5),each=n/5)
y = beta0 + rnorm(n,sd=0.8)
a = fusedlasso1d(y)
plot(a)

# Choose lambda by 5-fold cross-validation
cv = cv.trendfilter(a)
plot(cv)
plot(a,lambda=cv$lambda.min,main="Minimal CV error")
plot(a,lambda=cv$lambda.1se,main="One standard error rule")

\dontrun{
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
a = trendfilter(y,ord=3,maxsteps=150)
plot(a,nlam=5)

# Choose lambda by 5-fold cross-validation
cv = cv.trendfilter(a)
plot(cv)
plot(a,lambda=cv$lambda.min,main="Minimal CV error")
plot(a,lambda=cv$lambda.1se,main="One standard error rule")
}
}
\keyword{utilities}
