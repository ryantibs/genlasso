\name{plot.genlasso}
\alias{plot.genlasso}
\alias{plot.trendfilter}
\alias{plot.cv.trendfilter}
\title{
  Plotting methods for generalized lasso objects
}
\description{
  The function \code{plot.genlasso} produces a plot of the coordinate 
  paths for objects of class "genlasso". This can be helpful for
  visualizing the full solution path for small problems; however, for
  moderate or large problems, the plot produced can be quite dense and
  difficult to interpret. The function \code{plot.trendfilter} applies
  to objects of class "trendfilter", and plots trend filtering
  coefficients at a single value of lambda (or multiple 
  values, as specified by the user) as a function of the input positions
  (which, recall, are assumed to be evenly spaced if not specified).
  The function \code{plot.cv.trendfilter} plots the output of
  \code{cv.trendfilter}. 
}
\usage{
\method{plot}{genlasso}(x, type = c("primal", "dual", "both"), numbers = FALSE,
     vlines = TRUE, xlab, ylab, ...)  
\method{plot}{trendfilter}(x, style = c("trend", "path"), lambda, nlam, df, xlab,
     ylab, ...)
\method{plot}{cv.trendfilter}(x, legendpos = "top", xlab, ylab, ...) 
}
\arguments{
  \item{x}{
    an object of the appropriate class ("genlasso" or anything class
    inherits this for \code{plot.genlasso}, "trendfilter" for
    \code{plot.trendfilter}, and "cv.trendfilter" for
    \code{plot.cv.trendfilter}). 
  }
  \item{type}{
    for \code{plot.genlasso}, a character string, one of "primal",
    "dual", or "both", indicating which solution path system(s) should
    be plotted. Default is "primal". 
  }
  \item{numbers}{
    for \code{plot.genlasso},
    a logical variable indicating if coordinate paths should be labeled
    by their numbers. Default is \code{FALSE}.
  }
  \item{vlines}{
    for \code{plot.genlasso},
    a logical variable indicating if dashed lines should be drawn at
    knots in the path, with black lines for hitting events, and red
    lines for leaving events. Default is \code{TRUE}.
  }
  \item{style}{
    for \code{plot.trendfilter},
    a character string, either "trend" or "path". If "trend", then trend
    filtering coefficients are plotted according to their underlying 
    positions. If
    "path", then a plot of the coordinate paths is produced with the
    function \code{plot.genlasso}. Default is "trend".
  }
  \item{lambda, nlam, df}{
    for \code{plot.trendfilter},
    these arguments work exactly as they do in
    \code{\link{coef.genlasso}}, and they are used to specify which
    solutions should be extracted and plotted from the computed solution
    path stored in \code{x}. The only difference is, if all three
    are missing, then \code{nlam} is set to 10 (whereas in
    \code{coef.genlasso}, the default is to set \code{lambda} equal to
    the full set of knots along the solution path).
  }
  \item{legendpos}{
    for \code{plot.cv.trendfilter}, a character string indicating the
    position of the legend. Default is "top". 
  }
  \item{xlab}{
    an optional character string label for the x-axis. 
  }
  \item{ylab}{
    an optional character string label for the y-axis.
  }
  \item{\ldots}{
    additional arguments.
  }
}
\value{
  For \code{plot.trendfilter}, with \code{style} set to "trend", a
  coefficient object is silently returned as specified by \code{lambda},
  \code{nlam}, or \code{df}.
}
\seealso{
  \code{\link{genlasso}}, \code{\link{trendfilter}},
  \code{\link{cv.trendfilter}} 
}
\examples{
# Constant trend filtering (the 1d fused lasso)
set.seed(0)
n = 100
beta0 = rep(sample(1:10,5),each=n/5)
y = beta0 + rnorm(n,sd=0.8)
a = fusedlasso1d(y)
cv = cv.trendfilter(a)

plot(a,style="path")
plot(cv)
plot(a,lambda=cv$lambda.1se)
}
\keyword{hplot}
