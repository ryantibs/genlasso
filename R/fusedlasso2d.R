# We compute the solution path of the fused lasso problem:
#
# \hat{\beta}(\lambda) =
# \argmin_\beta \|y - X \beta|_2^2 + \lambda\|D \beta\|_1,
#
# where D is the incidence matrix of a 2d grid and X is a full
# column rank predictor matrix. The solution is piecewise constant
# over the grid.

fusedlasso2d <- function(y, X, dim1, dim2, gamma=0, approx=FALSE, 
                         maxsteps=2000, minlam=0, tol=1e-11, verbose=FALSE,
                         fileback=FALSE) {
  if (missing(y)) stop("y is missing.")
  if (!is.numeric(y) && !is.matrix(y)) stop("y must be numeric (if X is missing, then y can alternatively be a matrix).")
  if (length(y) == 0) stop("There must be at least one data point [must have length(y) > 1].")
  if (missing(X)) X = NULL
  if (missing(dim1) || missing(dim2)) {
    if (is.matrix(y) && is.null(X)) {
      dim1 = nrow(y)
      dim2 = ncol(y)
    }
    else {
      stop("Both dim1 and dim2 must be specified.")
    }
  }
  else if (dim1<0 || round(dim1)!=dim1 || dim2<0 || round(dim2)<0) {
    stop("Both dim1 and dim2 must be nonnegative integers.")
  }
  if (is.null(X) && length(y)!=dim1*dim2) {
    stop("Dimensions don't match [length(y) != dim1*dim2].")
  }
  if (!is.null(X) && ncol(X)!=dim1*dim2) {
    stop("Dimensions don't match [ncol(X) != dim1*dim2].")
  }
  
  y = as.numeric(y)
  D = getD2dSparse(dim1,dim2)

  out = fusedlasso(y,X,D,NULL,gamma,approx,maxsteps,minlam,tol,verbose,fileback)
  out$call = match.call()
  
  if (fileback==FALSE) return(out)
  else invisible(out)
}
