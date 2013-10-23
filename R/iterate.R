iterate <- function(object, maxsteps=2000, minlam=0, verbose=FALSE) {

  cl = match.call()

  # Cannot iterate if path completed:
  if(object$completepath) stop("Path has completed, nothing to iterate!")

  # Assuming for now we have dualpathTall
  implementedValues = c("tall", "wide", "wide_sparse", "fused", "fusedl1", "fusedx", "fusedl1x")
  if(!(object$dualpathObjs$dualpathType %in% implementedValues)) stop("not yet implemented")

  # Iterating dualpathTall:
  if(object$dualpathObjs$dualpathType == "tall") {
    
    output_dp = dualpathTall(object=object, maxsteps=maxsteps, minlam=minlam, verbose=verbose)

  }
  # Iterating dualpathWide:
  if(object$dualpathObjs$dualpathType == "wide") {

    output_dp = dualpathWide(object=object, maxsteps=maxsteps, minlam=minlam, verbose=verbose)

  }
  # Iterating dualpathWideSparse:
  if(object$dualpathObjs$dualpathType == "wide_sparse") {
    
    output_dp = dualpathWideSparse(object=object, maxsteps=maxsteps, minlam=minlam, verbose=verbose)

  }

  # Iterating for dualpathFused:
  if(object$dualpathObjs$dualpathType == "fused") {
    
    output_dp = dualpathFused(object=object, maxsteps=maxsteps, minlam=minlam, verbose=verbose)
    output_dp$call = cl
    class(output_dp) = c("fusedlasso", "genlasso", "list")
    return(output_dp)

  }
  # Iterating for dualpathFusedL1:
  if(object$dualpathObjs$dualpathType == "fusedl1") {
    
    output_dp = dualpathFusedL1(object=object, maxsteps=maxsteps, minlam=minlam, verbose=verbose)
    output_dp$call = cl
    class(output_dp) = c("fusedlasso", "genlasso", "list")
    return(output_dp)

  }
  # Iterating for dualpathFusedX:
  if(object$dualpathObjs$dualpathType == "fusedx") {
    
    output_dp = dualpathFusedX(object=object, maxsteps=maxsteps, minlam=minlam, verbose=verbose)
    output_dp$call = cl
    class(output_dp) = c("fusedlasso", "genlasso", "list")
    return(output_dp)

  }
  # Iterating for dualpathFusedL1X:
  if(object$dualpathObjs$dualpathType == "fusedl1x") {
    
    output_dp = dualpathFusedL1X(object=object, maxsteps=maxsteps, minlam=minlam, verbose=verbose)
    output_dp$call = cl
    class(output_dp) = c("fusedlasso", "genlasso", "list")
    return(output_dp)

  }


  # Create output object (work usually done by genlasso/dualpath/other
  #   high-level functions
  lambda = Xi = X = y2 = x = n0 = n = y=  y0 = D0 = NULL
  dualpathObjs = object$dualpathObjs
  for(j in 1:length(object)) assign(names(object)[j], object[[j]])
  for(j in 1:length(dualpathObjs)) assign(names(dualpathObjs)[j], dualpathObjs[[j]])
  lams = lambda
  if(is.null(n0)) n0 = n
  if(is.null(y0)) y0 = y
  if(is.null(D0)) D0 = D

  out = object
  out$lambda = output_dp$lambda
  out$u = output_dp$u
  out$hit = output_dp$hit
  out$df = output_dp$df

  if(sum(is.na(object$dualpathObjs$X)) == 0) {

    # If there is an explicit X matrix
    out$beta = Xi %*% output_dp$fit
    out$fit  = X %*% out$beta
    out$y    = object$y
    out$bls  = Xi %*% y2
    out$X    = X

    out$dualpathObjs = output_dp$dualpathObjs
    out$dualpathObjs$X  = X
    out$dualpathObjs$n0 = object$dualpathObjs$n0
    out$dualpathObjs$y0 = object$dualpathObjs$y0
    out$dualpathObjs$j  = object$dualpathObjs$j
    out$dualpathObjs$D0 = object$dualpathObjs$D0
    out$dualpathObjs$x  = x
    out$dualpathObjs$y2 = y2
    out$dualpathObjs$Xi = Xi

  } else {

    # If there is no explicit X matrix
    out$completepath = output_dp$completepath
    out$dualpathObjs = output_dp$dualpathObjs
    out$df = out$df + n0 - n
    beta = matrix(y0,n0,length(out$lambda))
    j = object$dualpathObjs$j
    if(!is.null(j)) beta[j,] = as.matrix(y0[j] - t(D0[,j])%*%out$u)  

    colnames(beta) = colnames(output_dp$u)
    out$beta = beta
    out$fit = beta
    out$y = y0
    out$bls = y0

    if(is.null(out$dualpathObjs$X)) out$dualpathObjs$X = NA
    out$dualpathObjs$n0 = object$dualpathObjs$n0
    out$dualpathObjs$y0 = object$dualpathObjs$y0
    out$dualpathObjs$j  = object$dualpathObjs$j
    out$dualpathObjs$D0 = object$dualpathObjs$D0

  }

  out$call = cl
  class(out) = class(object)

  return(out)
}
