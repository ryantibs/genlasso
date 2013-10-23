filebackout <- function(object, y, X, file) {
  if (missing(X)) X = NULL
  
  if (!missing(object)) {
    file = object$file
    y = object$y
    X = object$X
  }
  else if (missing(y) || missing(file)) {
    stop(paste("Must either specify a filebacked object, or y, X, and the file name (but",
               "X does not need to be passed if the predictor matrix is the identity)."))
  }
  
  fi = file.info(file)
  if(is.na(fi$isdir)) stop(sprintf("Invalid file name \"%s\": does not exist in working directory.",file))
  if(fi$isdir) stop(sprintf("Invalid file name \"%s\": points to a directory.",file))

  a = read.csv(file,header=FALSE,as.is=TRUE)
  k = nrow(a)
  m = as.integer(a[1,1])
  n = as.integer(a[1,2])
  gamma = as.numeric(a[1,3])
  lambda = as.numeric(a[-c(1,k),1])
  hit = as.logical(a[-c(1,k),2])
  df = as.integer(a[-c(1,k),3])
  u = t(as.matrix(a[-c(1,k),3+Seq(1,m)]))
  beta = t(as.matrix(a[-c(1,k),3+m+Seq(1,n)]))
  rownames(u) = rownames(beta) = NULL 
  colnames(u) = colnames(beta) = as.character(round(lambda,3))
  completepath = as.logical(a[k,1])
  if (as.logical(a[k,2])) bls = NULL
  else bls = as.numeric(a[k,2+Seq(1,n)])
  
  out = list(lambda=lambda,beta=beta,fit=beta,u=u,hit=hit,df=df,y=y,
    completepath=completepath,bls=bls,gamma=gamma,call=match.call())
  if (!is.null(X)) out$X = X
  class(out) = c("fusedlasso", "genlasso", "list")
  return(out)
}
