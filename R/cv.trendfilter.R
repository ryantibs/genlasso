cv.trendfilter <- function(object, k=10, mode=c("lambda", "df"),
                           approx=FALSE, tol=1e-11, verbose=FALSE) {
  cl = match.call()
  
  if (all(class(object)!="trendfilter")) {
    stop("Cross-validation can only be performed for trend filtering.")
  }
  if (!is.null(object$X)) {
    stop("Cross-validation for trend filtering can only be performed when X=I, the identity matrix.")
  }
#  if (!is.null(object$z)) {
#    stop("Cross-validation for trend filtering can only be performed on a unit spaced grid.")
#  }
  mode = mode[[1]]
  if (!(mode %in% c("lambda", "df"))) {
    stop("Invalid mode, must be \"lambda\" or \"df\".")
  }

  y = object$y
  n = length(y)

  if(k<2 || round(k)!=k || k>n-2) {
    stop("The number of folds must an integer between 2 and n-2.")
  }
  
  ord = object$trendorder
  foldid = c(0,rep(Seq(1,k),n-2)[Seq(1,n-2)],0)

  if (mode=="lambda") { 
    lambda = object$lambda
    cvall = matrix(0,k,length(lambda))

    for (i in Seq(1,k)) {
      cat(sprintf("Fold %i ... ",i))
      
      otr = which(foldid!=i)
      ntr = length(otr)
      ytr = y[otr]
      Dtr = getDtfSparse(ntr,ord)
      outtr = dualpathWideSparse(ytr,Dtr,NULL,approx,Inf,min(lambda)*ntr/n,tol,verbose)
      b = coef.genlasso(outtr,lambda=lambda*ntr/n)$beta
      
      ote = which(foldid==i)
      yte = matrix(y[ote],length(ote),length(lambda))
      pred = (b[which((Seq(1,n)%in%(ote-1))[otr]),] +
              b[which((Seq(1,n)%in%(ote+1))[otr]),])/2
      
      cvall[i,] = colMeans((yte-pred)^2)
    }

    cverr = colMeans(cvall)
    cvse = apply(cvall,2,sd)/sqrt(k)
    
    names(cverr) = names(cvse) = round(lambda,3)
    i0 = which.min(cverr)
    lam.min = lambda[i0]
    lam.1se = max(lambda[cverr<=cverr[i0]+cvse[i0]])
    i.min = which(lambda==lam.min)
    i.1se = which(lambda==lam.1se)
    
    out = list(err=cverr,se=cvse,mode="lambda",lambda=lambda,
      lambda.min=lam.min,lambda.1se=lam.1se,i.min=i.min,i.1se=i.1se,call=cl)
  }

  else {
    df = object$df
    # Only keep the dfs for which we have at least that
    # many points in any configuration of k-1 folds
    df = unique(df[df <= n-max(tabulate(foldid))])
    cvall = matrix(0,k,length(df))

    for (i in Seq(1,k)) {
      cat(sprintf("Fold %i ... ",i))
      
      otr = which(foldid!=i)
      ntr = length(otr)
      ytr = y[otr]
      Dtr = getDtfSparse(ntr,ord)
      outtr = dualpathWideSparse(ytr,Dtr,NULL,approx,Inf,0,tol,verbose)
      b = coef.genlasso(outtr,df=df)$beta

      ote = which(foldid==i)
      yte = matrix(y[ote],length(ote),length(df))
      pred = (b[which((Seq(1,n)%in%(ote-1))[otr]),] +
              b[which((Seq(1,n)%in%(ote+1))[otr]),])/2

      cvall[i,] = colMeans((yte-pred)^2)
    }

    cverr = colMeans(cvall)
    cvse = apply(cvall,2,sd)/sqrt(k)

    names(cverr) = names(cvse) = df
    i0 = which.min(cverr)
    df.min = df[i0]
    df.1se = min(df[cverr<=cverr[i0]+cvse[i0]])
    i.min = max(which(object$df==df.min))
    i.1se = max(which(object$df==df.1se))

    out = list(err=cverr,se=cvse,mode="df",df=df,
      df.min=df.min,df.1se=df.1se,i.min=i.min,i.1se=i.1se,call=cl)
  }

  cat("\n")
  class(out) = c("cv.trendfilter", "list")  
  return(out)
}
