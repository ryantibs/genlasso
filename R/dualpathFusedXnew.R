# We compute a solution path of the fused lasso dual problem:
#
# \hat{u}(\lambda) =
# \argmin_u \|y - (X^+)^T D^T u\|_2^2 \rm{s.t.} \|\u\|_\infty \leq \lambda
#
# where D is the incidence matrix of a given graph, and X^+ is the
# pseudoinverse of a full column rank predictor matrix X.
#
# Fortuitously, we never have to fully invert X (i.e. compute its 
# pseudoinverse).
#
# Note: the df estimates at each lambda_k can be thought of as the df
# for all solutions corresponding to lambda in (lambda_k,lambda_{k-1}),
# the open interval to the *right* of the current lambda_k.

dualpathFusedXnew <- function(y, X, D, approx=FALSE, maxsteps=2000, minlam=0,
                              tol=1e-11, verbose=FALSE, fileback=FALSE) {
  m = nrow(D)
  p = ncol(D)
  n = length(y)

  # Figure out whether or not we should be filebacking
  if (fileback!=FALSE) {
    if (fileback==TRUE) {
      fileback = paste("dualpathFused-output-",
        gsub(".","",format(Sys.time(),"%Y%m%d%H%M%OS2"),fixed=TRUE),
        ".csv", sep="")
    }
    if (!is.character(fileback) || length(fileback)>1) {
      stop("fileback must be either logical or a character string.")
    }
    zz = file(fileback, "w")
  }
  
  # Find the minimum 2-norm solution, using some linear algebra 
  # tricks and a little bit of graph theory
  L = abs(crossprod(D))
  diag(L) = 0
  gr = graph.adjacency(L,mode="undirected") # Underlying graph
  cl = clusters(gr)                         
  q = cl$no                                 # Number of clusters
  i = cl$membership                         # Cluster membership
  
  xy = t(X)%*%y
  A = matrix(0,n,q)
  z = numeric(q)
  for (j in Seq(1,q)) {
    oo = which(i==j)
    A[,j] = rowMeans(X[,oo,drop=FALSE])
    z[j] = mean(xy[oo])
  }

  # Compute a QR decomposition of A (we'll keep this)
  xx = qr(A)
  qrobj = list()
  Q = qr.Q(xx,complete=TRUE)           # n x n
  qrobj$Q1 = Q[,Seq(1,q),drop=FALSE]   # n x q
  qrobj$Q2 = Q[,Seq(q+1,n),drop=FALSE] # n x (n-q)
  qrobj$R = qr.R(xx)                   # q x q

  e = backsolve(qrobj$R,forwardsolve(qrobj$R,z,upper.tri=TRUE,transpose=TRUE))
  g = xy-t(X)%*%A%*%e

  x = f = numeric(p)
  
  # For efficiency, don't loop over singletons
  tab = tabulate(i)
  oo = which(tab[i]==1)
  if (length(oo)>0) {
    f[oo] = e[i][oo]
  }

  # Same for groups with two elements (doubletons?)
  oi = order(i)
  oo = which(tab[i][oi]==2)
  if (length(oo)>0) {
    f[oi][oo] = e[i][oi][oo]/2
    mm = colMeans(matrix(g[oi][oo],nrow=2))
    ii = oo[Seq(1,length(oo),by=2)]
    x[oi][ii] = g[oi][ii] - mm
  }

  # Now all groups with at least three elements
  cs = cumsum(tab)
  grps = which(tab>2)
  for (j in grps) {
    oo = oi[Seq(cs[j]-tab[j]+1,cs[j])]    
    f[oo] = e[j]/length(oo)
    gj = g[oo]
    Lj = crossprod(Matrix(D[,oo[-1]],sparse=TRUE))
    x[oo][-1] = as.numeric(solve(Lj,(gj-mean(gj))[-1]))
  }

  uhat = as.numeric(D%*%x)     # Dual solution
  betahat = f                  # Primal solution
  ihit = which.max(abs(uhat))  # Hitting coordinate
  hit = abs(uhat[ihit])        # Critical lambda
  s = sign(uhat[ihit])         # Sign

  if (verbose) {
    cat(sprintf("1. lambda=%.3f, adding coordinate %i, |B|=%i...",
                hit,ihit,1))
  }

  # Now iteratively find the new dual solution, and
  # the next critical point
  
  # Things to keep track of, and return at the end
  buf = min(maxsteps,1000)
  lams = numeric(buf)        # Critical lambdas
  h = logical(buf)           # Hit or leave?
  df = numeric(buf)          # Degrees of freedom
  lams[1] = hit
  h[1] = TRUE
  df[1] = q
  
  # We only record the solutions if there is no
  # filebacking
  if (fileback==FALSE) {
    u = matrix(0,m,buf)      # Dual solutions
    beta = matrix(0,p,buf)   # Primal solutions
    u[,1] = uhat
    beta[,1] = betahat
  }
  else {
    cat(m, p, file=zz, sep=",") 
    cat("\n", file=zz, sep="")
    cat(hit, TRUE, q, uhat, betahat, file=zz, sep=",")
    cat("\n", file=zz, sep="")
  }

  # Update our graph
  ed = which(D[ihit,]!=0)
  gr[ed[1],ed[2]] = 0             # Delete edge
  newcl = subcomponent(gr,ed[1])  # New cluster
  oldcl = which(i==i[ed[1]])      # Old cluster
  # If these two clusters aren't the same, update
  # the memberships
  if (length(newcl)!=length(oldcl) || any(sort(newcl)!=sort(oldcl))) {
    i[newcl] = q+1
    q = q+1
  }
  
  # Other things to keep track of, but not return
  r = 1                      # Size of boundary set
  B = ihit                   # Boundary set
  I = Seq(1,m)[-ihit]        # Interior set
  D1 = D[-ihit,,drop=FALSE]  # Matrix D[I,]
  D2 = D[ihit,,drop=FALSE]   # Matrix D[B,]
  k = 2                      # What step are we at?

  tryCatch({  
    while (k<=maxsteps && lams[k-1]>=minlam) {
      ##########
      # Check if we've reached the end of the buffer
      if (k > length(lams)) {
        buf = length(lams)
        lams = c(lams,numeric(buf))
        h = c(h,logical(buf))
        df = c(df,numeric(buf))
        if (fileback==FALSE) {
          u = cbind(u,matrix(0,m,buf))
          beta = cbind(beta,matrix(0,p,buf))
        }
      }
      
      ##########
      Ds = as.numeric(t(D2)%*%s)

      # Precomputation for the hitting times
  #    A = matrix(0,n,q)
      z = matrix(0,q,2)
      for (j in Seq(1,q)) {
        oo = which(i==j)
        A[,j] = rowMeans(X[,oo,drop=FALSE])
        z[j,1] = mean(xy[oo])
        z[j,2] = mean(Ds[oo])
      }
      e = backsolve(qrobj$R,forwardsolve(qrobj$R,z,upper.tri=TRUE,transpose=TRUE))
      ea = e[,1]
      eb = e[,2]
      ga = xy-t(X)%*%A%*%ea
      gb = Ds-t(X)%*%A%*%eb
      
      # If the interior is empty, then nothing will hit
      if (r==m) {
        fa = ea[i]
        fb = eb[i]
        hit = 0
      }
      
      # Otherwise, find the next hitting time
      else {            
        xa = xb = numeric(p)
        fa = fb = numeric(p)

        # For efficiency, don't loop over singletons
        tab = tabulate(i)
        oo = which(tab[i]==1)
        if (length(oo)>0) {
          fa[oo] = ea[i][oo]
          fb[oo] = eb[i][oo]
        }

        # Same for groups with two elements (doubletons?)
        oi = order(i)
        oo = which(tab[i][oi]==2)
        if (length(oo)>0) {
          fa[oi][oo] = ea[i][oi][oo]/2
          fb[oi][oo] = eb[i][oi][oo]/2
          ma = colMeans(matrix(ga[oi][oo],nrow=2))
          mb = colMeans(matrix(gb[oi][oo],nrow=2))
          ii = oo[Seq(1,length(oo),by=2)]
          xa[oi][ii] = ga[oi][ii] - ma
          xb[oi][ii] = gb[oi][ii] - mb
        }

        # Now all groups with at least three elements
        cs = cumsum(tab)
        grps = which(tab>2)
        for (j in grps) {
          oo = oi[Seq(cs[j]-tab[j]+1,cs[j])]    
          fa[oo] = ea[j]/length(oo)
          fb[oo] = eb[j]/length(oo)
          gaj = ga[oo]
          gbj = gb[oo]
          Lj = crossprod(Matrix(D1[,oo[-1]],sparse=TRUE))
          xa[oo][-1] = as.numeric(solve(Lj,(gaj-mean(gaj))[-1]))
          xb[oo][-1] = as.numeric(solve(Lj,(gbj-mean(gbj))[-1]))
        }
        
        a = as.numeric(D1%*%xa)
        b = as.numeric(D1%*%xb)
        shits = Sign(a)
        hits = a/(b+shits);

        # Make sure none of the hitting times are larger
        # than the current lambda (precision issue)
        hits[hits>lams[k-1]] = hits[k-1]
        
        ihit = which.max(hits)
        hit = hits[ihit]
        shit = shits[ihit]
      }
      
      ##########
      # If nothing is on the boundary, then nothing will leave
      # Also, skip this if we are in "approx" mode
      if (r==0 || approx) {
        leave = 0
      }

      # Otherwise, find the next leaving time
      else {
        c = as.numeric(s*(D2%*%fa))
        d = as.numeric(s*(D2%*%fb))
        leaves = c/d

        # c must be negative
        leaves[c>=0] = 0

        # Make sure none of the leaving times are larger
        # than the current lambda (precision issue)
        leaves[leaves>lams[k-1]] = lams[k-1]

        ileave = which.max(leaves)
        leave = leaves[ileave]
      }

      ##########
      # Stop if the next critical point is negative
      if (hit<=0 && leave<=0) break

      # If a hitting time comes next
      if (hit > leave) {
        # Record the critical lambda and properties
        lams[k] = hit
        h[k] = TRUE
        df[k] = q
        uhat = numeric(m)
        uhat[B] = hit*s
        uhat[I] = a-hit*b
        betahat = fa-hit*fb

        # Only record the solutions if we are not
        # filebacking
        if (fileback==FALSE) {
          u[,k] = uhat
          beta[,k] = betahat
        }
        else {
          cat(hit, TRUE, q, uhat, betahat, file=zz, sep=",")
          cat("\n", file=zz, sep="")
        }      

        # Update our graph
        ed = which(D1[ihit,]!=0)
        gr[ed[1],ed[2]] = 0             # Delete edge
        newcl = subcomponent(gr,ed[1])  # New cluster
        oldcl = which(i==i[ed[1]])      # Old cluster
        # If these two clusters aren't the same, update
        # the memberships
        if (length(newcl)!=length(oldcl) || any(sort(newcl)!=sort(oldcl))) {
          # Shift the two new clusters to the end
          oldno = i[ed[1]]
          i[i>oldno] = i[i>oldno]-1
          i[oldcl] = q
          i[newcl] = q+1

          # QR update: first remove one column, then add two
          qrobj = downdateW(qrobj$Q1,qrobj$Q2,qrobj$R,oldno)
          col1 = rowMeans(X[,which(i==q),drop=FALSE])
          col2 = rowMeans(X[,which(i==q+1),drop=FALSE])
          qrobj = updateW(qrobj$Q1,qrobj$Q2,qrobj$R,col1)
          qrobj = updateW(qrobj$Q1,qrobj$Q2,qrobj$R,col2)

          q = q+1
        }
        
        # Update all other variables
        r = r+1
        B = c(B,I[ihit])
        I = I[-ihit]
        s = c(s,shit)
        D2 = rBind(D2,D1[ihit,])
        D1 = D1[-ihit,,drop=FALSE]
          
        if (verbose) {
          cat(sprintf("\n%i. lambda=%.3f, adding coordinate %i, |B|=%i...",
                      k,hit,B[r],r))
        }
      }
                
      # Otherwise a leaving time comes next
      else {
        # Record the critical lambda and properties
        lams[k] = leave
        h[k] = FALSE
        df[k] = q
        uhat = numeric(m)
        uhat[B] = leave*s
        uhat[I] = a-leave*b
        betahat = fa-leave*fb
        
        # Only record the solutions if we are not
        # filebacking
        if (fileback==FALSE) {
          u[,k] = uhat
          beta[,k] = betahat
        }
        else {
          cat(leave, FALSE, q, uhat, betahat, file=zz, sep=",")
          cat("\n", file=zz, sep="")
        }

        # Update our graph
        ed = which(D2[ileave,]!=0)
        gr[ed[1],ed[2]] = 1             # Add edge
        newcl = subcomponent(gr,ed[1])  # New cluster
        oldcl = which(i==i[ed[1]])      # Old cluster
        # If these two clusters aren't the same, update
        # the memberships
        if (length(newcl)!=length(oldcl) || !all(sort(newcl)==sort(oldcl))) {
          # Shift the new cluster to the end
          oldno1 = max(i[ed[1]],i[ed[2]])
          oldno2 = min(i[ed[1]],i[ed[2]])
          i[i==oldno1] = oldno2
          i[i>oldno1] = i[i>oldno1]-1
          i[i==oldno2] = q
          i[i>oldno2] = i[i>oldno2]-1
          
          # QR update: first remove two columns, then add one
          qrobj = downdateW(qrobj$Q1,qrobj$Q2,qrobj$R,oldno1)
          qrobj = downdateW(qrobj$Q1,qrobj$Q2,qrobj$R,oldno2)
          col = rowMeans(X[,which(i==q-1),drop=FALSE])
          qrobj = updateW(qrobj$Q1,qrobj$Q2,qrobj$R,col)
          
          q = q-1
        }
        
        # Update all other variables
        r = r-1
        I = c(I,B[ileave])
        B = B[-ileave]
        s = s[-ileave]
        D1 = rBind(D1,D2[ileave,])
        D2 = D2[-ileave,,drop=FALSE]

        if (verbose) {
          cat(sprintf("\n%i. lambda=%.3f, deleting coordinate %i, |B|=%i...",
                      k,leave,I[m-r],r))
        }
      }
      
      # Step counter
      k = k+1
    }
  }, error = function(err) {
    err$message = paste(err$message,"\n(Path computation has been terminated;",
      " partial path is being returned.)",sep="")
    warning(err)})

  # Trim 
  lams = lams[Seq(1,k-1)]
  h = h[Seq(1,k-1)]
  df = df[Seq(1,k-1)]
  if (fileback==FALSE) {
    u = u[,Seq(1,k-1),drop=FALSE]
    beta = beta[,Seq(1,k-1),drop=FALSE]
  }

  # If we reached the maximum number of steps
  if (k>maxsteps) {
    if (verbose) {
      cat(sprintf("\nReached the max number of steps (%i),",maxsteps))
      cat(" skipping the rest of the path.")
    }
    completepath = FALSE
  }

  # If we reached the minimum lambda
  else if (lams[k-1]<minlam) {
    if (verbose) {
      cat(sprintf("\nReached the min lambda (%.3f),",minlam))
      cat(" skipping the rest of the path.")
    }
    completepath = FALSE
  }

  # Otherwise, note that we completed the path
  else completepath = TRUE
  
  if (verbose) cat("\n")

  if (fileback==FALSE) {
    colnames(u) = as.character(round(lams,3))
    colnames(beta) = as.character(round(lams,3))
    return(list(lambda=lams,beta=beta,fit=X%*%beta,u=u,hit=h,df=df,y=y,X=X,
                completepath=completepath,bls=if (completepath) fa else NULL))
  }
  else {
    close(zz)
    cat("\nPath object can be read in by calling filebackout on the function output.\n")
    invisible(list(y=y,X=X,file=fileback))
  }
}
