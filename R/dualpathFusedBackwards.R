# We compute a solution path of the fused lasso dual problem:
#
# \hat{u}(\lambda) =
# \argmin_u \|y - D^T u\|_2^2 \rm{s.t.} \|\u\|_\infty \leq \lambda
#
# where D is the incidence matrix of a given graph.
# 
# The main difference here is that we compute the path backwards,
# i.e. starting from lambda=0 rather than lambda=infinity.
#
# Note: the df estimates at each lambda_k can be thought of as the df
# for all solutions corresponding to lambda in (lambda_k,lambda_{k+1}),
# the open interval to the *right* of the current lambda_k.

dualpathFusedBackwards <- function(y, D, approx=FALSE, maxsteps=2000, maxlam=Inf,
                                   tol=1e-11, verbose=FALSE, fileback=FALSE) {
  m = nrow(D)
  n = ncol(D)

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

  # Find the first hitting time (all of the terminology here
  # refers to the idea of decreasing lambda, even though we're
  # actually increasing lambda; so a hitting event is one in
  # which a coordinate would have hit the boundary if we think
  # of decreasing lambda, so it would have left the boundary if
  # we think of increasing lambda)
  s = as.numeric(sign(D%*%y))
  Ds = as.numeric(t(D)%*%s)
  c = as.numeric(s*(D%*%y))
  d = as.numeric(s*(D%*%Ds))
  hits = c/d

  # d must be positive
  hits[d<=0] = Inf
  
  # Make sure none of the hitting times are smaller
  # than the current lambda (precision issue)
  hits[hits<0] = 0
  
  ihit = which.min(hits)
  hit = hits[ihit]

  if (verbose) {
    cat(sprintf("1. lambda=%.3f, adding coordinate %i, |B|=%i...",
                hit,ihit,m-1))
  }

  # Build a graph
  e = which(D[ihit,]!=0)
  gr = graph(edges=e,n=n,directed=FALSE)  # Underlying graph
  cl = clusters(gr)                         
  q = cl$no                               # Number of clusters
  i = cl$membership                       # Cluster membership

  uhat = hit*s                  # Dual solution
  betahat = y-hit*Ds            # Primal solution
  betahat[e] = mean(betahat[e])
  
  # Things to keep track of, and return at the end
  buf = min(maxsteps,1000)
  lams = numeric(buf)        # Critical lambdas
  h = logical(buf)           # Hit or leave?
  df = numeric(buf)          # Degrees of freedom
  
  lams[1] = hit
  h[1] = TRUE
  df[1] = n-1

  # We only record the solutions if there is no
  # filebacking
  if (fileback==FALSE) {
    u = matrix(0,m,buf)      # Dual solutions
    beta = matrix(0,n,buf)   # Primal solutions
    u[,1] = uhat
    beta[,1] = betahat
  }
  else {
    cat(m, n, file=zz, sep=",")
    cat("\n", file=zz, sep="")
    cat(hit, TRUE, n-1, uhat, betahat, file=zz, sep=",")
    cat("\n", file=zz, sep="")
  }
  
  # Other things to keep track of
  r = m-1                    # Size of boundary set
  B = Seq(1,m)[-ihit]        # Boundary set
  I = ihit                   # Interior set
  s = s[-ihit]               # Boundary signs
  D1 = D[ihit,,drop=FALSE]   # Matrix D[I,]
  D2 = D[-ihit,,drop=FALSE]  # Matrix D[B,]
  k = 2                      # What step are we at?

  tryCatch({  
    while (k<=maxsteps && lams[k-1]<=maxlam) {
      ##########
      # Check if we've reached the end of the buffer
      if (k > length(lams)) {
        buf = length(lams)
        lams = c(lams,numeric(buf))
        h = c(h,logical(buf))
        df = c(df,numeric(buf))
        if (fileback==FALSE) {
          u = cbind(u,matrix(0,m,buf))
          beta = cbind(beta,matrix(0,n,buf))
        }
      }

      ##########
      Ds = as.numeric(t(D2)%*%s)
     
      # If the boundary set is full, then nothing could have left
      # Also, skip this if we are in "approx" mode
      if (r==m || approx) {
        fa = y
        fb = Ds
        leave = Inf
      }
      
      # Otherwise, find the next leaving time
      else {
        xa = xb = numeric(n)
        fa = fb = numeric(n)
        
        # For efficiency, don't loop over singletons
        tab = tabulate(i)
        oo = which(tab[i]==1)
        if (length(oo)>0) {
          fa[oo] = y[oo]
          fb[oo] = Ds[oo]
        }
        
        # Same for groups with two elements (doubletons?)
        oi = order(i)
        oo = which(tab[i][oi]==2)
        if (length(oo)>0) {
          ma = colMeans(matrix(y[oi][oo],nrow=2))
          mb = colMeans(matrix(Ds[oi][oo],nrow=2))
          fa[oi][oo] = rep(ma,each=2)
          fb[oi][oo] = rep(mb,each=2)
          ii = oo[Seq(1,length(oo),by=2)]
          xa[oi][ii] = y[oi][ii] - ma
          xb[oi][ii] = Ds[oi][ii] - mb
        }

        # Now all groups with at least three elements
        cs = cumsum(tab)
        grps = which(tab>2)
        for (j in grps) {
          oo = oi[Seq(cs[j]-tab[j]+1,cs[j])]
          yj = y[oo]
          Dsj = Ds[oo]
          fa[oo] = mean(yj)
          fb[oo] = mean(Dsj)
          Lj = crossprod(Matrix(D1[,oo[-1]],sparse=TRUE))
          xa[oo][-1] = as.numeric(solve(Lj,(yj-mean(yj))[-1]))
          xb[oo][-1] = as.numeric(solve(Lj,(Dsj-mean(Dsj))[-1]))
        }
        
        a = as.numeric(D1%*%xa)
        b = as.numeric(D1%*%xb)
        sleaves = Sign(-b)
        leaves = a/(b+sleaves);

        # Make sure none of the leaving times are smaller
        # than the current lambda (not just precision!)
        leaves[abs(b)<=1] = Inf
        leaves[leaves<lams[k-1]] = lams[k-1]
        
        ileave = which.min(leaves)
        leave = leaves[ileave]
        sleave = sleaves[ileave]
      }

      ##########
      # If nothing is on the boundary, then nothing could have hit
      if (r==0) {
        hit = Inf
      }
      
      # Otherwise find the next hitting time
      else {
        c = as.numeric(s*(D2%*%fa))
        d = as.numeric(s*(D2%*%fb))
        hits = c/d

        # d must be positive
        hits[d<=0] = Inf
      
        # Make sure none of the leaving times are smaller
        # than the current lambda (precision issue)
        hits[hits<lams[k-1]] = lams[k-1]
      
        ihit = which.min(hits)
        hit = hits[ihit]
      }
     
      ##########
      # Stop if the next critical point is infinite
      if (hit==Inf && leave==Inf) break

      # If a hitting time comes next
      if (hit < leave) {
        # Update our graph
        e = which(D2[ihit,]!=0)
        gr[e[1],e[2]] = 1             # Add edge
        newcl = subcomponent(gr,e[1]) # New cluster
        oldcl = which(i==i[e[1]])     # Old cluster
        # If these two clusters aren't the same, update
        # the memberships
        if (length(newcl)!=length(oldcl) || any(sort(newcl)!=sort(oldcl))) {
          newno = i[e[2]]
          oldno = i[e[1]]
          i[oldcl] = newno
          i[i>oldno] = i[i>oldno]-1
          q = q-1
        }

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

        # Update all other variables
        r = r-1
        I = c(I,B[ihit])
        B = B[-ihit]
        s = s[-ihit]
        D1 = rBind(D1,D2[ihit,])
        D2 = D2[-ihit,,drop=FALSE]

        if (verbose) {
          cat(sprintf("\n%i. lambda=%.3f, adding coordinate %i, |B|=%i...",
                      k,hit,I[m-r],r))
        }
      }
                
      # Otherwise a leaving time comes next
      else {
        # Update our graph
        e = which(D1[ileave,]!=0)
        gr[e[1],e[2]] = 0             # Delete edge
        newcl = subcomponent(gr,e[1]) # New cluster
        oldcl = which(i==i[e[1]])     # Old cluster
        # If these two clusters aren't the same, update
        # the memberships
        if (length(newcl)!=length(oldcl) || any(sort(newcl)!=sort(oldcl))) {
          i[newcl] = q+1
          q = q+1
        }
        
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
        
        # Update all other variables
        r = r+1
        B = c(B,I[ileave])
        I = I[-ileave]
        s = c(s,sleave)
        D2 = rBind(D2,D1[ileave,])
        D1 = D1[-ileave,,drop=FALSE]
          
        if (verbose) {
          cat(sprintf("\n%i. lambda=%.3f, deleting coordinate %i, |B|=%i...",
                      k,leave,B[r],r))
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
      cat(sprintf("\nReached the maximum number of steps (%i),",maxsteps))
      cat(" skipping the rest of the path.")
    }
    completepath = FALSE
  }

  # If we reached the maximum lambda
  else if (lams[k-1]>maxlam) {
    if (verbose) {
      cat(sprintf("\nReached the maximum lambda (%.3f),",maxlam))
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
    return(list(lambda=lams,beta=beta,fit=beta,u=u,hit=h,df=df,y=y,
                completepath=completepath,bls=y))
  }
  else {
    close(zz)
    cat("\nPath object can be read in by calling filebackout on the function output.\n")
    invisible(list(y=y,file=fileback))
  }
}
