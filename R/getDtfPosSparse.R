getDtfPosSparse <- function(n,ord,z) {
  weights = 1/diff(c(z,z[n]+1))
  D = bandSparse(n, m=n, c(0,1), diagonals=list(-weights,weights))
  D0 = bandSparse(n, m=n, c(0,1), diagonals=list(rep(-1,n),rep(1,n-1)))
  for (i in Seq(1,ord)) {
    D = D0 %*% D
  }
  return(D[1:(n-ord-1),])
}
