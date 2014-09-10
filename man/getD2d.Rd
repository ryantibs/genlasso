\name{getDxx}
\alias{getD1dSparse}
\alias{getD2dSparse}
\alias{getDgSparse}
\alias{getDtfSparse}
\alias{getDtfPosSparse}
\alias{getD1d}
\alias{getD2d}
\alias{getDg}
\alias{getDtf}
\alias{getDtfPos}
\alias{getGraph}
\title{
  Helper functions for constructing generalized lasso penalty matrices 
}
\description{	
  These are utility functions for creating penalty matrices for the
  fused lasso and trend filtering problems. Most users will not need to
  explicitly construct these as they are created internally by the
  \code{fusedlasso} or \code{trendfilter} functions. The sparse
  variants output sparse matrices, which should be used whenever
  possible because of a significant savings in both construction speed
  and memory usage.

  The function getGraph is an inverse function for fused lasso problems,
  returning an \code{igraph} object (from the \code{igraph} package), the
  graph corresponding to the passed penalty matrix. 
}
\usage{
getD1d(n)
getD1dSparse(n)
getD2d(dim1, dim2)
getD2dSparse(dim1, dim2)
getDg(graph)
getDgSparse(graph)
getDtf(n, ord)
getDtfSparse(n, ord) 
getDtfPos(n, ord, pos)
getDtfPosSparse(n, ord, pos)
getGraph(D)
}
\arguments{
  The arguments for the sparse variants are identical to those for the
  regular variants, which are described below.
  \item{n}{
    for \code{getD1d}, \code{getDtf}, and \code{getDtfPos}, 
    the number of points in the 1d sequence.
  }
  \item{dim1, dim2}{
    for \code{getD2d}, the number of rows and columns in the 2d grid,
    respectively. 
  }
  \item{graph}{
    for \code{getDg}, an \code{igraph} object from the \code{igraph}
    package, upon which the penalty matrix should be based (the penalty
    matrix is the oriented incidence matrix of the graph, with arbitrary
    orientations assigned to each edge).
  }
  \item{ord}{
    for \code{getDtf}, and \code{getDtfPos}, the order of the
    polynomial. E.g., ord=0 is the 1d fused lasso and ord=1 is linear 
    trend filtering. Hence the returned matrix is the discrete (ord+1)st
    derivative operator. 
  }
  \item{pos}{
    for \code{getDtfPos}, a numeric vector giving the positions of
    points in the 1d sequence. Must have length \code{n}. 
  }
  \item{D}{
    for \code{getGraph}, a fused lasso penalty matrix, the incidence
    matrix of an undirected graph, with arbitrary edge orientations.
  }
}
\value{
  All functions except \code{getGraph} return a penalty matrix, either
  in standard R matrix format or as a sparse matrix of class
  \code{dgCMatrix} via the \code{Matrix} package. The function
  \code{getGraph} returns an \code{igraph} object from the \code{igraph}
  package. 
}
\seealso{
  \code{\link{fusedlasso}}, \code{\link{trendfilter}}
}
\examples{
getD1d(9)
getDtfSparse(10,2)

graph = getGraph(getD2dSparse(4,4))
plot(graph)
}
\keyword{utilities}
