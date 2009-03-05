\name{gchol-class}
\docType{class}
\alias{gchol-class}
\alias{coerce,gchol,matrix-method}
\alias{diag,gchol-method}
\alias{dim,gchol-method}
\alias{dimnames,gchol-method}
\alias{show,gchol-method}
\alias{gchol,matrix-method}

\title{Class "gchol"}
\description{The result of a generalized Cholesky decomposition A=LDL'
where A is a symmetric matrix, L is lower triangular with 1s on the
diagonal, and D is a diagonal matrix.}
\section{Objects from the Class}{
These objects are created by the \code{gchol} function.
}
\section{Slots}{
	 \describe{
    \item{\code{.Data}:}{A numeric vector containing the results of the
decompostion}
    \item{\code{Dim}:}{An integer vector of length 2, the dimension of 
the matrix}
    \item{\code{Dimnames}:}{A list of length 2 containing the dimnames.  These
default to the dimnames of the matrix A}
    \item{\code{rank}:}{The rank of the matrix}
  }
}
\section{Methods}{
  \describe{
    \item{coerce}{\code{signature(from = "gchol", to = "matrix")}: Use of
the \code{as.matrix} function will return L}
    \item{diag}{\code{signature(x = "gchol")}: Use of the \code{diag} function
will return D}
    \item{dim}{\code{signature(x = "gchol")}: returns the dimension of the 
matrix }
    \item{dimnames}{\code{signature(x = "gchol")}: returns the dimnames}
    \item{show}{\code{signature(object = "gchol")}: By default a triangular
      matrix is printed showing D on the diagonal and L off the diagonal
    }
    \item{gchol}{\code{signature(x= "matrix")}: create a generalized
      Cholesky decompostion of the matrix}
	 }
}
\author{Terry Therneau}
\note{The primary advantages of the genearlized decomposition, as compared to
the standard \code{chol function}, has to do with redundant columns 
and generalized inverses (g-inverse).
The lower triangular matrix L is always of full rank.  The diagonal matrix
D has a 0 element at position j if and only if the jth column of A is
linearly dependent on columns 1 to j-1 preceding it. 
The g-inverse of A involves the inverse of L and a g-inverse of D.
The g-inverse of D retains the zeros and inverts non-zero elements
of D.   
All this makes manipulation of the result particularly convenient.}
\seealso{\code{\link{gchol}}}
\examples{
showClass("gchol")
}
\keyword{classes}
