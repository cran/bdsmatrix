\name{bdsmatrix-class}
\docType{class}
\alias{bdsmatrix-class}
\alias{\%*\%,matrix,bdsmatrix-method}
\alias{\%*\%,numeric,bdsmatrix-method}
\alias{\%*\%,bdsmatrix,matrix-method}
\alias{\%*\%,bdsmatrix,numeric-method}
\alias{Math2,bdsmatrix-method}
\alias{Math,bdsmatrix-method}
\alias{Ops,bdsmatrix,numeric-method}
\alias{Ops,bdsmatrix,bdsmatrix-method}
\alias{Ops,bdsmatrix,matrix-method}
\alias{Ops,numeric,bdsmatrix-method}
\alias{Ops,matrix,bdsmatrix-method}
\alias{[,bdsmatrix-method}
\alias{all,bdsmatrix-method}
\alias{any,bdsmatrix-method}
\alias{coerce,bdsmatrix,matrix-method}
\alias{coerce,bdsmatrix,vector-method}
\alias{diag,bdsmatrix-method}
\alias{diag<-,bdsmatrix-method}
\alias{dim,bdsmatrix-method}
\alias{dimnames,bdsmatrix-method}
\alias{dimnames<-,bdsmatrix-method}
\alias{gchol,bdsmatrix-method}
\alias{max,bdsmatrix-method}
\alias{min,bdsmatrix-method}
\alias{prod,bdsmatrix-method}
\alias{range,bdsmatrix-method}
\alias{show,bdsmatrix-method}
\alias{sum,bdsmatrix-method}

\title{Class "bdsmatrix"}
\description{Representation for a Block Diagonal Sparse matrix}
\section{Objects from the Class}{
Objects of this class are usually created using the \code{bdsmatrix},
\code{bdsI} or \code{bdsBlock} functions.
The result is a symmetrix matrix whose upper left portion is block-diagonal,
with an optional border on the right and bottom that is dense.
The matrices were originally created to represent familial correlation
structures, which have a block for each family but no connection between
families.
}
\section{Slots}{
	 \describe{
    \item{\code{blocksize}:}{An integer vector containing the
sizes of the diagonal blocks}
    \item{\code{blocks}:}{A numeric vector containing the contents of the
block portion.  Only the lower triangle of each block is stored.}
    \item{\code{rmat}:}{An optional numeric matrix containing the dense 
portion}
    \item{\code{offdiag}:}{A single numeric element, default zero, which is
the value for elements off the block-diagonal }
    \item{\code{Dim}:}{The dimension of the matrix, an integer vector of length 2}
    \item{\code{Dimnames}:}{The dimnames of the matrix, a list with 2 elements}
  }
}
\section{Methods}{
  \describe{
    \item{\%*\%}{\code{signature(x = "matrix", y = "bdsmatrix")}: the result
will be an ordinary matrix}
    \item{\%*\%}{\code{signature(x = "numeric", y = "bdsmatrix")}: the result
will be a vector}
    \item{\%*\%}{\code{signature(x = "bdsmatrix", y = "matrix")}:  the result
will be an ordinary matrix}
    \item{\%*\%}{\code{signature(x = "bdsmatrix", y = "numeric")}:  the result
will be a vector}
    \item{Math2}{\code{signature(x = "bdsmatrix")}:  }
    \item{Math}{\code{signature(x = "bdsmatrix")}:  }
    \item{Ops}{\code{signature(e1 = "bdsmatrix", e2 = "numeric")}: }
    \item{Ops}{\code{signature(e1 = "bdsmatrix", e2 = "bdsmatrix")}: }
    \item{Ops}{\code{signature(e1 = "bdsmatrix", e2 = "matrix")}: }
    \item{Ops}{\code{signature(e1 = "numeric", e2 = "bdsmatrix")}:  }
    \item{Ops}{\code{signature(e1 = "matrix", e2 = "bdsmatrix")}: }
    \item{[}{\code{signature(x = "bdsmatrix")}: if the subscripts are a
set of increasing integers, and the row and column subscripts are identical,
then the result is aslo a bdsmatrix.  This is useful for example to create
the kinship matrix for all females from an overall kinship matrix.  If the
subscripts do not match, then an ordinary matrix is created}
    \item{all}{\code{signature(x = "bdsmatrix")}: ... }
    \item{any}{\code{signature(x = "bdsmatrix")}: ... }
    \item{coerce}{\code{signature(from = "bdsmatrix", to = "matrix")}: ... }
    \item{coerce}{\code{signature(from = "bdsmatrix", to = "vector")}: ... }
    \item{diag}{\code{signature(x = "bdsmatrix")}: retrieve the diagonal of
the matrix}
    \item{diag<-}{\code{signature(x = "bdsmatrix")}: set the diagonal of the
matrix to a given value}
    \item{dim}{\code{signature(x = "bdsmatrix")}: dimension of the matrix}
    \item{dimnames}{\code{signature(x = "bdsmatrix")}: dimnames of the
matrix}
    \item{dimnames<-}{\code{signature(x = "bdsmatrix")}: set the dimnames of
the matrix}
    \item{gchol}{\code{signature(x = "bdsmatrix")}: generalized cholesky 
decomposition of the matrix }
    \item{max}{\code{signature(x = "bdsmatrix")}: maximum of the matrix}
    \item{min}{\code{signature(x = "bdsmatrix")}: minimum of the matrix}
    \item{prod}{\code{signature(x = "bdsmatrix")}: }
    \item{range}{\code{signature(x = "bdsmatrix")}: }
    \item{show}{\code{signature(object = "bdsmatrix")}: print out the matrix}
    \item{sum}{\code{signature(x = "bdsmatrix")}: }
	 }
}
\note{
Many of the actions above will result in conversion to an ordinary matrix
object, including \code{print}, addition to an ordinary matrix, etc.
This can easily create objects that are too large for system memory.
By default the value of options('bdsmatrixsize') is consulted first, and
if the resulting object would be have a length greater than this option
the conversion an error is generated and conversion is not attempted.
The default value for the option is 1000.
}
\author{Terry Therneau}
\seealso{\code{\link{gchol}}
}
\examples{
showClass("bdsmatrix")
}
\keyword{classes}
