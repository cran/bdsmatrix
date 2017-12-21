\name{gchol.bdsmatrix-class}
\docType{class}
\alias{gchol.bdsmatrix-class}
\alias{\%*\%,gchol.bdsmatrix,matrix-method}
\alias{\%*\%,gchol.bdsmatrix,numeric-method}
\alias{\%*\%,matrix,gchol.bdsmatrix-method}
\alias{\%*\%,numeric,gchol.bdsmatrix-method}
\alias{[,gchol.bdsmatrix-method}
\alias{coerce,gchol.bdsmatrix,matrix-method}
\alias{diag,gchol.bdsmatrix-method}
\alias{dim,gchol.bdsmatrix-method}
\alias{show,gchol.bdsmatrix-method}

\title{Class "gchol.bdsmatrix"}
\description{Generalized cholesky decomposition of a \code{bdsmatrix} object,
A= LDL' where A is symmetric, L is lower triangular with 1 on the diagonal,
and D is diagonal.}
\section{Objects from the Class}{
These are created by the \code{gchol} function.
}
\section{Slots}{
	 \describe{
    \item{\code{blocksize}:}{Integer vector of block sizes}
    \item{\code{blocks}:}{Numeric vector containing the blocks}
    \item{\code{rmat}:}{Dense portion of the decomposition}
    \item{\code{rank}:}{The rank of A}
    \item{\code{Dim}:}{Integer vector of length 2 containing the dimension}
    \item{\code{Dimnames}:}{List of length 2 containing the dimnames}
  }
}
\section{Methods}{
  \describe{
    \item{\%*\%}{\code{signature(x = "gchol.bdsmatrix", y = "matrix")}: ... }
    \item{\%*\%}{\code{signature(x = "gchol.bdsmatrix", y = "numeric")}: ... }
    \item{\%*\%}{\code{signature(x = "matrix", y = "gchol.bdsmatrix")}: ... }
    \item{\%*\%}{\code{signature(x = "numeric", y = "gchol.bdsmatrix")}: ... }
    \item{[}{\code{signature(x = "gchol.bdsmatrix")}: ... }
    \item{coerce}{\code{signature(from = "gchol.bdsmatrix", to = "matrix")}: ... }
    \item{diag}{\code{signature(x = "gchol.bdsmatrix")}: ... }
    \item{dim}{\code{signature(x = "gchol.bdsmatrix")}: ... }
    \item{show}{\code{signature(object = "gchol.bdsmatrix")}: ... }
	 }
}
\author{Terry Therneau}
\note{
The Cholesky decompostion of a block diagonal symmetric matrix is also
block diagonal symmetric, so is stored in the same manner as a
\code{bdsmatrix} object}
\seealso{\code{\link{bdsmatrix}}, \code{\link{gchol}}}

\examples{
showClass("gchol.bdsmatrix")
}
\keyword{classes}
