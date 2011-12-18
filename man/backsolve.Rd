\name{backsolve}
\alias{backsolve-methods}
\title{Solve an Upper or Lower Triangular System}
\alias{backsolve}
\alias{backsolve,gchol-method}
\alias{backsolve,gchol.bdsmatrix-method}
\description{
  Solves a system of linear equations where the coefficient matrix is
  upper (or \sQuote{right}, \sQuote{R}) or lower (\sQuote{left},
  \sQuote{L}) triangular.\cr 

  \code{x <- backsolve(R, b)} solves \eqn{R x = b}.
}
\usage{
  backsolve(r, \dots)
  \S4method{backsolve}{gchol}(r, x, k=ncol(r), upper.tri=TRUE, \dots)
  \S4method{backsolve}{gchol.bdsmatrix}(r, x, k=ncol(r), upper.tri=TRUE, \dots)
}
\arguments{
  \item{r}{a matrix or matrix-like object}
  \item{x}{a vector or a matrix whose columns give the right-hand sides for
    the equations.}
  \item{k}{The number of columns of \code{r} and rows of \code{x} to use.}
  \item{upper.tri}{logical; if \code{TRUE} (default), the \emph{upper}
    \emph{tri}angular part of \code{r} is used.  Otherwise, the lower
    one.}
  \item{\dots}{further arguments passed to other methods}
}
\value{
  The solution of the triangular system.  The result will be a vector if
  \code{x} is a vector and a matrix if \code{x} is a matrix.

  Note that \code{forwardsolve(L, b)} is just a wrapper for
  \code{backsolve(L, b, upper.tri=FALSE)}.
}
\section{Methods}{
  Use \code{\link{showMethods}(backsolve)} to see all the defined methods;
    the two created by the bdsmatrix library are described here:
    \describe{
      \item{bdsmatrix}{\code{signature=(r= "gchol")} for a generalized
	cholesky decomposition}
      \item{bdsmatrix}{\code{signature=(r= "gchol.bdsmatrix")} for the
	generalize cholesky decomposition of a bdsmatrix object}
    }
  }
\details{
  The generalized Cholesky decompostion of a symmetric matrix A is
  \eqn{A = LDL'}{A= LD t(L)} where D is diagonal, L is lower triangular,
  and \eqn{L'}{t(L)} is the transpose of L.
  These functions solve either \eqn{L\sqrt{D} x =b}{L sqrt(D) x=b}
  (when \code{upper.tri=FALSE}) or \eqn{\sqrt{D}L' x=b}{sqrt(D) t(L)
    x=b}.
  }
\note{
  The \code{bdsmatrix} package promotes the base R \code{backsolve}
  function to a
  generic.
  To see the full documentation for the default method, view \code{backsolve}
  from the \code{base} package.
}
\seealso{
\code{\link{forwardsolve}}, \code{\link{gchol}}
}
\keyword{ array }
\keyword{ algebra }
