\name{backsolve}
\alias{backsolve}
\alias{backsolve.default}
\alias{backsolve.gchol}
\alias{backsolve.gchol.bdsmatrix}
\title{Solve an Upper or Lower Triangular System}
\description{
  Solves a system of linear equations where the coefficient matrix is
  upper (or \sQuote{right}, \sQuote{R}) or lower (\sQuote{left},
  \sQuote{L}) triangular.\cr 

  \code{x <- backsolve(R, b)} solves \eqn{R x = b}.
}
\usage{
  backsolve(r, x, ...)
  \method{backsolve}{default}(r, x, k=ncol(r), upper.tri=TRUE,
         transpose=FALSE, \dots)
  \method{backsolve}{gchol}(r, x, k=ncol(r), upper.tri=TRUE, \dots)
  \method{backsolve}{gchol.bdsmatrix}(r, x, k=ncol(r), upper.tri=TRUE, \dots)
}
\arguments{
  \item{r}{a matrix or matrix-like object}
  \item{x}{a vector or a matrix whose columns give the right-hand sides for
    the equations.}
  \item{k}{The number of columns of \code{r} and rows of \code{x} to use.}
  \item{upper.tri}{logical; if \code{TRUE} (default), the \emph{upper}
    \emph{tri}angular part of \code{r} is used.  Otherwise, the lower
    one.}
  \item{transpose}{logical; if \code{TRUE}, solve \eqn{r' * y = x} for
    \eqn{y}, i.e., \code{t(r) \%*\% y == x}.}
  \item{\dots}{further arguments passed to other methods}
}
\value{
  The solution of the triangular system.  The result will be a vector if
  \code{x} is a vector and a matrix if \code{x} is a matrix.

  Note that \code{forwardsolve(L, b)} is just a wrapper for
  \code{backsolve(L, b, upper.tri=FALSE)}.
}
\note{
  The \code{bdsmatrix} package promotes the base R function to a
  generic.
  To see the full documentation for the default method, view \code{backsolve}
  from the \code{base} package.
}
\seealso{
\code{\link{forwardsolve}}
}
\keyword{ array }
\keyword{ algebra }
