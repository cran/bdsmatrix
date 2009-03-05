\name{as.matrix.bdsmatrix}
\alias{as.matrix.bdsmatrix}
\title{Convert a bdsmatrix to a ordinary (dense) matrix}
\description{Method to convert from a Block Diagonal Sparse (bdsmatrix)
matrix representation to an ordinary one}

\usage{
\method{as.matrix}{bdsmatrix}(x, ...)}
\arguments{
	\item{x}{a bdsmatrix object}
	\item{...}{other arguments are ignored (necessary to match the 
\code{as.matrix} template)}
}

\value{ a matrix}
\details{
  Note that the conversion of a large bdsmatrix can easily exceed memory.
}
\seealso{bdsmatrix
}
\keyword{array}
% docclass is function
% Converted by Sd2Rd version 43267.
