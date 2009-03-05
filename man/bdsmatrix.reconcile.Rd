\name{bdsmatrix.reconcile}
\Rdversion{1.1}
\alias{bdsmatrix.reconcile}
\title{
Ensure alignment of two bdsmatrix objects
}
\description{
  This function is used by coxme.  When a random effect is expressed
  as a sum of variance terms (matrices), it is important that all of
  them have the same row/column order and the same block structure.
  This does so, while retaining as much sparsity in the result as
  possible.
}
\usage{
bdsmatrix.reconcile(varlist, group)
}
\arguments{
  \item{varlist}{
a list, each element of which is a matrix or bdsmatrix object
}
  \item{group}{
a vector of dimnames, the target match for matrice's dimnames
}
}
\value{
a varlist, whose individual elements may have had row/column  rearrangment.
}
\author{Terry Therneau}

\seealso{
\code{\link{bdsmatrix}}
}
\keyword{array}
