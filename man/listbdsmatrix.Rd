\name{listbdsmatrix}
\Rdversion{1.1}
\alias{listbdsmatrix}
\title{List out a bdsmatrix as row/col/value triplets}
\description{
  This routine is the inverse of the bdsmatrix.ibd function
  found in the kinship library.
}
\usage{listbdsmatrix(x, id = TRUE, diag = FALSE)}
\arguments{
  \item{x}{a \code{bdsmatrix} object}
  \item{id}{if true, the dimnames of the object are used as the row and
    column identifiers in the output, if false integer row and column
    numbers are used}
  \item{diag}{include the diagonal elements in the output}
}
\details{
  The non-zero elements of the matrix are listed out as row-col-value
  triplets, one per line, in a data frame.  Since the matrix is known to
  be symmetric, only elements with row >= col are listed.  When familial
  correlation data is represented in a bdsmatrix, e.g. kinship or
  identity-by-descent information, the diagonal is a known value and can
  be omitted from the listing.
  Genetic software often produces matrices in the list form; this
  routine is the inverse of the bdsmatrix.ibd routine, found in the
  kinship library, which converts list form to bdsmatrix form.
}
\value{a data frame with variables \code{row}, \code{col}, and
  \code{value}.
  }
\author{Terry Therneau}
\seealso{\code{\link{bdsmatrix}}}

