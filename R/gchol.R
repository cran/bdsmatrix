#
# Code for the generalized cholesky  A = LDL', where L is lower triangular
#   with 1's on the diagonal, and D is diagonal.
# The decompostions exists for any square symmetric matrix.  
# If A is positive definite, then all elements of D will be positve.
# If A is not full rank, then 0's on the diagonal of D signal the redundant 
#   columns.  
# Note that gchol is both a class (setClass) and a generic function.
#
setClass('gchol', 
	 representation(.Data= 'numeric',
			Dim = 'integer',
			Dimnames = 'list',
			rank = 'integer'))
			
setGeneric('gchol', function(x, tolerance=1e-10) standardGeneric('gchol'),
           useAsDefault=FALSE)

as.matrix.gchol <- function(x, ones=TRUE, ...) {
    temp <- matrix(x@.Data, x@Dim[1], dimnames=x@Dimnames, byrow=TRUE)
    if (ones) diag(temp) <- 1
    temp
    }

setAs('gchol', 'matrix', function(from) as.matrix.gchol(from))

setMethod('gchol', signature(x='matrix'),
    function(x,  tolerance) {
	d <- dim(x)
	if (d[1] != d[2]) 
		stop("Cholesky decomposition requires a square matrix")
	if (!is.logical(all.equal(as.vector(x), as.vector(t(x)))))
		stop("Cholesky decomposition requires a symmetric matrix")
	temp <- .C("gchol", as.integer(d[1]),
		   x =   as.double(x),
		   rank= as.double(tolerance))
      
        dnames <- dimnames(x)
        if (is.null(dnames)) dnames <- list(NULL, NULL)

	new('gchol', .Data= temp$x  , Dim=d, 
	    Dimnames= dnames, rank=as.integer(temp$rank))
	})

setMethod('diag', signature(x='gchol'),
    function(x, nrow, ncol) {
	d <- x@Dim[1]
	x@.Data[ seq(1, length=d, by=d+1)]
	})


setMethod('show', 'gchol', function(object) show(as.matrix(object, F)))
setMethod('dim', 'gchol',  function(x) x@Dim)
setMethod('dimnames', 'gchol', function(x) x@Dimnames)
 	

