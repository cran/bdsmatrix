#
# A temporary hack for K Scurrah, to replace some elements in a bdsmatrix
#
#  This should be a more general [<- method
#
# Given a bdsmatrix "xt", and a set of locations, replace those elements.
#   if irow= 1,2,3, icol= 4,6 2, new=10,11,12, then set
#     x[1,4]=10, x[2,6]=11 and x[3,2]=12
#
 
subhack <- function(x, irow, icol, new) {
    dd <- dim(x)
    n <- length(irow)
    if (length(icol) !=n || length(new) !=n) 
        stop("Mismatched lengths")

    if (any(irow <1 | irow > dd[1] | irow !=floor(irow)))
        stop("Invalid irow index")

    if (any(icol <1 | icol > dd[1] | icol !=floor(icol)))
        stop("Invalid icolumn index")
    
    if (class(x) != 'bdsmatrix') stop("Invalid matrix")

    # This routine gives the row and col number of each element found in
    #   the sparse packing.
    indx <- .C('bdsmatrix_index2',
               as.integer(length(x@blocksize)),
               as.integer(x@blocksize),
               rows= integer(length(x@blocks)),
               cols= integer(length(x@blocks)))

    temp <- x@blocks

    # Only the lower diagonal is stored, so change e.g. x[1,4] to x[4,1]
    ir2 <- pmax(irow, icol)  
    icol<- pmin(irow, icol)

    for (i in 1:n) {
        location <- which(indx$rows==ir2[i] & indx$cols==icol[i])
        if (length(location) ==0) 
            stop("Tried to replace an element in the sparse portion")
        
        temp[location] <- new[i]
        }
    x@blocks <- temp
    x
    }



               
