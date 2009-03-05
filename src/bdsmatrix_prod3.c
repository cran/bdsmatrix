/*
**  Product of a gchol.bdsmatrix object and a regular one
**
** nr           dimension of the  bdsmatrix (it is square)  L
** nb		number of blocks       for the bdsmatrix
** bsize	the block sizes             ""
** bmat		the vector of blocks	    ""
** rmat         right hand matrix           ""
** rhs          "right hand side?" if=0, then compute yL, if 1 then Ly
** ny	        number of rows of y (if rhs=0), or number of cols (if rhs=1)
** y	        the y vector or matrix,
** temp		a scratch vector of doubles of length nr, only needed if
**                rhs=1
**
** The S code has already verified that the dimensions of y match those of L
*/
#include "bdsS.h"
#include "bdsmatrix.h"
#include <stdio.h>
double sqrt(double);

void bdsmatrix_prod3(Sint *nr,     Sint *nb,	  Sint *bsize,     
		     double *bmat, double *rmat,  Sint *rhs,
		     Sint *ny2,	   double *y,     double *temp) {

    int nblock;
    int nrow, ny ;
    int brow, rrow;
    int i,j, k, col, yrow;
    int itemp;
    int nk;
    int icol;
    int blocksize, offset, irow, n, block;
    double sum, scale;
    double *x, *rx;
    
    nblock = *nb;
    nrow   = *nr;
    ny     = *ny2;
    
    brow =0;	/* number of rows in the block diagonal portion */
    for (i=0; i<nblock; i++) brow += bsize[i];
    rrow = nrow -brow;   /* this will be 0 if there is no rmat */

    if (*rhs==1) {
	/*
	**  I think that this will be the much less frequent case
	** Do one colum of the matrix at a time
	*/
	for (col=0; col<ny  ; col++) {
	    offset = col*nrow;
	    /*
	    ** The indexing below is a little opaque, but fast when done right
	    ** irow = the current row of bmat being processed
	    ** n    = the index of the [irow, irow] element
	    ** k    = irow of the start of this block
	    ** nk   = the index of the [k,k] element of bmat
	    ** itemp = the indices of the elements representing bmat[irow,]
	    **    Say the first block is 4 by 4
	    **     	irow=0, itemp= 0
	    **    	irow=1, itemp= 1 4
	    **   	irow=2, itemp= 2 5 7
	    **	irow=3, itemp= 3 6 8 9
	    **
	    **  The product of bmat[irow,] and y[,col] is placed in temp[irow]
	    */
	    irow =0;
	    n=0;
	    for (block=0; block < nblock; block++) {
		blocksize = bsize[block];
		k = irow;
		nk =n;
		for (i=0; i<blocksize; i++) {  /* march down the rows */
		    sum = 0;
		    y[offset+irow] *= sqrt(bmat[n]); /* y = sqrt(D) y */ 
		    sum = y[offset+irow];              /* L[i,i] * y[i] */
		    itemp = nk +i;
		    for (j=0; j<i; j++) {
			sum += bmat[itemp] * y[offset+j + k];
			itemp += blocksize - (j+1);
		        }
		    temp[irow] = sum;
		    irow++;
		    n += blocksize -i;
		    }
	        }
	
	    /* 
	    ** Add in the rmat part, if present 
	    **  n is now the index of the first element of the row of rmat
	    */
	    n=0;
	    for (i=0; i<rrow; i++) {
		y[offset + irow] *= sqrt(rmat[n+irow]);
		sum = y[offset + irow];
		for (j=0; j< irow; j++) 
		    sum += rmat[n+j] * y[offset+j];
		temp[irow] = sum;
		irow++;
		n += nrow;
	        }
	    /*
	    **  Copy the temp vector back over the top of y
	    */
	    for (i=0; i<nrow; i++) y[offset+i] = temp[i];
	    offset += nrow;
	    }
	}

    else {
	/*
	** The y matrix is on the left (the more common case)
	**   The indexing is much easier in this case, as well.
	** 
	** Each row of the y matrix (left hand side) is done separately
	**  For that row, we replace the first element of the row by the
	**  new result, then the second, etc.  The current element of focus,
	**  which is also the current column of L, is "icol".
	** The multiply/add is done in exactly the same order you would do
	**  it with non-sparse matrices on paper, with "sum" as a the 
	**  temporary sum, 
	*/
	for (yrow=0; yrow < ny; yrow++) {  /* separate compute for each */
	    x = bmat ;   /* this walks throught the bdsmatrix */
	    icol =0;     /* we will create the [yrow, icol] element of answer*/
	    for (block=0; block < nblock; block++) {
		for (j= bsize[block]; j>0; j--) {
		    offset = yrow + icol*ny;  /* starting with y[offset] */
		    scale = sqrt(*x);
		    sum = y[offset] *scale;  /* mult by implicit 1 on diag */
		    x++;
		    offset += ny;
		    /* First the sparse rows, beyond diag, for this col of L */
		    for (k=1; k<j; k++) {
			sum += *x *  scale * y[offset];
			x++;
			offset += ny;
			}
		    /*
		    ** Now the dense rows
		    */
		    offset = yrow + brow*ny;  /* bounce ahead */
		    rx = rmat + icol;
		    for (k=0; k<rrow; k++) {
			sum += *rx * scale * y[offset];
			offset += ny;
			rx += nrow;
			}

		    /* save the result */
		    y[yrow + icol*ny] = sum;  /* write the result over y */
		    icol++;
		    }
		}
	    /*
	    ** Now, the product that involves only "pure dense" rows
	    **  of the bdsmatrix
	    */
	    for (i=0; i<rrow; i++) {
		offset = yrow + icol*ny;
		rx = rmat + icol + i*nrow;
		scale = sqrt(*rx);
		sum = scale * y[offset];

		for (j=i+1; j<rrow; j++) {
		    rx += nrow;
		    offset += ny;
		    sum += *rx * scale * y[offset];
                    }
		y[yrow + icol*ny] = sum;  /* write the result over y */
		icol++;
		}
	    } /* on to the next row of y */
	}
    }

