/*
**  Product of a gchol.bdsmatrix object and a vector
**   Variant of bdsmatrix_prod3, for internal calls by C
**   instead of from S.  The multiplication vector is (b, beta),
**   where beta is known to be zero, so only the frailties b
**   are passed in.
**
** nrow         total number of rows in the bdsmatrix
** nblock	number of blocks       for the bdsmatrix
** bsize	the block sizes             ""
** bmat		the vector of blocks	    ""
** rmat         right hand matrix           ""
** nfrail       number of rows that we are using (length of y)
** y	        the left hand matrix, which will be overwritten
**
*/
#include "bdsmatrix.h"
double sqrt(double);

void bdsmatrix_prod4(int nrow,  int nblock,  int *bsize, 
		    double *bmat, double *rmat,    
		    int nfrail, double *y) {

    int brow, rrow;
    int i,j, k;
    int block;
    int icol, offset;
    double *x, *rx;
    double sum, scale;
    
    brow =0;	/* number of rows in the block diagonal portion */
    for (i=0; i<nblock; i++) brow += bsize[i];
    rrow = nfrail - brow; /* number or rows of r that we care about */

    x = bmat ;   /* this walks throught the bdsmatrix */
    icol =0;     /* we will create the [icol] element of answer*/
    for (block=0; block < nblock; block++) {
	for (j= bsize[block]; j>0; j--) {
	    offset = icol;  /* starting with y[offset] */
	    scale = sqrt(*x);
	    sum = y[icol] *scale;  /* mult by implicit 1 on diag */
	    x++;
	    offset++;
	    /* First the sparse rows, beyond diag, for this col of L */
	    for (k=1; k<j; k++) {
		sum += *x *  scale * y[offset];
		x++;
		offset++;
	        }
	    /*
	    ** Now the dense rows
	    */
	    offset = brow;  /* bounce ahead */
	    rx = rmat + icol;
	    for (k=0; k<rrow; k++) {
		sum += *rx * scale * y[offset];
		offset++;
		rx += nrow;
	        }

	    /* save the result */
	    y[icol] = sum;  /* write the result over y */
	    icol++;
	    }
	}
    /*
    ** Now, the product that involves only "pure dense" rows
    **  of the bdsmatrix
    */	
    for (i=0; i<rrow; i++) {
	offset = icol;
	rx = rmat + icol + i*nrow;
	scale = sqrt(*rx);
	sum = scale * y[offset];

	for (j=i+1; j<rrow; j++) {
	    rx += nrow;
	    offset++;
	    sum += *rx * scale * y[offset];
	    }
	y[icol] = sum;  /* write the result over y */
	icol++;
	}
    }
