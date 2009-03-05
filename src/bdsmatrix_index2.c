/*
** These are support for the "rows" and "cols" functions for 
**  block diagonal matrices.  Returns the list of row numbers and
**  the list of column numbers
*/
#include "bdsS.h"
#include "bdsmatrix.h"
void bdsmatrix_index2(Sint *nblock, Sint *bsize, 
		      Sint *rows,   Sint *cols) { 

    int i, j;
    int blocksize;
    int block;           /* block currently being processed */
    int n;               /* indexes the return list */
    int irow;            /* row number of the upper corner of the block */

    n =0;  
    irow=0;
    for (block=0; block < *nblock; block++) {
	blocksize = bsize[block];
	for (i=0; i<blocksize; i++) {
	    for (j=i; j<blocksize; j++) {
		rows[n] = 1 + irow +j -i;
		cols[n] = 1 + irow;
		n++;
		}
	    irow++;
	    }
	}
    }
