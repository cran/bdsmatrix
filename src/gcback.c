/*
** The backsolve function for a gchol matrix
*/
#include "bdsS.h"
#include "bdsmatrix.h"

SEXP gcback(SEXP sr,   SEXP  sx,    SEXP supper, SEXP sk) {
    int i;
    double *r;
    double **rmat;

    int nr;   /*number of rows in x, or rows/cols in r */
    int nc;   /*number of columns in x */
    int upper;
    int flag;
    int k;    /*number of colums to use in the solution */
    /*
    ** structure for the returned object
    */
    SEXP sy;
    double *y;

    sy = PROTECT(duplicate(sx)); /* this cause row/col names to duplicate */
    y = REAL(sy);
    r = REAL(sr);

    nr = nrows(sx);
    nc=  ncols(sx);
    k =  asInteger(sk);
    upper = asLogical(supper);
    flag = 1+ upper;   /* for chsolve5, 2=upper and 1=lower */
    rmat = dmatrix(r, nr, nr);

    for (i=0; i<nc; i++) {
	    chsolve5(rmat, k, y, flag);
	    y += nr;
	    }
    
    UNPROTECT(1);
    return(sy);
}
	
/*
** Same for a gchol.bsdmatrix
*/
SEXP gcback2(SEXP sblocksize,   SEXP  sblocks,    SEXP srmat,
	      SEXP sx,           SEXP  supper) {
    int i;
    double **rmat;
    double *blocks;
    int *blocksize;

    int nr;   /*number of rows in x, or rows/cols in r */
    int nc;   /*number of columns in x */
    int upper;
    int flag;
    /*
    ** structure for the returned object
    */
    SEXP sy;
    double *y;
    
    blocksize = INTEGER(sblocksize);
    blocks = REAL(sblocks);
    if (ncols(srmat) >0) 
	rmat = dmatrix(REAL(srmat), ncols(srmat), nrows(srmat));
    else rmat=NULL;

    sy = PROTECT(duplicate(sx)); /* this cause row/col names to duplicate */
    y = REAL(sy);
    
    nr = nrows(sx);
    nc=  ncols(sx);
    upper = asLogical(supper);
    flag = 1+upper;   /* for chsolve4, 2=lower and 1=upper */
    rmat = dmatrix(REAL(srmat), nr, nr);
    
    for (i=0; i<nc; i++) {
	    chsolve4(rmat, nr, LENGTH(sblocksize), blocksize,
		     blocks, y, flag);  
	    y += nr;
	    }
    
    UNPROTECT(1);
    return(sy);
}
	
		
	    
	
