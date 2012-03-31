/*
** Cholesky operations for block-diagonal-symmetric matrices
**
** As a further speedup, we allow the solve and inv routines to be
**  called with either a pre-done cholesky decompostion or the raw matrix.
** 
*/
#include "bdsS.h"
#include "bdsmatrix.h"

/*
** Create a Cholesky decompostion
*/
void gchol_bds(Sint   *nb,     Sint   *bs2,  Sint *n2,
	       double *dmat,   double *rmat, double toler[]) {
    int i,j;

    int *bsize,
	bsum,
	n, 
	nblock;
    double **mat;

    nblock = *nb;
    n = *n2;
    /* 
    ** copy over arguments from long to int form, if needed 
    ** (it will be needed in Splus, not needed currently in R)
    ** recant -- even when sizeof(long) == sizeof(int), crashes linux/R
    */
    /* if (sizeof(Sint) != sizeof(int)) {  */
    bsize = (int *) ALLOC(nblock, sizeof(int));
    bsum =0;
    for (i=0; i<nblock; i++) {
	bsize[i] = bs2[i];
	bsum += bsize[i];
	}

    /* create indices for the right-hand side matrix, if it is present */
    if (n > bsum) {
	mat = dmatrix(rmat, n, n-bsum);
	}
    else {
	mat = (double **) 0;
	/* in this case cholesky4 will never touch the "mat" argument 
           doing nothing here certainly LOOKS like a bug, however  */
    }

    i = cholesky4(mat, n, nblock, bsize, dmat, *toler);
    *toler = i; 

    /* zero out the upper triangle */
    for (i=0; i<n-bsum; i++){  /*columns of R */
	for (j= 1+i+bsum; j<n; j++) mat[i][j] =0;  /* zero below diagonal */
	}
    }

/*
** The cholesky decomp is A=LDL'
**    flag= 0: input is A (original matrix)  
**       or 1: input is LD
**    plus 
**          0: return inverse of A
**       or 2: return inverse of L
**
** The "inverse" of A is not a real inverse.  The inverse of L is sparse,
**   but the inverse of A is not.  Since A is stored in sparsely, and we
**   overwrite it with it's inverse, only some of the elements of the inverse
**   are returned.
** By default solve.bdsmatrix uses flag=0+2, and then multiplies things out
**   itself to get the real inverse.  solve.gchol.bdsmatrix uses 1+0 (full=T)
**   or 1+2 (full=F).  The 0+0 case is used by solve.bdsmatrix with the full=F
**   argument -- for degrees of freedom, coxme only needs the block diagonal
**   elements returned.
*/
void gchol_bdsinv(Sint   *nb,     Sint   *bs2,  Sint *n2,
		  double *dmat,   double *rmat, double *toler,
                  Sint   *flag) {
    int i,j;

    int *bsize,
	bsum,
	n, 
	nblock;
    double **mat;

    /* copy over arguments from Sint to int form */
    nblock = *nb;
    n = *n2;
    bsize = (int *) ALLOC(nblock, sizeof(int));
    bsum =0;
    for (i=0; i<nblock; i++) {
	bsize[i] = bs2[i];
	bsum += bsize[i];            /* also total up the block sizes */
	}

    /* create indices for the right-hand side matrix, if it is present */
    if (n > bsum) {
	mat = dmatrix(rmat, n, n-bsum);
	}

    if (*flag==0 || *flag==2) {
	i = cholesky4(mat, n, nblock, bsize, dmat, *toler);
	*toler = i;
	for (i=0; i<n-bsum; i++) {  /*columns of R */
	    for (j= 1+i+bsum; j<n; j++) mat[i][j] =0; /*zero out */
	    }
	}
    if (*flag>=2) chinv4(mat, n, nblock, bsize, dmat, 0);
    else          chinv4(mat, n, nblock, bsize, dmat, 1);
    }

/*
** Solve Ab = y for an input vector y.  y is overwritten with b.
**  The decompostion is A=LDL'
**    flag= 0: input is A (original matrix)  
**       or 1: input is LD
**    plus 
**          0: return solution to Ab=y
**       or 2: return solution to sqrt(D)L'b =y
**
*/
void gchol_bdssolve(Sint   *nb,     Sint   *bs2,  Sint *n2,
		    double *blocks, double *rmat, double *toler,
		    double *y,      Sint   *flag) {
    int i,j;

    int *bsize,
	bsum,
	n, 
	nblock;
    double **mat;

    /* copy over arguments from Sint to int form */
    nblock = *nb;
    n = *n2;
    bsize = (int *) ALLOC(nblock, sizeof(int));
    bsum =0;
    for (i=0; i<nblock; i++) {
	bsize[i] = bs2[i];
	bsum += bsize[i];
	}

    /* create indices for the right-hand side matrix, if it is present */
    if (n > bsum) {
	mat = dmatrix(rmat, n, n-bsum);
	}

    if (*flag==0 || *flag==2) {
	i = cholesky4(mat, n, nblock, bsize, blocks, *toler);
	for (i=0; i<n-bsum; i++) {  /*columns of R */
	    for (j= 1+i+bsum; j<n; j++) mat[i][j] =0; /*zero out */
	    }
	}

    if (*flag >1) chsolve4(mat, n, nblock, bsize, blocks, y, 1);
    else          chsolve4(mat, n, nblock, bsize, blocks, y, 0);
    }

    
	
    
