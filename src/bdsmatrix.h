/*
** This contains the prototype calls for all the .c functions that
**    are called by another C function, or by R
**  It stops errors due to having things declared differently
**    in different routines, and serves as input to R_init_bdsmatrix
*/
void bdsmatrix_index1(Sint *nblock, Sint *bsize, Sint *flag,
		     Sint *nrow,   Sint *rows,  Sint *indexa,
		     Sint *indexb, Sint *indexc);

void bdsmatrix_index2(Sint *nblock, Sint *bsize, 
		      Sint *rows,   Sint *cols);
void bdsmatrix_index3(Sint *nblock, Sint *bsize, Sint *index);

void bdsmatrix_prod(Sint *nb,	  Sint *bsize,     Sint *ydim,
		    double *bmat, double *rmat,    double *offdiag,
		    double *temp, Sint *itemp,     double *y);

void bdsmatrix_prod2(int nblock,     int *bsize,     int nrow,
		     double *bmat,   double *rmat,  
		     double *y,      double *result, int *itemp);

void bdsmatrix_prod3(Sint *nr,     Sint *nb,	  Sint *bsize,     
		     double *bmat, double *rmat,  Sint *rhs,
		     Sint *ny2,	   double *y,     double *temp);

void bdsmatrix_prod4(int nrow,  int nblock,  int *bsize, 
		    double *bmat, double *rmat,    
		    int nfrail, double *y);

void chinv4(double **matrix, 	int n, 	 	int nblock, 	int *bsize, 
	    double *bd,	        int flag) ;

void chinv5(double **matrix ,	int n,		int flag);

int cholesky4(double **matrix,	int n,		int nblock, 	int *bsize,
	      double *bd, 	double toler) ;

int cholesky5(double **matrix, 	int n, 		double toler);

void chsolve4(double **matrix, 	int n, 		int nblock, 	int *bsize,
	      double *bd, 	double *y, 	int flag);

void chsolve5(double **matrix, 	int n, 		double *y, 	int flag);

double **dmatrix(double *array, int ncol, 	int nrow);

void gchol(Sint *n2, double *matrix, double *toler);
void gchol_inv(Sint *n2, double *matrix, Sint *flag2);
void gchol_solve(Sint *n2, double *matrix, double *y, Sint *flag2);
void gchol_bds(Sint   *nb,     Sint   *bs2,  Sint *n2,
	       double *dmat,   double *rmat, double toler[]) ;
void gchol_bdsinv(Sint   *nb,     Sint   *bs2,  Sint *n2,
		  double *dmat,   double *rmat, double *toler,
                  Sint   *flag);
void gchol_bdssolve(Sint   *nb,     Sint   *bs2,  Sint *n2,
		    double *blocks, double *rmat, double *toler,
		    double *y,      Sint   *flag);
SEXP gcback(SEXP sr,   SEXP  sx,    SEXP supper, SEXP sk);
SEXP gcback2(SEXP sblocksize,   SEXP  sblocks,    SEXP srmat,
	      SEXP sx,           SEXP  supper);

