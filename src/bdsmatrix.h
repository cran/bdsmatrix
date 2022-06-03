/*
** This contains the prototype calls for all the .c functions that
**    are called by another C function, or by R
**  It stops errors due to having things declared differently
**    in different routines, and serves as input to R_init_bdsmatrix
*/
void bdsmatrix_index1(int *nblock, int *bsize, int *flag,
		     int *nrow,   int *rows,  int *indexa,
		     int *indexb, int *indexc);

void bdsmatrix_index2(int *nblock, int *bsize, 
		      int *rows,   int *cols);
void bdsmatrix_index3(int *nblock, int *bsize, int *index);

void bdsmatrix_prod(int *nb,	  int *bsize,     int *ydim,
		    double *bmat, double *rmat,    double *offdiag,
		    double *temp, int *itemp,     double *y);

void bdsmatrix_prod2(int nblock,     int *bsize,     int nrow,
		     double *bmat,   double *rmat,  
		     double *y,      double *result, int *itemp);

void bdsmatrix_prod3(int *nr,     int *nb,	  int *bsize,     
		     double *bmat, double *rmat,  int *rhs,
		     int *ny2,	   double *y,     double *temp);

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

void gchol(int *n2, double *matrix, double *toler);
void gchol_inv(int *n2, double *matrix, int *flag2);
void gchol_solve(int *n2, double *matrix, double *y, int *flag2);
void gchol_bds(int   *nb,     int   *bs2,  int *n2,
	       double *dmat,   double *rmat, double toler[]) ;
void gchol_bdsinv(int   *nb,     int   *bs2,  int *n2,
		  double *dmat,   double *rmat, double *toler,
                  int   *flag);
void gchol_bdssolve(int   *nb,     int   *bs2,  int *n2,
		    double *blocks, double *rmat, double *toler,
		    double *y,      int   *flag);
SEXP gcback(SEXP sr,   SEXP  sx,    SEXP supper, SEXP sk);
SEXP gcback2(SEXP sblocksize,   SEXP  sblocks,    SEXP srmat,
	      SEXP sx,           SEXP  supper);

