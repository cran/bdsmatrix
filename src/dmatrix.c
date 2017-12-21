/*
** set up ragged arrays, with #of columns and #of rows
*/
#include "bdsS.h"
#include "bdsmatrix.h"

double **dmatrix(double *array, int ncol, int nrow)
    {
    int i;
    double **pointer;

    pointer = (double **) ALLOC(nrow, sizeof(double *));
    for (i=0; i<nrow; i++) {
	pointer[i] = array;
	array += ncol;
	}
    return(pointer);
    }
