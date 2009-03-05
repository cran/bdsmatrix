/*
** This routine registers all of the C routines called
**  by other libraries, in particular by coxme, so that
**  they can be found.
*/
#include "bdsS.h"
#include <R_ext/Rdynload.h>
#include "bdsmatrix.h"

void R_init_bdsmatrix(DllInfo *info) {
    R_RegisterCCallable("bdsmatrix","bdsmatrix_prod2", 
			(DL_FUNC) &bdsmatrix_prod2);
    R_RegisterCCallable("bdsmatrix","bdsmatrix_prod4", 
			(DL_FUNC) &bdsmatrix_prod4);
    R_RegisterCCallable("bdsmatrix","cholesky4", (DL_FUNC) &cholesky4);
    R_RegisterCCallable("bdsmatrix","cholesky5", (DL_FUNC) &cholesky5);
    R_RegisterCCallable("bdsmatrix","chinv4", (DL_FUNC) &chinv4);
    R_RegisterCCallable("bdsmatrix","chinv5", (DL_FUNC) &chinv5);
    R_RegisterCCallable("bdsmatrix","chsolve4", (DL_FUNC) &chsolve4);
    R_RegisterCCallable("bdsmatrix","chsolve5", (DL_FUNC) &chsolve5);
    }

/*
** From the R include file R_ext/Rdynload.h
**   typedef void * (*DL_FUNC)();
**
** (DL_FUNC) is defined as "pointer to a function returning a pointer to void", 
**  which is a fib for 3 of the above; but that is what the RegisterCCallable
**  routine wants.  The bsdmatrix.h file contains the truth, which is what
**  users will want.  "Pointer to a function" is the important part.
*/
