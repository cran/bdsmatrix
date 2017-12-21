/*
** This routine registers all of the C routines called
**  by other libraries, in particular by coxme, so that
**  they can be found.
*/
#include "bdsS.h"
#include <R_ext/Rdynload.h>
#include <Rversion.h>
#include "bdsmatrix.h"

/* Define routines that can be called from R, or by another package */
static const R_CMethodDef Centries[] = {
    {"Cbdsmatrix_index1",  (DL_FUNC) &bdsmatrix_index1,  8},
    {"Cbdsmatrix_index2",  (DL_FUNC) &bdsmatrix_index2,  4},
    {"Cbdsmatrix_index3",  (DL_FUNC) &bdsmatrix_index3,  3},
    {"Cbdsmatrix_prod",  (DL_FUNC) &bdsmatrix_prod,  9},
    {"Cbdsmatrix_prod2", (DL_FUNC) &bdsmatrix_prod2, 8},
    {"Cbdsmatrix_prod3", (DL_FUNC) &bdsmatrix_prod3, 9},
    {"Cbdsmatrix_prod4", (DL_FUNC) &bdsmatrix_prod4, 7},
    {"Ccholesky4",       (DL_FUNC) &cholesky4,       6},
    {"Ccholesky5",       (DL_FUNC) &cholesky5,       3},
    {"Cchinv4",          (DL_FUNC) &chinv4,          6},
    {"chinv5",           (DL_FUNC) &chinv5,          3},
    {"Cchsolve4",        (DL_FUNC) &chsolve4,        7},
    {"Cchsolve5",        (DL_FUNC) &chsolve5,        4},
    {"Cgchol_bds",       (DL_FUNC) &gchol_bds,       6},
    {"Cgchol",           (DL_FUNC) &gchol,           3},
    {"Cgchol_bdsinv",    (DL_FUNC) &gchol_bdsinv,    7},
    {"Cgchol_bdssolve",  (DL_FUNC) &gchol_bdssolve,  8},
    {"Cgchol_inv",       (DL_FUNC) &gchol_inv,       3},
    {"Cgchol_solve",     (DL_FUNC) &gchol_solve,     4},
    {NULL, NULL, 0}
};

static const R_CallMethodDef Callentries[] = {
    {"Cgcback",    (DL_FUNC) &gcback,  4},
    {"Cgcback2",   (DL_FUNC) &gcback2, 5},
    {NULL, NULL, 0}
};


/* The callable routines can be used by other packages */
void R_init_bdsmatrix(DllInfo *dll) {
    R_RegisterCCallable("bdsmatrix","bdsmatrix_prod2", 
			(DL_FUNC) &bdsmatrix_prod2);
    R_RegisterCCallable("bdsmatrix","bdsmatrix_prod4", 
			(DL_FUNC) &bdsmatrix_prod4);
    R_RegisterCCallable("bdsmatrix","cholesky4", (DL_FUNC) &cholesky4);
    R_RegisterCCallable("bdsmatrix","cholesky5", (DL_FUNC) &cholesky5);
    R_RegisterCCallable("bdsmatrix","chinv4",    (DL_FUNC) &chinv4);
    R_RegisterCCallable("bdsmatrix","chinv5",    (DL_FUNC) &chinv5);
    R_RegisterCCallable("bdsmatrix","chsolve4",  (DL_FUNC) &chsolve4);
    R_RegisterCCallable("bdsmatrix","chsolve5",  (DL_FUNC) &chsolve5);

    /* register the interal routines.  We have no .Fortran or .External
        call in the code, hence the NULL, NULL at the end */
    R_registerRoutines(dll, Centries, Callentries, NULL, NULL);

    /* The following line makes only those routines defined above
       available to outside packages, i.e., internal things like
       dmatrix() are now invisible.
    */
    R_useDynamicSymbols(dll, FALSE); 

    /*
    ** This line makes them only available via the symbols above
    **  i.e., .Call("tmerge", ) won't work but .Call(Ctmerge, )  will
    ** This feature was added in version 3.0.0
    */
#if defined(R_VERSION) && R_VERSION >= R_Version(3, 0, 0)
    R_forceSymbols(dll, TRUE);
#endif
}
