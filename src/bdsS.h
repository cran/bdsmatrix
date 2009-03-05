#include "R.h"
#include "Rinternals.h"
#ifdef USING_R
#define S_EVALUATOR    /* Turn this into a "blank line" in R */
#else
/*
** Splus definitions
*/
typedef long Sint;
#endif

/*
** Memory defined with ALLOC is removed automatically by S.
**  That with "Calloc" I have to remove myself.  Use the
**  latter for objects that need to to persist between calls.
*/
#ifdef USING_R
#define ALLOC(a,b)  R_alloc(a,b)
#else
#define ALLOC(a,b)  S_alloc(a,b)
#endif

