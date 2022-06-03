#include "R.h"
#include "Rinternals.h"

/*
** Memory defined with ALLOC is removed automatically by S.
**  That with "Calloc" I have to remove myself.  Use the
**  latter for objects that need to to persist between calls.
** int was used by Splus, the typedef allowed common code
*/
#ifdef USING_R
#define ALLOC(a,b)  R_alloc(a,b)
#else
#define ALLOC(a,b)  S_alloc(a,b)
#endif

