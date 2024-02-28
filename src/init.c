#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .C calls */
extern void count_letters(void *, void *);

static const R_CMethodDef CEntries[] = {
    {"count_letters", (DL_FUNC) &count_letters, 2},
    {NULL, NULL, 0}
};

void R_init_canprot(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
/*
    R_registerRoutines(dll, CEntries, NULL, FortranEntries, NULL);
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
*/
    R_useDynamicSymbols(dll, FALSE);
}

