#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include "pheno.h"

static const R_CMethodDef CEntries[] = {
    {"Cconnectivity", (DL_FUNC) &Cconnectivity, 5},
    {"Cdate2jul1",    (DL_FUNC) &Cdate2jul1,    3},
    {"Cdate2jul2",    (DL_FUNC) &Cdate2jul2,    4},
    {"Cdaylength",    (DL_FUNC) &Cdaylength,    4},
    {"Cdaysbetween",  (DL_FUNC) &Cdaysbetween,  3},
    {"Cjul2date1",    (DL_FUNC) &Cjul2date1,    3},
    {"Cjul2date2",    (DL_FUNC) &Cjul2date2,    4},
    {"Cleapyear",     (DL_FUNC) &Cleapyear,     2},
    {"Cmaxdaylength", (DL_FUNC) &Cmaxdaylength, 2},
    {NULL, NULL, 0}
};

void R_init_pheno(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

